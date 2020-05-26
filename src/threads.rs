use std::{thread, error, result};
use std::sync::{Arc, Mutex, mpsc};
use std::io::{Read, Write};
use classify::{count_kmers_in_read, calc_scaling_factors};
use seq;
use kmer;

// number of Records that can be kept in the channel buffer at once
const BUFFER_SIZE: usize = 100;

struct CountResult {
    hap_a_count: u32,
    hap_b_count: u32,
    record: seq::SeqRecord,
}

/// message to be passed to the counter pool
enum RecordMessage {
    Record(seq::SeqRecord), // "here's a new record to analyze"
    Done, // "there are no more records to analyze, so wrap it up now"
}

enum ResultMessage {
    Count(CountResult), // "here's a new result for you to print out"
    Done, // "there is nothing more to print out, so wrap it up"
}

/// Create a new Reader that has a SeqReader to read records from, a sender
/// to send records to a pool, and a thread to do the work
fn spawn_reader<T>(seq_reader: seq::SeqReader<T>,
                   record_sender: mpsc::SyncSender<RecordMessage>)
    -> thread::JoinHandle<mpsc::SyncSender<RecordMessage>>
    where T: Read + Send + 'static {
    thread::spawn(move || {
        for result in seq_reader {
            // TODO proper error handling
            record_sender.send(RecordMessage::Record(result.unwrap())).unwrap();
        }
        return record_sender; // give ownership back
    })
}

fn spawn_counter(hap_a_kmers: Arc<kmer::KmerSet>,
                 hap_b_kmers: Arc<kmer::KmerSet>,
                 record_receiver: Arc<Mutex<mpsc::Receiver<RecordMessage>>>,
                 result_sender: Arc<Mutex<mpsc::Sender<ResultMessage>>>,
                 k: usize)
    -> thread::JoinHandle<()> {
    thread::spawn(move || {
        loop {
            // TODO proper error handling
            match record_receiver.lock().unwrap().recv().unwrap() {
                RecordMessage::Record(record) => {
                    let (hap_a_count, hap_b_count) = count_kmers_in_read(
                            &hap_a_kmers, &hap_b_kmers, &record, k).unwrap();

                    result_sender.lock().unwrap().send(ResultMessage::Count(
                        CountResult {
                            hap_a_count,
                            hap_b_count,
                            record,
                        }
                    )).unwrap();
                },
                RecordMessage::Done => break,
            }
        }
    })
}

fn spawn_writer<T>(mut hap_a_out: T, mut hap_b_out: T, mut hap_u_out: T,
                   result_receiver: mpsc::Receiver<ResultMessage>,
                   scaling_factor_a: f32, scaling_factor_b: f32)
    -> thread::JoinHandle<()>
    where T: Write + Send + 'static {

    thread::spawn(move || {
        loop {
            match result_receiver.recv().unwrap() {
                ResultMessage::Count(results) => {
                    let hap_a_score = (results.hap_a_count as f32)
                        * scaling_factor_a;
                    let hap_b_score = (results.hap_b_count as f32)
                        * scaling_factor_b;

                    let haplotype;
                    if hap_a_score > hap_b_score {
                        hap_a_out.write(results.record.entry_string
                                        .as_bytes()).unwrap();
                        haplotype = "A";
                    } else if hap_b_score > hap_a_score {
                        hap_b_out.write(results.record.entry_string
                                        .as_bytes()).unwrap();
                        haplotype = "B";
                    } else {
                        hap_u_out.write(results.record.entry_string
                                        .as_bytes()).unwrap();
                        haplotype = "U";
                    }
                    println!("{}\t{}\t{}\t{}", results.record.id, haplotype,
                             hap_a_score, hap_b_score);
                },
                ResultMessage::Done => break,
            }
        }
    })
}

type Result<T> = result::Result<T, Box<dyn error::Error>>;
pub fn count_kmers_multithreaded<T,U>(hap_a_kmers: kmer::KmerSet,
                                      hap_b_kmers: kmer::KmerSet,
                                      reads_in: seq::SeqReader<T>,
                                      hap_a_out: U,
                                      hap_b_out: U,
                                      hap_u_out: U,
                                      num_threads: usize,
                                      k: usize) -> Result<()>
    where T: Read + Send + 'static,
          U: Write + Send + 'static {

    // calculate scaling factors
    let (scaling_factor_a, scaling_factor_b) =
        calc_scaling_factors(&hap_a_kmers, &hap_b_kmers);

    // open up two channels: one for sending SeqRecords from the Reader to
    // the thread pool and one for sending results from the thread pool to
    // the writer. The first channel is a synchronous channel because the
    // reader thread will likely be able to read records faster than the
    // counter pool can count them, so we don't want to the channel to get
    // full of records, taking up lots of memory.
    let (record_sender, record_receiver) = mpsc::sync_channel(BUFFER_SIZE);
    let (result_sender, result_receiver) = mpsc::channel();

    // set up a reader thread that can send records through the channel we
    // just opened
    let reader = spawn_reader(reads_in, record_sender);

    // all the Counters must have access to record_receiver and
    // result_sender, so we have to wrap these in atomic smart pointers
    let record_receiver = Arc::new(Mutex::new(record_receiver));
    let result_sender = Arc::new(Mutex::new(result_sender));

    // we also need atomic smart pointers for the kmer sets
    let hap_a_kmers = Arc::new(hap_a_kmers);
    let hap_b_kmers = Arc::new(hap_b_kmers);

    // spawn num_threads-2 threads to count kmers
    let mut counters: Vec<thread::JoinHandle<()>> =
        Vec::with_capacity(num_threads);

    for _ in 0..(num_threads-2) {
        counters.push(spawn_counter(Arc::clone(&hap_a_kmers),
                                    Arc::clone(&hap_b_kmers),
                                    Arc::clone(&record_receiver),
                                    Arc::clone(&result_sender), k));
    }

    // spawn a writer thread to write the results of the counting
    let writer = spawn_writer(hap_a_out, hap_b_out, hap_u_out, result_receiver,
                              scaling_factor_a, scaling_factor_b);

    // wait for the reader to finish, which means EOF has been reached, and get
    // ownership of the record sender back
    let record_sender = reader.join().unwrap();

    // send a bunch of "Done" messages down the channel, one for each counter
    for _ in 0..counters.len() {
        record_sender.send(RecordMessage::Done).unwrap();
    }

    // wait for all the counters to finish up
    for counter in counters {
        counter.join().unwrap();
    }

    // send a "Done" message to the writer and wait for it to finish up
    result_sender.lock().unwrap().send(ResultMessage::Done).unwrap();
    writer.join().unwrap();

    Ok(())
}
