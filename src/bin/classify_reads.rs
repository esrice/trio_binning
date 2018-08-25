extern crate trio_binning;
extern crate clap;

use trio_binning::kmer::*;
use trio_binning::classify::*;
use std::fs::File;
use clap::{Arg, App, ArgGroup, ArgMatches};
use std::{process, error, thread, fmt};

type BoxResult<T> = Result<T, Box<error::Error>>;

#[derive(Debug)]
struct SimpleError {
    message: String,
}

impl SimpleError {
    fn new(message: String) -> Box<SimpleError> {
        Box::new(SimpleError {
            message: message,
        })
    }
}

impl fmt::Display for SimpleError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.message)
    }
}

impl error::Error for SimpleError {}

fn parse_args() -> ArgMatches<'static> {
    App::new("classify_reads")
        .version("0.1.0")
        .author("Edward S. Rice <erice11@unl.edu>")
        .about("Classify reads based on haplotype")
        .arg(Arg::with_name("threads")
             .short("t")
             .long("threads")
             .takes_value(true)
             .default_value("1")
             .help("Number of threads to use"))
        .arg(Arg::with_name("hapA-kmers")
             .short("a")
             .long("hapA-kmers")
             .required(true)
             .takes_value(true)
             .help("File with kmers in canonical representation, \
                   one per line, unique to haplotype A"))
        .arg(Arg::with_name("hapB-kmers")
             .short("b")
             .long("hapB-kmers")
             .required(true)
             .takes_value(true)
             .help("File with kmers in canonical representation, \
                   one per line, unique to haplotype B"))
        .group(ArgGroup::with_name("input-reads")
               .args(&["input-unpaired", "input-paired-end"])
               .required(true))
        .arg(Arg::with_name("input-unpaired")
             .short("u")
             .long("input-unpaired")
             .takes_value(true)
             .help("Fasta/q/bam file containing unpaired reads to classify, \
                   e.g. PacBio"))
        .arg(Arg::with_name("input-paired-end")
             .short("p")
             .long("input-paired-end")
             .takes_value(true)
             .number_of_values(2)
             .help("A pair of fastq files containing paired reads to classify"))
        .arg(Arg::with_name("hapA-out-prefix")
             .short("A")
             .long("hapA-out-prefix")
             .takes_value(true)
             .default_value("hapA")
             .help("Prefix for haplotype A output"))
        .arg(Arg::with_name("hapB-out-prefix")
             .short("B")
             .long("hapB-out-prefix")
             .takes_value(true)
             .default_value("hapB")
             .help("Prefix for haplotype B output"))
        .arg(Arg::with_name("hapU-out-prefix")
             .short("U")
             .long("hapU-out-prefix")
             .takes_value(true)
             .default_value("hapU")
             .help("Prefix for haplotype U output"))
        .arg(Arg::with_name("compress-output")
             .short("c")
             .long("compress-output")
             .help("Output gz-compressed files"))
        .get_matches()
}

//fn asdf(filename: String) -> Result<> {
//    let f = File::open(filename)?;
//    let s = read_kmers_into_set(f)?;
//
//}

fn simple_error<E: 'static + error::Error>(e: E) -> SimpleError {
    SimpleError {
        message: e.to_string(),
    }
}

fn run() -> BoxResult<()> {
    let args = parse_args();

    // get the number of threads to use
    let num_threads = match args.value_of("threads").unwrap().parse::<u32>() {
        Ok(t) => {
            if t >= 1 { t } else {
                return Err(SimpleError::new(
                    format!("Number of threads must be >= 1: {}", t)
                ))
            }
        },
        Err(_) =>
            return Err(SimpleError::new(
                format!("--threads argument not an integer: {}",
                        args.value_of("threads").unwrap())
            )),
    };

    // figure out k by looking at the first line of one of the kmers file
    let k = get_kmer_size(File::open(args.value_of("hapA-kmers").unwrap())?)?;

    // read k-mers into HashSets
    let (hap_a_kmers, hap_b_kmers);
    if num_threads > 1 { // trying out some concurrency!
        let hap_a_kmers_filename = args.value_of("hapA-kmers")
            .unwrap().to_string();

        // read the kmers from haplotype A in a spawned thread
        // error::Error does not implement Send, so the thread has to return a
        // concrete error type, in this case, SimpleError.
        let handle = thread::spawn(move ||
            File::open(hap_a_kmers_filename).map_err(simple_error)
                .and_then(|f| read_kmers_into_set(f).map_err(simple_error)));

        // read the kmers from haplotype B in the main thread
        hap_b_kmers = read_kmers_into_set(File::open(
                args.value_of("hapB-kmers").unwrap())?)?;

        // wait to continue until the spawned thread is done
        hap_a_kmers = handle.join().unwrap()?;
    } else {
        hap_a_kmers = read_kmers_into_set(File::open(
                args.value_of("hapA-kmers").unwrap())?)?;
        hap_b_kmers = read_kmers_into_set(File::open(
                args.value_of("hapB-kmers").unwrap())?)?;
    }

    // call the correct function depending on whether the input is unpaired
    // reads or paired-end reads
    match args.value_of("input-unpaired") {
        Some(input_reads_filename) => {
            classify_unpaired(&hap_a_kmers, &hap_b_kmers, input_reads_filename,
                              args.value_of("hapA-out-prefix").unwrap(),
                              args.value_of("hapB-out-prefix").unwrap(),
                              args.value_of("hapU-out-prefix").unwrap(),
                              args.is_present("compress-output"), k)?;
        }
        None => {
            let filenames: Vec<&str> = args.values_of("input-paired")
                .unwrap().collect();
            classify_paired(&hap_a_kmers, &hap_b_kmers, filenames[0],
                            filenames[1],
                            args.value_of("hapA-out-prefix").unwrap(),
                            args.value_of("hapB-out-prefix").unwrap(),
                            args.value_of("hapU-out-prefix").unwrap())?;
        }
    }
    Ok(())
}

fn main() {
    if let Err(e) = run() {
        println!("fatal error: {}", e);
        process::exit(1);
    }
}
