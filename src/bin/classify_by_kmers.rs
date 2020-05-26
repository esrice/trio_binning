extern crate trio_binning;
extern crate clap;
extern crate ansi_term;

use trio_binning::kmer::*;
use trio_binning::classify::*;
use std::fs::File;
use clap::{Arg, App, ArgMatches};
use std::{process, error, thread, fmt};
use ansi_term::Colour;

type BoxResult<T> = Result<T, Box<dyn error::Error>>;

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
    App::new("classify_by_kmers")
        .version("0.2.0")
        .author("Edward S. Rice <erice11@unl.edu>")
        .about("Classify reads based on haplotype using k-mers")
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
        .arg(Arg::with_name("input-reads")
             .short("i")
             .long("input-reads")
             .required(true)
             .takes_value(true)
             .help("Fasta/q/bam file containing unpaired reads to classify, \
                   e.g. PacBio"))
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
        .get_matches()
}

fn simple_error<E: 'static + error::Error>(e: E) -> SimpleError {
    SimpleError {
        message: e.to_string(),
    }
}

fn run() -> BoxResult<()> {
    let args = parse_args();

    // get the number of threads to use
    let num_threads = match args.value_of("threads").unwrap().parse::<usize>() {
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
    eprintln!("{}", Colour::Blue.bold().paint("Reading k-mers into sets..."));
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
    eprintln!("{}", Colour::Blue.bold().paint("Classifying reads..."));
    classify_unpaired(hap_a_kmers, hap_b_kmers,
                      args.value_of("input-reads").unwrap(),
                      args.value_of("hapA-out-prefix").unwrap(),
                      args.value_of("hapB-out-prefix").unwrap(),
                      args.value_of("hapU-out-prefix").unwrap(),
                      false, k, num_threads)?;

    Ok(())
}

fn main() {
    if let Err(e) = run() {
        println!(
            "{}\n{}",
            Colour::Red.bold().paint("--------FATAL ERROR:---------"),
            Colour::Red.paint(e.to_string()));
        process::exit(1);
    }
}
