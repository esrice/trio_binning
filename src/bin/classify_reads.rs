extern crate trio_binning;
extern crate clap;

use trio_binning::kmer::*;
use trio_binning::file::*;
use clap::{Arg, App};
use std::process;
use std::error::Error;

fn run() -> Result<(), Box<Error>> {
    let matches = App::new("classify_reads")
        .version("0.1.0")
        .author("Edward S. Rice <erice11@unl.edu>")
        .about("Classify reads based on haplotype")
        .arg(Arg::with_name("hapA-kmers")
             .short("a")
             .long("hapA-kmers")
             .required(true)
             .takes_value(true)
             .help("File with kmers, one per line, unique to haplotype A"))
        .arg(Arg::with_name("hapB-kmers")
             .short("b")
             .long("hapB-kmers")
             .required(true)
             .takes_value(true)
             .help("File with kmers, one per line, unique to haplotype B"))
        .get_matches();

    // unwrapping is safe because these arguments are required
    let hap_a_kmers_filename = matches.value_of("hapA-kmers").unwrap();
    let hap_b_kmers_filename = matches.value_of("hapB-kmers").unwrap();

    // figure out k by looking at the first line of one of the kmers file
    let k = get_kmer_size(try_open(hap_a_kmers_filename)?)?;

    let hap_a_kmers = read_kmers_into_set(try_open(hap_a_kmers_filename)?)?;
    let hap_b_kmers = read_kmers_into_set(try_open(hap_b_kmers_filename)?)?;

    for kmer_bits in &hap_a_kmers {
        println!("Hap A: {}", bits_to_kmer(*kmer_bits, k)?);
    }

    for kmer_bits in &hap_b_kmers {
        println!("Hap B: {}", bits_to_kmer(*kmer_bits, k)?);
    }

    Ok(())
}

fn main() {
    if let Err(e) = run() {
        println!("fatal error: {}", e);
        process::exit(1);
    }
}
