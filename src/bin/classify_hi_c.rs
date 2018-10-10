extern crate trio_binning;
extern crate clap;
extern crate rust_htslib;
extern crate textwrap;
extern crate ansi_term;

use clap::{Arg, App, ArgMatches};
use std::{process, error, fmt};
use ansi_term::Colour;
use self::rust_htslib::bam;
use self::rust_htslib::prelude::*;

type BoxResult<T> = Result<T, Box<error::Error>>;

#[derive(Debug)]
struct TagError {
    message: String,
}

impl TagError {
    fn new() -> Box<TagError> {
        Box::new(TagError {
            message: "Cannot find alignment score (AS) tag in \
                alignment. Try using bwa mem or another aligner \
                that outputs these tags.".to_string(),
        })
    }
}

impl fmt::Display for TagError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", textwrap::fill(&self.message, textwrap::termwidth()))
    }
}

impl error::Error for TagError {}

#[derive(Debug)]
struct SortError {
    message: String,
}

impl SortError {
    fn new() -> Box<SortError> {
        Box::new(SortError {
            message: "The two input files do not have primary alignments \
                in the same order. The bam files must contain the same \
                reads in the same order. Make sure the bams are sorted by \
                read name and that you didn't do anything like remove \
                duplicates that could cause the two bams to have different \
                reads. If you did all this but are still getting this error, \
                please file a bug report at
                https://github.com/esrice/trio_binning/issues".to_string(),
        })
    }
}

impl fmt::Display for SortError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", textwrap::fill(&self.message, textwrap::termwidth()))
    }
}

impl error::Error for SortError {}

fn parse_args() -> ArgMatches<'static> {
    App::new("classify_hi_c")
        .version("0.1.0")
        .author("Edward S. Rice <erice11@unl.edu>")
        .about("Classify Hi-C reads based on alignments")
        .arg(Arg::with_name("hapA-in")
             .short("a")
             .long("hapA-in")
             .required(true)
             .takes_value(true)
             .help("Input bam file for haplotype A, sorted by read name"))
        .arg(Arg::with_name("hapB-in")
             .short("b")
             .long("hapB-in")
             .required(true)
             .takes_value(true)
             .help("Input bam file for haplotype B, sorted by read name"))
        .arg(Arg::with_name("hapA-out")
             .short("A")
             .long("hapA-out")
             .required(true)
             .takes_value(true)
             .help("Place to put haplotype A output bam"))
        .arg(Arg::with_name("hapB-out")
             .short("B")
             .long("hapB-out")
             .required(true)
             .takes_value(true)
             .help("Place to put haplotype B output bam"))
        .get_matches()
}

/// The only way to know if we've reached EOF using the `read` function (which
/// I'm using because it's way more efficient than using an iterator) is by
/// whether it returns a ``ReadError`. However, some `ReadError` types indicate
/// actual problems, as opposed to EOF, so this is a wrapper function that
/// returns:
/// - `Ok(true)` if the next record was successfully read into `record`
/// - `Ok(false)` if we've reached EOF
/// - `Err(ReadError)` if there was a problem reading the bam file
fn read_or_eof(reader: &mut bam::Reader,
               record: &mut bam::record::Record) -> BoxResult<bool> {
    match reader.read(record) {
        Ok(_) => return Ok(true),
        Err(e) => {
            if let bam::ReadError::NoMoreRecord = e {
                return Ok(false);
            } else {
                return Err(Box::new(e));
            }
        }
    }
}

/// Keep reading until either the next primary alignment or EOF.
/// Returns:
/// - `Ok(true)` if `record` now contains the next primary alignment
/// - `Ok(false)` if we reached EOF before finding a primary alignment
/// - `Err(e)` if there was a problem reading the bam
fn next_primary(reader: &mut bam::Reader,
                record: &mut bam::record::Record)
        -> BoxResult<bool> {
    loop {
        if !read_or_eof(reader, record)? {
            return Ok(false);
        }
        if !(record.is_secondary() || record.is_supplementary()) {
            return Ok(true);
        }
    }
}

/// Compare the alignment scores of each read in two bam files and output all
/// alignments with score greater than or equal to the alignment of the same
/// read in the opposite bam file. Return the number of reads output to only
/// the first haplotype, only the second haplotype, and both haplotypes.
fn classify_hi_c(in_bam_a: &mut bam::Reader,
                 in_bam_b: &mut bam::Reader,
                 out_bam_a: &mut bam::Writer,
                 out_bam_b: &mut bam::Writer) -> BoxResult<(u64, u64, u64)> {

    // allocate new empty bam records to store actual records
    let mut current_record_a = bam::Record::new();
    let mut current_record_b = bam::Record::new();

    let mut score_a: i64;
    let mut score_b: i64;

    // counters to keep track of the number of reads classified to each bin
    let mut count_a: u64 = 0;
    let mut count_b: u64 = 0;
    let mut count_unclassified: u64 = 0;

    // continue reading records until we reach EOF in one of the files
    let mut eof_a = !next_primary(in_bam_a, &mut current_record_a)?;
    let mut eof_b = !next_primary(in_bam_b, &mut current_record_b)?;
    while !eof_a && !eof_b {
        // check if bam is sorted wrong, and exit with an error if so.
        if current_record_a.qname() != current_record_b.qname() {
            return Err(SortError::new());
        }

        // now that we have alignments of the same read in current_record_a
        // and current_record_b, we can compare them. First, get the
        // alignment scores from the "AS" tag of the bam record:
        score_a = current_record_a.aux(b"AS")
            .ok_or(TagError::new())
            .map(|s| s.integer())?;
        score_b = current_record_b.aux(b"AS")
            .ok_or(TagError::new())
            .map(|s| s.integer())?;

        // then, output the higher-scoring alignment to its corresponding
        // output file, or both if the scores are equal.
        if score_a > score_b {
            out_bam_a.write(&current_record_a)?;
            count_a += 1;
        } else if score_b > score_a {
            out_bam_b.write(&current_record_b)?;
            count_b += 1;
        } else { // tie, so write both records
            out_bam_a.write(&current_record_a)?;
            out_bam_b.write(&current_record_b)?;
            count_unclassified += 1;
        }

        eof_a = !next_primary(in_bam_a, &mut current_record_a)?;
        eof_b = !next_primary(in_bam_b, &mut current_record_b)?;
    }

    Ok((count_a, count_b, count_unclassified))
}

fn run() -> BoxResult<()> {
    let args = parse_args();

    // open bam readers
    let mut in_bam_a = bam::Reader::from_path(
        args.value_of("hapA-in").unwrap())?;
    let mut in_bam_b = bam::Reader::from_path(
        args.value_of("hapB-in").unwrap())?;

    // get headers from input files so we can copy them to output files
    let header_a = bam::Header::from_template(in_bam_a.header());
    let header_b = bam::Header::from_template(in_bam_b.header());

    // open bam writers
    let mut out_bam_a = bam::Writer::from_path(
        args.value_of("hapA-out").unwrap(), &header_a)?;
    let mut out_bam_b = bam::Writer::from_path(
        args.value_of("hapB-out").unwrap(), &header_b)?;

    let (count_a, count_b, count_unclassified) = classify_hi_c(
        &mut in_bam_a, &mut in_bam_b, &mut out_bam_a, &mut out_bam_b)?;
    eprintln!("{} {}",
              Colour::Green.paint("# reads classified to haplotype A:"),
              count_a);
    eprintln!("{} {}",
              Colour::Green.paint("# reads classified to haplotype B:"),
              count_b);
    eprintln!("{} {}", Colour::Green.paint("# reads with equal alignment \
                                           scores in both haplotypes:"),
              count_unclassified);

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
