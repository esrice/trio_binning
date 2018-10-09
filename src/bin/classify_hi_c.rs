extern crate trio_binning;
extern crate clap;
extern crate rust_htslib;

use clap::{Arg, App, ArgMatches};
use std::{process, error, fmt};
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
        write!(f, "{}", self.message)
    }
}

impl error::Error for TagError {}

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
             .default_value("hapA")
             .help("Place to put haplotype A output bam"))
        .arg(Arg::with_name("hapB-out")
             .short("B")
             .long("hapB-out")
             .required(true)
             .takes_value(true)
             .default_value("hapB")
             .help("Place to put haplotype B output bam"))
        .get_matches()
}

/// If the read ID's of current_record_a and current_record_b are not the same,
/// advance the reader that is behind the other one until they are the same.
/// This function can return three different things:
/// 1. Ok(false) if everything went fine and neither reader is at EOF
/// 2. Ok(true) if one of the readers is at EOF
/// 3. Err if there's an error reading or writing the bam files
fn advance_until(in_bam_a: &mut bam::Reader,
                 in_bam_b: &mut bam::Reader,
                 out_bam_a: &mut bam::Writer,
                 out_bam_b: &mut bam::Writer,
                 current_record_a: &mut bam::Record,
                 current_record_b: &mut bam::Record) -> BoxResult<bool> {

    while current_record_a.qname() > current_record_b.qname() {
        // if our iteration of A is ahead of B, advance B until it catches
        // up, outputting the unpaired records
        out_bam_b.write(current_record_b)?;
        if read_or_eof(in_bam_b, current_record_b)? {
            // make current_record_b empty so we don't print it twice
            current_record_b.set_qname(b"EOF");
            return Ok(true)
        }
    }

    while current_record_b.qname() > current_record_a.qname() {
        // if our iteration of B is ahead of A, advance A until it catches
        // up, outputting the unpaired records
        out_bam_a.write(&current_record_a)?;
        if read_or_eof(in_bam_a, current_record_a)? {
            // set qname of current record to "EOF" so we know we are at the
            // end of the file ugh ugh ugh this is so ugly sorry
            current_record_a.set_qname(b"EOF");
            return Ok(true)
        }
    }

    return Ok(false)
}

/// Read the next record in `in_bam` into `record`. Three things can happen:
/// 1. Next record is successfully read. Return Ok(false).
/// 2. Can't read next record because EOF! Return Ok(true).
/// 3. Can't read next record because of some other problem. Return the error.
fn read_or_eof(in_bam: &mut bam::Reader,
               record: &mut bam::Record) -> BoxResult<bool> {
    match in_bam.read(record) {
        Ok(_) => return Ok(false),
        Err(e) => {
            if let bam::ReadError::NoMoreRecord = e {
                return Ok(true)
            } else {
                return Err(Box::new(e))
            }
        }
    }
}

fn classify_hi_c(in_bam_a: &mut bam::Reader,
                 in_bam_b: &mut bam::Reader,
                 out_bam_a: &mut bam::Writer,
                 out_bam_b: &mut bam::Writer) -> BoxResult<()> {

    // allocate new empty bam records to store actual records
    let mut current_record_a = bam::Record::new();
    let mut current_record_b = bam::Record::new();

    let mut score_a: i64;
    let mut score_b: i64;

    // TODO have some check that the file is sorted by read name, or else the
    // program may exit without error but output empty files, which would be bad

    let mut eof_a = read_or_eof(in_bam_a, &mut current_record_a)?;
    let mut eof_b = read_or_eof(in_bam_b, &mut current_record_b)?;

    // continue reading records until we reach EOF in one of the files
    while !eof_a && !eof_b {
        // we can only compare two alignments if they are of the same read, so
        // if we are not looking at records describing the same read, we need to
        // fix that. `advance_until` returns true if we've reached EOF of one of
        // the bam files, so only do the score comparing stuff if it returns
        // false.
        if !advance_until(in_bam_a, in_bam_b,
                          out_bam_a, out_bam_b,
                          &mut current_record_a, &mut current_record_b)? {

            // now that we have alignments of the same read in current_record_a
            // and current_record_b, we can compare them. First, get the
            // alignment scores from the "AS" tag of the bam record:
            score_a = current_record_a.aux(b"AS")
                .ok_or(TagError::new())
                .map(|s| s.integer())?;
            score_b = current_record_a.aux(b"AS")
                .ok_or(TagError::new())
                .map(|s| s.integer())?;

            // then, output the higher-scoring alignment to its corresponding
            // output file, or both if the scores are equal.
            if score_a >= score_b {
                out_bam_a.write(&current_record_a)?;
            }

            if score_b >= score_a {
                out_bam_b.write(&current_record_b)?;
            }

            eof_a = read_or_eof(in_bam_a, &mut current_record_a)?;
            eof_b = read_or_eof(in_bam_b, &mut current_record_b)?;
        } else { // advance_until reached EOF for one of the files, but which?
            if current_record_a.qname() == b"EOF" {
                eof_a = true;
            }

            if current_record_b.qname() == b"EOF" {
                eof_b = true;
            }
        }
    }

    // now that one of the files has reached EOF, we make sure both have
    while !eof_a {
        out_bam_a.write(&current_record_a)?;
        eof_a = read_or_eof(in_bam_a, &mut current_record_a)?;
    }

    while !eof_b {
        out_bam_b.write(&current_record_b)?;
        eof_b = read_or_eof(in_bam_b, &mut current_record_b)?;
    }

    Ok(())
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

    classify_hi_c(&mut in_bam_a, &mut in_bam_b, &mut out_bam_a, &mut out_bam_b)
}

fn main() {
    if let Err(e) = run() {
        println!("fatal error: {}", e);
        process::exit(1);
    }
}
