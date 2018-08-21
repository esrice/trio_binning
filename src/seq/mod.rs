use std::fmt;
use std::path::Path;
use std::fs::File;
use std::io::{self, Read};
use std::{result,error};

pub mod fasta;
pub mod fastq;

/// wrapper for various errors that can occur here
#[derive(Debug)]
pub enum ReaderError {
    Fasta(fasta::FastaError), // errors coming from the fasta module
    Fastq(fastq::FastqError), // errors coming from the fastq module
    Io(io::Error), // errors coming from the io module
    Other(String), // errors generated in this module
}

impl fmt::Display for ReaderError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            ReaderError::Fasta(e) => write!(f, "{}", e),
            ReaderError::Fastq(e) => write!(f, "{}", e),
            ReaderError::Io(e) => write!(f, "{}", e),
            ReaderError::Other(s) => write!(f, "{}", s),
        }
    }
}

impl error::Error for ReaderError {}

type Result<T> = result::Result<T, ReaderError>;

/// wrapper for records containing an ID and a sequence, so that the user does
/// not need to worry about what kind of file is being read, just that it
/// contains sequences with IDs
pub enum Record {
    Fasta(fasta::Record),
    Fastq(fastq::Record),
}

impl Record {
    pub fn id(&self) -> &str {
        match self {
            Record::Fasta(r) => r.id(),
            Record::Fastq(r) => r.id(),
        }
    }

    pub fn seq(&self) -> &str {
        match self {
            Record::Fasta(r) => r.seq(),
            Record::Fastq(r) => r.seq(),
        }
    }
}

impl fmt::Display for Record {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Record::Fasta(r) => write!(f, "{}", r.to_string()),
            Record::Fastq(r) => write!(f, "{}", r.to_string()),
        }
    }
}

/// Gets the extension from a filename.
///
/// # Errors
/// `ReaderError::Other` if filename has no extension
fn get_extension(filename: &str) -> Result<&str> {
    Path::new(filename).extension().and_then(|e| e.to_str())
        .ok_or(ReaderError::Other("Filename has no extension!".to_owned()))
}

/// wrapper for various file types containing sets of sequences, so that the
/// user can just call next() and get generic Records back without worrying
/// about what kind of file is being parsed
pub enum Reader<T> {
    Fasta(fasta::Reader<T>),
    Fastq(fastq::Reader<T>),
}

impl Reader<File> {
    pub fn from_filename(filename: &str) -> Result<Reader<File>> {
        let extension = get_extension(filename)?;
        let file = File::open(filename).map_err(|e| ReaderError::Io(e))?;

        match extension {
            "fasta" | "fa" => Ok(Reader::Fasta(fasta::Reader::new(file))),
            "fastq" | "fq" => Ok(Reader::Fastq(fastq::Reader::new(file))),
            _ => Err(ReaderError::Other(
                    "Do not recognize filename.".to_owned())),
        }
    }
}

impl<T: Read> Iterator for Reader<T> {
    type Item = Result<Record>;

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            Reader::Fasta(r) => r.next().map(|s| s.map(|t| Record::Fasta(t))
                .map_err(|e| ReaderError::Fasta(e))),
            Reader::Fastq(r) => r.next().map(|s| s.map(|t| Record::Fastq(t))
                .map_err(|e| ReaderError::Fastq(e))),
        }
    }
}

