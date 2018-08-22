extern crate seq_io;
extern crate flate2;

use std::io::Read;
use std::fs::File;
use std::{error, result, fmt};
use self::seq_io::{fasta, fastq};
use self::seq_io::fasta::Record as FastaRecord;
use self::seq_io::fastq::Record as FastqRecord;
use self::flate2::read::GzDecoder;

#[derive(Debug)]
pub struct ExtensionError {
    path: String,
}

impl fmt::Display for ExtensionError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}: cannot determine file type from extension", self.path)
    }
}

impl error::Error for ExtensionError {}

type Result<T> = result::Result<T, Box<error::Error>>;

pub struct SeqRecord {
    pub id: String,
    pub seq: String,
    pub entry_string: String,
}

pub enum SeqReader<T: Read> {
    Fasta(fasta::Reader<T>),
    Fastq(fastq::Reader<T>),
}

fn open_gz_or_uncompressed(filename: &str) -> Result<Box<dyn Read>> {
    if filename.ends_with(".gz") {
        Ok(Box::new(GzDecoder::new(File::open(filename)?)))
    } else {
        Ok(Box::new(File::open(filename)?))
    }
}

impl SeqReader<Box<dyn Read>> {
    pub fn from_path(path: &str) -> Result<SeqReader<Box<dyn Read>>> {
        let read = open_gz_or_uncompressed(path)?;
        let subpath = if path.ends_with(".gz") { &path[..path.len()-3] }
                      else { path };

        if subpath.ends_with(".fasta") || subpath.ends_with(".fa") {
            Ok(SeqReader::Fasta(fasta::Reader::new(read)))
        } else if subpath.ends_with(".fastq") || subpath.ends_with(".fq") {
            Ok(SeqReader::Fastq(fastq::Reader::new(read)))
        } else {
            Err(Box::new(ExtensionError { path: path.to_string() }))
        }
    }
}

impl<T: Read> Iterator for SeqReader<T> {
    type Item = Result<SeqRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            SeqReader::Fasta(r) => {
                let record = match r.next() {
                    Some(s) => match s {
                        Ok(o) => o,
                        Err(e) => return Some(Err(Box::new(e))),
                    },
                    None => return None,
                };

                let mut entry_utf: Vec<u8> = Vec::new();
                record.write(&mut entry_utf).unwrap();

                Some(Ok(SeqRecord {
                    // I'm not handling utf8 errors because bleh. Maybe later.
                    id: record.id().unwrap().to_string(),
                    seq: String::from_utf8(record.seq().to_vec()).unwrap(),
                    entry_string: String::from_utf8(entry_utf).unwrap(),
                }))
            },

            SeqReader::Fastq(r) => {
                let record = match r.next() {
                    Some(s) => match s {
                        Ok(o) => o,
                        Err(e) => return Some(Err(Box::new(e))),
                    },
                    None => return None,
                };

                let mut entry_utf: Vec<u8> = Vec::new();
                record.write(&mut entry_utf).unwrap();

                Some(Ok(SeqRecord {
                    // I'm not handling utf8 errors because bleh. Maybe later.
                    id: record.id().unwrap().to_string(),
                    seq: String::from_utf8(record.seq().to_vec()).unwrap(),
                    entry_string: String::from_utf8(entry_utf).unwrap(),
                }))
            },
        }
    }
}
