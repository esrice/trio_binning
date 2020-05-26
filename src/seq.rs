extern crate seq_io;
extern crate flate2;
extern crate rust_htslib;

use std::io::Read;
use std::fs::File;
use std::{error, result, fmt};
use self::seq_io::{fasta, fastq};
use self::seq_io::fasta::Record as FastaRecord;
use self::seq_io::fastq::Record as FastqRecord;
use self::flate2::read::GzDecoder;
use self::rust_htslib::bam::{self, Read as BamRead, ReadError};

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

type Result<T> = result::Result<T, Box<dyn error::Error>>;

pub struct SeqRecord {
    pub id: String,
    pub seq: String,
    pub entry_string: String,
}

pub enum SeqReader<T: Read + Send + 'static> {
    Fasta(fasta::Reader<T>),
    Fastq(fastq::Reader<T>),
    Bam(bam::Reader, bam::record::Record),
}

fn open_gz_or_uncompressed(filename: &str) ->
    Result<Box<dyn Read + Send + 'static>> {
    if filename.ends_with(".gz") {
        Ok(Box::new(GzDecoder::new(File::open(filename)?)))
    } else {
        Ok(Box::new(File::open(filename)?))
    }
}

impl SeqReader<Box<dyn Read + Send + 'static>> {
    pub fn from_path(path: &str)
        -> Result<SeqReader<Box<dyn Read + Send + 'static>>> {
        if path.ends_with(".bam") {
            return Ok(SeqReader::Bam(
                    bam::Reader::from_path(path)?,
                    bam::record::Record::new()))
        }

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

impl<T: Read + Send + 'static> Iterator for SeqReader<T> {
    type Item = Result<SeqRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        match *self {
            SeqReader::Fasta(ref mut r) => {
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

            SeqReader::Fastq(ref mut r) => {
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

            SeqReader::Bam(ref mut reader, ref mut record) => {
                // can return bam::ReadError
                if let Err(e) = reader.read(record) {
                    return match e {
                        ReadError::NoMoreRecord => None,
                        _ => Some(Err(Box::new(e))),
                    }
                }

                let id = String::from_utf8(record.qname().to_vec()).unwrap();
                let seq = String::from_utf8(record.seq().as_bytes()).unwrap();
                let qual = String::from_utf8(
                    record.qual().iter().map(|q| q + 33).collect()).unwrap();

                // write bam entries as fastq because there's really no good
                // reason to output unaligned reads as bam.
                let entry_string = format!("@{}\n{}\n+\n{}\n",
                                           &id, &seq, &qual);

                Some(Ok(SeqRecord {
                    id: id,
                    seq: seq,
                    entry_string: entry_string,
                }))
            },
        }
    }
}
