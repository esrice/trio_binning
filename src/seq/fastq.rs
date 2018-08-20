use std::io::{self, BufRead, BufReader, Lines, Read};
use std::error;
use std::fmt;
use std::iter::Enumerate;

#[derive(Debug)]
pub enum FastqError {
    Parse,
    MissingLine,
    DefLine,
    ParseLine(usize),
    Io(io::Error),
}

impl fmt::Display for FastqError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            FastqError::Parse => write!(f, "Parsing error!"),
            FastqError::MissingLine => write!(f, "Missing line!"),
            FastqError::DefLine => write!(f, "Error parsing defline!"),
            FastqError::ParseLine(line) =>
                write!(f, "Problem with entry starting on line {}", line),
            FastqError::Io(e) => write!(f, "{}", e),
        }
    }
}

impl error::Error for FastqError {}

#[derive(Clone, Debug)]
pub struct Record {
    id: String,
    seq: String,
    qual: String,
    entry_string: String,
}

impl Record {
    /// Creates a new Record from a &String containing a fastq entry.
    /// Returns None if the string is empty.
    pub fn new(entry_string: &String) -> Result<Record, FastqError> {
        let lines: Vec<&str> = entry_string.split('\n').collect();

        Ok(Record {
            id: lines.get(0).ok_or(FastqError::MissingLine)
                .and_then(|l| get_id_from_defline(&l))?.to_string(),
            seq: lines.get(1).ok_or(FastqError::MissingLine)
                .map(|l| l.trim().to_string())?,
            qual: lines.get(3).ok_or(FastqError::MissingLine)
                .map(|l| l.trim().to_string())?,
            entry_string: entry_string.to_owned(),
        })
    }

    pub fn id(&self) -> &str { &self.id }
    pub fn seq(&self) -> &str { &self.seq }
    pub fn qual(&self) -> &str { &self.qual }
    pub fn to_string(&self) -> &str { &self.entry_string }
}

/// Takes a fasta defline (e.g., "@seqID sequence desccription") and returns the
/// ID of the entry (e.g., "SeqID")
///
/// # Errors
/// Returns Err("Parsing error!") if an ID cannot be found in the defline, e.g.,
/// if the defline is empty or there is a space after the ">"
fn get_id_from_defline(defline: &str) -> Result<&str, FastqError> {
    defline.split_whitespace().next() // get the first word
        .ok_or(FastqError::DefLine)
        .map(|w| w.trim_left_matches('@')) // trim the '@' delimiter
}

pub struct Reader<T> {
    lines_enum_iter: Enumerate<Lines<BufReader<T>>>,
    line_num: usize,
}

impl<T: Read> Reader<T> {
    /// Creates a new fasta Reader that reads from `file`. `file` can be
    /// anything that implements `std::io::Read`, e.g., `std::io::File`
    pub fn new(file: T) -> Reader<T> {
        Reader {
            lines_enum_iter: BufReader::new(file).lines().enumerate(),
            line_num: 0,
        }
    }
}

impl<T: Read> Iterator for Reader<T> {
    type Item = Result<Record, FastqError>;

    fn next(&mut self) -> Option<Self::Item> {
        // if we're at EOF before starting a new entry, return None
        let (line_num, mut entry_string) = match self.lines_enum_iter.next() {
            Some((i,r)) => match r {
                Ok(l) => (i,l+"\n"),
                Err(e) => return Some(Err(FastqError::Io(e))),
            },
            None => return None,
        };

        self.line_num = line_num + 1; // lines are 1-indexed

        // if we reach EOF before finishing the entry, that's an incomplete
        // entry and indicative of a misformatted fastq file
        for _ in 0..3 {
            let next_line = match self.lines_enum_iter.next() {
                Some((_i,r)) => match r {
                    Ok(l) => l,
                    Err(e) => return Some(Err(FastqError::Io(e))),
                },
                None => return Some(Err(FastqError::ParseLine(self.line_num))),
            };

            entry_string.push_str(&(next_line + "\n"));
        }

        Some(Record::new(&entry_string))
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn fastq_record() {
        let entry_string = "@id\nACTG\n+\nQQQQ".to_string();
        let rec = Record::new(&entry_string).unwrap();
        assert_eq!(rec.id(), "id".to_string());
        assert_eq!(rec.seq(), "ACTG".to_string());
        assert_eq!(rec.qual(), "QQQQ".to_string());
        assert_eq!(rec.to_string(), entry_string);
    }
}
