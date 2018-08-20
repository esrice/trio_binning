use std::io::{self, BufRead, BufReader, Lines, Read};
use std::error;
use std::fmt;

/// a variant of try!/? to use in functions/methods that return
/// Option<Result<T,E>>. For example, if you call a function that returns a
/// Result from inside an iterator's next() method and want to propogate any
/// errors that come from that function, you can't use try! because next() has
/// to return an Option but try! returns a Result.
macro_rules! try_or_some_err {
    ($x:expr) => (match $x {
        Ok(val) => val,
        Err(err) => return Some(Err(err)),
    });
}

#[derive(Debug)]
pub enum FastaError {
    Parse,
    ParseLine(usize),
    Io(io::Error),
}

impl fmt::Display for FastaError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            FastaError::Parse => write!(f, "Parsing error!"),
            FastaError::ParseLine(line) =>
                write!(f, "Parsing error on line {}", line),
            FastaError::Io(e) => write!(f, "{}", e),
        }
    }
}

impl error::Error for FastaError {}

#[derive(Clone, Debug)]
pub struct Record {
    id: String,
    seq: String,
    entry_string: String,
}

impl Record {
    /// Creates a new Record from a &String containing a fasta entry.
    /// Returns None if the string is empty.
    pub fn new(entry_string: &String) -> Result<Record, FastaError> {
        let mut lines_iter = entry_string.split('\n');

        let id = lines_iter.next()
            .ok_or(FastaError::Parse)
            .and_then(|l| get_id_from_defline(&l))?
            .to_string();

        let mut seq = String::new();

        for line in lines_iter {
            seq.push_str(line);
        }

        Ok(Record {
            id: id,
            seq: seq,
            entry_string: entry_string.to_owned(),
        })
    }

    pub fn id(&self) -> &str { &self.id }
    pub fn seq(&self) -> &str { &self.seq }
    pub fn to_string(&self) -> &str { &self.entry_string }
}

/// Takes a fasta defline (e.g., ">seqID sequence desccription") and returns the
/// ID of the entry (e.g., "SeqID")
///
/// # Errors
/// Returns Err("Parsing error!") if an ID cannot be found in the defline, e.g.,
/// if the defline is empty or there is a space after the ">"
fn get_id_from_defline(defline: &str) -> Result<&str, FastaError> {
    defline.split_whitespace().next() // get the first word
        .ok_or(FastaError::Parse)
        .map(|w| w.trim_left_matches('>')) // trim the '>' delimiter
}

pub struct Reader<T> {
    lines_iter: Lines<BufReader<T>>,
    current_entry: Record,
    current_line_number: usize,
}

impl<T: Read> Reader<T> {
    /// Creates a new fasta Reader that reads from `file`. `file` can be
    /// anything that implements `std::io::Read`, e.g., `std::io::File`
    pub fn new(file: T) -> Reader<T> {
        Reader {
            lines_iter: BufReader::new(file).lines(),
            current_entry: Record {
                id: String::new(),
                seq: String::new(),
                entry_string: String::new(),
            },
            current_line_number: 0,
        }
    }
}

impl<T: Read> Iterator for Reader<T> {
    type Item = Result<Record, FastaError>;

    fn next(&mut self) -> Option<Result<Record, FastaError>> {
        while let Some(result) = self.lines_iter.next() {
            self.current_line_number += 1;

            let line = try_or_some_err!(result.map_err(|e| FastaError::Io(e)));

            if line.starts_with(">") {
                if self.current_entry.entry_string != "" {
                    // we have reached the beginning of a new entry, so we move
                    // the instance of Record representing the current one to a
                    // new variable, start a new instance of Record for the new
                    // one, and then return the completed one.
                    let finished_entry = self.current_entry.clone();
                    self.current_entry = Record {
                        id: try_or_some_err!(get_id_from_defline(&line)
                            .map_err(|_| FastaError::ParseLine(
                                self.current_line_number)))
                            .to_string(),
                        seq: String::new(),
                        entry_string: String::from(line),
                    };
                    return Some(Ok(finished_entry));
                } else {
                    // we're on the first line, so don't return anything; just
                    // update the entry string and id.
                    self.current_entry.entry_string.push_str(&line);
                    self.current_entry.id = try_or_some_err!(
                        get_id_from_defline(&line)
                        .map_err(|_| FastaError::ParseLine(
                                self.current_line_number)))
                        .to_string();
                }
            } else { // line is not the defline
                if self.current_entry.id == "" {
                    // must start the file with a defline!
                    return Some(Err(FastaError::ParseLine(
                                self.current_line_number)));
                } else {
                    self.current_entry.entry_string.push_str(&line);
                    self.current_entry.seq.push_str(&line.trim());
                }
            }
        }
       
        // we've reached EOF, so return the final entry, or None if we already
        // did that
        if self.current_entry.entry_string != "" {
            let finished_entry = self.current_entry.clone();

            // change current_entry.entry_string to an empty String so that the
            // next time next() is called, we know to return None
            self.current_entry.entry_string = String::new();

            Some(Ok(finished_entry))
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn fasta_record() {
        let entry_string = ">id\nACTG\nAAAA\nACGT".to_string();
        let rec = Record::new(&entry_string).unwrap();
        assert_eq!(rec.id(), "id".to_string());
        assert_eq!(rec.seq(), "ACTGAAAAACGT".to_string());
        assert_eq!(rec.to_string(), entry_string);
    }
}
