use std::io::{BufRead, BufReader, Lines, Read};
use std::error::Error;

pub struct Record {
    id: String,
    seq: String,
    entry_string: String,
}

type BoxResult<T> = Result<T, Box<Error>>;

impl Record {
    /// Creates a new Record from a &String containing a fasta entry.
    /// Returns None if the string is empty.
    pub fn new(entry_string: &String) -> BoxResult<Record> {
        let mut lines_iter = entry_string.split('\n');

        let first_line = lines_iter.next().ok_or("Parsing error!")?;
        let first_word = first_line.split_whitespace().next()
            .ok_or("Parsing error!")?;
        let id = &first_word[1..];
        let mut seq = String::new();

        for line in lines_iter {
            seq.push_str(line);
        }

        Ok(Record {
            id: id.to_string(),
            seq: seq.to_string(),
            entry_string: entry_string.to_string(),
        })
    }

    pub fn id(&self) -> &str { &self.id }
    pub fn seq(&self) -> &str { &self.seq }
    pub fn to_string(&self) -> &str { &self.entry_string }
}

pub struct Reader<T> {
    lines_iter: Lines<BufReader<T>>,
    current_entry: String,
}

impl<T: Read> Reader<T> {
    pub fn new(file: T) -> Reader<T> {
        Reader {
            lines_iter: BufReader::new(file).lines(),
            current_entry: String::new(),
        }
    }
}

impl<T: Read> Iterator for Reader<T> {
    type Item = BoxResult<Record>;

    fn next(&mut self) -> Option<BoxResult<Record>> {
        while let Some(result) = self.lines_iter.next() {
            let line = match result {
                Ok(r) => r,
                Err(e) => return Some(Err(Box::new(e))),
            };
            if line.starts_with(">") {
                if self.current_entry != "" {
                    let new_entry = self.current_entry.clone();
                    self.current_entry = String::from(line);
                    return Some(Record::new(&new_entry));
                } else {
                    self.current_entry.push_str(&line);
                }
            } else {
                self.current_entry.push_str(&line);
            }
        }
        None
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
