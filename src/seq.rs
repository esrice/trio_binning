use std::io::{BufRead, BufReader, Lines, Read};

pub struct FastaRecord {
    id: String,
    seq: String,
    entry_string: String,
}

impl FastaRecord {
    pub fn new(entry_string: &String) -> Option<FastaRecord> {
        let mut lines_iter = entry_string.split('\n');
        // TODO replace unwraps with actual error handling
        let first_line = lines_iter.next().unwrap();
        let first_word = first_line.split_whitespace().next().unwrap();
        let id = &first_word[1..];
        let mut seq = String::new();

        for line in lines_iter {
            seq.push_str(line);
        }

        Some(FastaRecord {
            id: id.to_string(),
            seq: seq.to_string(),
            entry_string: entry_string.to_string(),
        })
    }

    pub fn id(&self) -> &str { &self.id }
    pub fn seq(&self) -> &str { &self.seq }
    pub fn to_string(&self) -> &str { &self.entry_string }
}

pub struct FastaReader<T> {
    lines_iter: Lines<BufReader<T>>,
    current_entry: String,
}

impl<T: Read> FastaReader<T> {
    pub fn new(file: T) -> FastaReader<T> {
        FastaReader {
            lines_iter: BufReader::new(file).lines(),
            current_entry: String::new(),
        }
    }
}

impl<T: Read> Iterator for FastaReader<T> {
    type Item = FastaRecord;

    // TODO replace unwraps with actual error handling
    fn next(&mut self) -> Option<FastaRecord> {
        while let Some(result) = self.lines_iter.next() {
            let line = format!("{}\n", result.unwrap());
            if line.starts_with(">") {
                if self.current_entry != "" {
                    let new_entry = self.current_entry.clone();
                    self.current_entry = String::from(line);
                    return Some(FastaRecord::new(&new_entry).unwrap());
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
        let rec = FastaRecord::new(&entry_string).unwrap();
        assert_eq!(rec.id(), "id".to_string());
        assert_eq!(rec.seq(), "ACTGAAAAACGT".to_string());
        assert_eq!(rec.to_string(), entry_string);
    }
}
