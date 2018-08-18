use std::error::Error;

/// Wrapper for reads from various file types, so that the user does not need
/// to worry about what the file type is.
pub trait SeqRecord {
    fn id(&self) -> &str;
    fn seq(&self) -> &str;
    fn to_string(&self) -> &str;
}

pub struct FastaRecord {
    id: String,
    seq: String,
    entry_string: String,
}

type BoxResult<T> = Result<T, Box<Error>>;

impl FastaRecord {
    pub fn new(entry_string: &String) -> BoxResult<FastaRecord> {
        let mut lines_iter = entry_string.split('\n');
        let first_line = lines_iter.next().unwrap();
        let first_word = first_line.split_whitespace().next().unwrap();
        let id = &first_word[1..];
        let mut seq = String::new();

        for line in lines_iter {
            seq.push_str(line);
        }

        Ok(FastaRecord {
            id: id.to_string(),
            seq: seq.to_string(),
            entry_string: entry_string.to_string(),
        })
    }
}

impl SeqRecord for FastaRecord {
    fn id(&self) -> &str { &self.id }
    fn seq(&self) -> &str { &self.seq }
    fn to_string(&self) -> &str { &self.entry_string }
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
