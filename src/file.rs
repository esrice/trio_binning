use std::fs::File;
use std::error::Error;
use std::fmt;
use std::io;

#[derive(Debug)]
pub enum IOError {
    ReadError(String, io::Error),
}

impl fmt::Display for IOError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            IOError::ReadError(filename, err) =>
                write!(f, "{}: {}", filename, err),
        }
    }
}

impl Error for IOError {}

pub fn try_open(filename: &str) -> Result<File, IOError> {
    File::open(filename).map_err(|err| IOError::ReadError(
            String::from(filename), err))
}
