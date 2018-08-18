extern crate trio_binning;

use trio_binning::seq;
use std::fs::File;

fn main() {
    let reader = seq::FastaReader::new(File::open("test.fa").unwrap());

    for read in reader {
        println!("ID:{}\tSEQ:{}", read.id(), read.seq());
    }
}
