extern crate trio_binning;

use trio_binning::seq::{self, fastq};

fn main() {
    let reader = seq::Reader::from_filename("test.fa").unwrap();

    for result in reader {
        let read = result.unwrap();
        println!("ID:{}\tSEQ:{}", read.id(), read.seq());
    }

    match fasta::Record::new(&String::from("")) {
        Ok(_) => panic!("Should not get Ok for an empty record!"),
        Err(e) => eprintln!("ERROR: {}", e),
    }

    let mut bad_reader = seq::Reader::from_filename("test_bad.fastq").unwrap();
    let good_read = bad_reader.next().unwrap().unwrap();
    println!("ID:{}\tSEQ:{}", good_read.id(), good_read.seq());
    match bad_reader.next().unwrap() {
        Ok(_) => panic!("This read should not be Ok!"),
        Err(e) => eprintln!("test_bad.fastq: {}", e),
    };
}
