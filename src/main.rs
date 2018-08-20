extern crate trio_binning;

use trio_binning::seq::{fasta, fastq};
use std::fs::File;

fn main() {
    let reader = fasta::Reader::new(File::open("test.fa").unwrap());

    for result in reader {
        let read = result.unwrap();
        println!("ID:{}\tSEQ:{}", read.id(), read.seq());
    }

    match fasta::Record::new(&String::from("")) {
        Ok(_) => panic!("Should not get Ok for an empty record!"),
        Err(e) => eprintln!("ERROR: {}", e),
    }

    let s = "ABCDEFG\nABCDEFGsdf";
    let lines: Vec<&str> = s.split('\n').collect();
    println!("{:?}", lines);

    let bad_reader = fasta::Reader::new(File::open("test_bad.fa").unwrap());
    for result in bad_reader {
        match result {
            Ok(r) => 
                panic!("Should not get Ok for a misformatted fasta file! {:?}",
                       r),
            Err(e) => eprintln!("test_bad.fa: {}", e),
        }
    }

    let fastq_reader = fastq::Reader::new(File::open("test.fastq").unwrap());

    for result in fastq_reader {
        let read = result.unwrap();
        println!("ID:{}\tSEQ:{}\tQUAL:{}", read.id(), read.seq(), read.qual());
    }

    let mut bad_fq_reader = fastq::Reader::new(
        File::open("test_bad.fastq").unwrap());
    let good_read = bad_fq_reader.next().unwrap().unwrap();
    println!("ID:{}\tSEQ:{}\tQUAL:{}", good_read.id(), good_read.seq(),
            good_read.qual());
    match bad_fq_reader.next().unwrap() {
        Ok(_) => panic!("This read should not be Ok!"),
        Err(e) => eprintln!("test_bad.fastq: {}", e),
    };
}
