extern crate trio_binning;

use trio_binning::fasta;
use std::fs::File;

fn main() {
    let reader = fasta::Reader::new(File::open("test.fa").unwrap());

    for result in reader {
        let read = result.unwrap();
        println!("ID:{}\tSEQ:{}", read.id(), read.seq());
    }

    if let Ok(bad_record) = fasta::Record::new(&String::from("")) {
        println!("ID:{}\tSEQ:{}", bad_record.id(), bad_record.seq());
    } else { eprintln!("Bad record!!!"); }

    let s = "ABCDEFG\nABCDEFGsdf";
    let lines: Vec<&str> = s.split('\n').collect();
    println!("{:?}", lines);
}
