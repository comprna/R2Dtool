#![allow(non_snake_case)]
use clap::{Arg,Command};

// modules
pub mod parse_annotation;
pub mod annotate;
pub mod liftover;
pub mod parse_gtf;

fn main() {
    let matches = Command::new("R2Dtool")
        .version("1.0.0")
        .author("AJ Sethi, compRNA <aditya.sethi@anu.edu.au; CompRNA@ANU365.onmicrosoft.com>")
        .about("R2Dtool")
        .arg_required_else_help(true) 
        .subcommand(
            Command::new("liftover")
                .about("Converts transcriptomic to genomic coordinates")
                .arg(
                    Arg::new("gtf")
                    .short('g')
                    .long("gtf")
                    .value_name("GTF_FILE")
                    .help("Path to GTF gene structure annotation")
                    .required(true)
                )
                .arg(
                    Arg::new("input")
                    .short('i')
                    .long("input")
                    .value_name("INPUT_FILE")
                    .help("Path to input file with transcriptomic coordinates")
                    .required(true)
                )
                .arg(
                    Arg::new("header")
                    .short('H')
                    .long("header")
                    .help("Indicates that the input file has a header in line 1")
                    .action(clap::ArgAction::SetTrue)
                )
                .arg(
                    Arg::new("output")
                    .short('o')
                    .long("output")
                    .value_name("OUTPUT_FILE")
                    .help("Path to output file")
                )
                // .arg(
                //     Arg::new("format")
                //     .short('f')
                //     .long("format")
                //     .value_name("FORMAT")
                //     .help("Specify the gene strutcture annotation format: gtf or gff (default: gff)")
                // )
        )
        .subcommand(
            Command::new("annotate")
                .about("Annotates transcriptomic sites with genomic cooridnates")
                .arg(
                    Arg::new("gtf")
                    .short('g')
                    .long("gtf")
                    .value_name("GTF_FILE")
                    .help("Path to GTF gene structure annotation")
                    .required(true)
                )
                .arg(
                    Arg::new("input")
                    .short('i')
                    .long("input")
                    .value_name("INPUT_FILE")
                    .help("Path to input file with transcriptomic coordinates")
                    .required(true)
                )
                .arg(
                    Arg::new("header")
                    .short('H')
                    .long("header")
                    .help("Indicates that the input file has a header in line 1")
                    .action(clap::ArgAction::SetTrue)
                )
                .arg(
                    Arg::new("output")
                    .short('o')
                    .long("output")
                    .value_name("OUTPUT_FILE")
                    .help("Path to output file")
                )
                // .arg(
                //     Arg::new("format")
                //     .short('f')
                //     .long("format")
                //     .value_name("FORMAT")
                //     .help("Specify the gene strutcture annotation format: gtf or gff (default: gff)")
                // )
        )
        .get_matches();

    // Handle the liftover subcommand
    if let Some(liftover_matches) = matches.subcommand_matches("liftover") {
        let has_header = liftover_matches.get_flag("header");
        
        
        eprintln!("Running liftover...");

        if let Err(e) = liftover::run_liftover(liftover_matches, has_header) {
            eprintln!("Error running liftover: {}", e);
        }
    }
    

    // Handle the annotate subcommand
    if let Some(annotate_matches) = matches.subcommand_matches("annotate") {
        let has_header = annotate_matches.get_flag("header");
        
        eprintln!("Running annotate...");
        
        if let Err(e) = annotate::run_annotate(annotate_matches, has_header) {
            eprintln!("Error running annotate: {}", e);
        }
    }
}
