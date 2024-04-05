#![allow(non_snake_case)]
use clap::{Arg,Command};

// modules
pub mod parse_annotation;
pub mod annotate;
pub mod liftover;

fn main() {
    let matches = Command::new("R2Dtool")
        .version("1.0.0")
        .author("AJ Sethi, compRNA <aditya.sethi@anu.edu.au; CompRNA@ANU365.onmicrosoft.com>")
        .about("R2Dtool")
        .subcommand(
            Command::new("liftover")
                .about("Converts transcriptomic to genomic coordinates")
                .arg(
                    Arg::new("gff")
                    .short('g')
                    .long("gff")
                    .value_name("GFF_FILE")
                    .help("Path to gene structure annotation")
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
                )
                .arg(
                    Arg::new("output")
                    .short('o')
                    .long("output")
                    .value_name("OUTPUT_FILE")
                    .help("Path to output file")
                )
                .arg(
                    Arg::new("format")
                    .short('f')
                    .long("format")
                    .value_name("FORMAT")
                    .help("Specify the gene strutcture annotation format: gtf or gff (default: gff)")
                )
        )
        .subcommand(
            Command::new("Annotate")
                .about("Annotates transcriptomic sites")
                .arg(
                    Arg::new("gff")
                    .short('g')
                    .long("gff")
                    .value_name("GFF_FILE")
                    .help("Path to gene structure annotation")
                    .required(true)
                )
                .arg(
                    Arg::new("input transcriptome-mapped sites")
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
                )
                .arg(
                    Arg::new("output")
                    .short('o')
                    .long("output")
                    .value_name("OUTPUT_FILE")
                    .help("Path to output file")
                )
                .arg(
                    Arg::new("format")
                    .short('f')
                    .long("format")
                    .value_name("FORMAT")
                    .help("Specify the gene strutcture annotation format: gtf or gff (default: gff)")
                )
        )
        .get_matches();

    // Handle the liftover subcommand
    if let Some(liftover_matches) = matches.subcommand_matches("liftover") {
        // Call your liftover function with liftover_matches as an argument
        eprintln!("Running liftover...");
        liftover::run_liftover(&liftover_matches);
    }

    // Handle the annotate subcommand
    if let Some(annotate_matches) = matches.subcommand_matches("annotate") {
        // Call your annotate function with annotate_matches as an argument
        eprintln!("Running annotate...");
        annotate::run_annotate(&annotate_matches);
    }

}

//
// fn main() {
//     let matches = App::new("R2Dtool")
//     .version("0.1.0")
//     .author("Your Name your.email@example.com")
//     .about("Converts transcriptomic coordinates to genomic coordinates")
//     .arg(
//         Arg::new("gff")
//         .short("g")
//         .long("gff")
//         .value_name("GFF_FILE")
//         .help("Path to gene structure annotation")
//         .required(true)
//         .takes_value(true),
//     )
//     .arg(
//         Arg::new("input")
//         .short("i")
//         .long("input")
//         .value_name("INPUT_FILE")
//         .help("Path to input file with transcriptomic coordinates")
//         .required(true)
//         .takes_value(true),
//     )
//     .arg(
//         Arg::new("header")
//         .short("H")
//         .long("header")
//         .help("Indicates that the input file has a header in line 1")
//         .takes_value(false),
//     )
//     .arg(
//         Arg::new("output")
//         .short("o")
//         .long("output")
//         .value_name("OUTPUT_FILE")
//         .help("Path to output file")
//         .takes_value(true),
//     )
//     .arg(
//         Arg::new("format")
//         .short("f")
//         .long("format")
//         .value_name("FORMAT")
//         .help("Specify the gene strutcture annotation format: gtf or gff (default: gff)")
//         .takes_value(true),
//     )
//     .arg(
//         Arg::new("subcommand")
//         .help("Subcommand to run: liftover or annotate")
//         .required(true)
//         .index(1),
//     )
//     .get_matches();
//
//     // Check the subcommand provided and call the appropriate function
//     match matches.subcommand() {
//         ("liftover", Some(sub_m)) => liftover::run_liftover(sub_m),
//         ("annotate", Some(_sub_m)) => annotate::run_annotate(),
//         _ => {
//             eprintln!("Invalid subcommand provided. Please choose 'liftover' or 'annotate'.");
//             std::process::exit(1);
//         }
//     }
// }