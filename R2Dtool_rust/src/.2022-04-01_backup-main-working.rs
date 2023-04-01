use bio::io::gff;
use clap::{App, Arg};
use multimap::MultiMap;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

#[derive(Debug, PartialEq)]
struct Exon {
    seq_id: String,
    source: String,
    feature_type: String,
    start: u64,
    end: u64,
    score: Option<f64>,
    strand: String,
    frame: Option<char>,
    attributes: HashMap<String, String>,
}

fn read_gff_file(gff_file: &str) -> MultiMap<String, Exon> {
    let mut annotations = MultiMap::new();
    let mut reader = gff::Reader::from_file(gff_file, gff::GffType::GFF3).unwrap();

    for record in reader.records() {
        let record = record.unwrap();
        let attributes = record.attributes().clone();

        if record.feature_type() == "exon" {
            let exon = Exon {
                seq_id: record.seqname().to_string(),
                source: record.source().to_string(),
                feature_type: record.feature_type().to_string(),
                start: *record.start(),
                end: *record.end(),
                score: record.score().map(|s| s as f64),
                strand: record.strand().expect("Missing strand information").strand_symbol().to_string(),
                frame: record.frame().chars().next(),
                attributes: multimap_to_hashmap(&attributes),
            };
            if let Some(transcript_id) = exon.attributes.get("Parent") {
                annotations.insert(transcript_id.to_string(), exon);
            }
        }
    }

    annotations
}

fn multimap_to_hashmap(multimap: &MultiMap<String, String>) -> HashMap<String, String> {
    let mut hashmap = HashMap::new();
    for (key, value) in multimap.iter() {
        hashmap.insert(key.clone(), value.clone());
    }
    hashmap
}

fn parse_gff_attributes(attributes: &MultiMap<String, String>) -> HashMap<String, String> {
    let mut attr_map = HashMap::new();
    for (key, value) in attributes.iter() {
        attr_map.insert(key.clone(), value.clone());
    }
    attr_map
}

fn preview_annotations(annotations: &MultiMap<String, Exon>) {
    println!("Previewing the first transcript:");

    for (key, exons) in annotations.iter_all() {
        println!("Transcript ID: {}", key);
        // for exon in exons {
        //     println!("  Exon: {:?}", exon);
        // }
        // break; // Add this line to exit the loop after the first transcript
    }
}

fn print_exon_info(annotations: &MultiMap<String, Exon>) {
    // Print table headers
    println!(
        "{:<10} {:<10} {:<15} {:<15} {:<7}",
        "Start", "End", "Transcript ID", "Gene ID", "Strand"
    );
    println!("{}", "-".repeat(57));

    let na_gene_id = String::from("N/A");

    for (transcript_id, exon) in annotations.iter() {
        let gene_id_str = exon.attributes.get("gene_id").unwrap_or(&na_gene_id);

        println!(
            "{:<10} {:<10} {:<15} {:<15} {:<7}",
            exon.start, exon.end, transcript_id, gene_id_str, exon.strand
        );
    }
}

fn read_gtf_file(gtf_file: &str) -> MultiMap<String, Exon> {
    let mut annotations = MultiMap::new();
    let mut reader = gff::Reader::from_file(gtf_file, gff::GffType::GTF2).unwrap();

    for record in reader.records() {
        let record = record.unwrap();

        if record.feature_type() == "#" {
            continue;
        }

        if record.feature_type() == "exon" {
            let attributes = parse_gff_attributes(record.attributes());

            let exon = Exon {
                seq_id: record.seqname().to_string(),
                source: record.source().to_string(),
                feature_type: record.feature_type().to_string(),
                start: *record.start(),
                end: *record.end(),
                score: record.score().map(|s| s as f64),
                strand: record.strand().expect("Missing strand information").strand_symbol().to_string(),
                frame: record.frame().chars().next(),
                attributes: attributes,
            };
            if let Some(transcript_id) = exon.attributes.get("transcript_id") {
                annotations.insert(transcript_id.to_string(), exon);
            }
        }
    }

    annotations
}


fn convert_transcriptomic_to_genomic_coordinates(
    site_fields: &[&str],
    annotations: &MultiMap<String, Exon>,
) -> Option<String> {
    if site_fields.len() < 4 {
        return None;
    }

    let transcript_id_with_version = site_fields[0];
    let transcript_id = transcript_id_with_version.split('.').next().unwrap(); // Add this line to remove version from transcript ID
    let position: u64 = site_fields[1].parse().unwrap();
    let mut current_position = 0;

    if let Some(exons) = annotations.get_vec(transcript_id) {
        for exon in exons {
            let exon_length = exon.end - exon.start + 1;
            if current_position + exon_length >= position { // Change this line
                let genomic_position = if exon.strand == "+" {
                    position - current_position + exon.start
                } else {
                    exon.end - (position - current_position) - 2
                };
                let chrom = &exon.seq_id;
                let genomic_strand = &exon.strand;
                let additional_columns = site_fields[2..].join("\t");
                return Some(format!(
                    "{}\t{}\t{}\t\t\t{}\t{}",
                    chrom, genomic_position, genomic_position + 1, genomic_strand, additional_columns
                ));
            }
            current_position += exon_length;
        }
    }

    None
}



fn main() {
    let matches = App::new("R2Dtool")
    .version("0.1.0")
    .author("Your Name <your.email@example.com>")
    .about("Converts transcriptomic coordinates to genomic coordinates")
    .arg(
        Arg::with_name("gff")
        .short("g")
        .long("gff")
        .value_name("GFF_FILE")
        .help("Sets the GFF annotation file")
        .required(true)
        .takes_value(true),
    )
    .arg(
        Arg::with_name("input")
        .short("i")
        .long("input")
        .value_name("INPUT_FILE")
        .help("Sets the input file with transcriptomic coordinates")
        .required(true)
        .takes_value(true),
    )
    .arg(
        Arg::with_name("header")
        .short("H")
        .long("header")
        .help("Indicates that the input file has a header")
        .takes_value(false),
    )
    .arg(
        Arg::with_name("output")
        .short("o")
        .long("output")
        .value_name("OUTPUT_FILE")
        .help("Sets the output file with genomic coordinates")
        .takes_value(true),
    )
    .arg(
        Arg::with_name("format")
        .short("f")
        .long("format")
        .value_name("FORMAT")
        .help("Specify the file format: gtf or gff (default: g
            ff)")
            .takes_value(true),
        )
        .get_matches();

        let gff_file = matches.value_of("gff").unwrap();
        let input_file = matches.value_of("input").unwrap();
        let output_file = matches.value_of("output");
        let format = matches.value_of("format").unwrap_or("gff");

        let annotations = if format == "gtf" {
            read_gtf_file(gff_file)
        } else {
            read_gff_file(gff_file)
        };

        // preview_annotations(&annotations);
        // print_exon_info(&annotations);
        // std::process::exit(0);

        let has_header = matches.is_present("header");

        let mut input_reader = BufReader::new(File::open(input_file).unwrap());

        let mut output_writer: Box<dyn Write> = match output_file {
            Some(file_name) => Box::new(File::create(file_name).unwrap()),
            None => Box::new(std::io::stdout()),
        };

        let mut header = String::new();
        if has_header {
            input_reader.read_line(&mut header).unwrap();
            let header_fields: Vec<&str> = header.trim().split('\t').collect();
            // println!("Header: {:?}", header_fields);

            // Update the header for the output file
            let output_header = format!(
                "Chromosome\tStart\tEnd\tName\tScore\tStrand\t{}",
                header_fields[2..].join("\t")
            );
            writeln!(output_writer, "{}", output_header).unwrap();
        }

        let mut line = String::new();
        while input_reader.read_line(&mut line).unwrap() > 0 {
            let site_fields: Vec<&str> = line.trim().split('\t').collect();
            if has_header && site_fields[0] == "transcript" {
                continue;
            }
            if let Err(e) = site_fields[1].parse::<u64>() {
                eprintln!("Error parsing position from line: '{}'\nError: {}", line.trim(), e);
            } else if let Some(genomic_coordinates) =
            convert_transcriptomic_to_genomic_coordinates(&site_fields, &annotations)
            {
                writeln!(output_writer, "{}", genomic_coordinates).unwrap();
            }
            line.clear();
        }
    }
