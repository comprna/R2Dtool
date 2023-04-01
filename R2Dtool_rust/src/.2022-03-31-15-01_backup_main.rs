use bio::io::gff;
use clap::{App, Arg};
use multimap::MultiMap;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

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

    // Print the site fields for inspection while debugging
    println!("Site fields: {:?}", site_fields);

    let mut transcript_id = site_fields[0];
    let position: u64 = site_fields[1].parse().unwrap();
    let mut current_position = 0;

    if let Some(exons) = annotations.get_vec(transcript_id) {
        for exon in exons {
            let exon_length = exon.end - exon.start + 1;
            if current_position + exon_length >= position {
                let genomic_position = position - current_position + exon.start - 1;
                let chrom = &exon.seq_id;
                let genomic_strand = &exon.strand;
                return Some(format!(
                    "{}\t{}\t{}\t{}",
                    chrom, genomic_position, site_fields[2], genomic_strand
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
        Arg::with_name("output")
        .short("o")
        .long("output")
        .value_name("OUTPUT_FILE")
        .help("Sets the output file with genomic coordinates")
        .required(true)
        .takes_value(true),
    )
    .arg(
        Arg::with_name("format")
        .short("f")
        .long("format")
        .value_name("FORMAT")
        .help("Specify the file format: gtf or gff (default: gff)")
        .takes_value(true),
    )
    .get_matches();

    let gff_file = matches.value_of("gff").unwrap();
    let input_file = matches.value_of("input").unwrap();
    let output_file = matches.value_of("output").unwrap();

    let format = matches.value_of("format").unwrap_or("gff");


    let annotations = match format.to_lowercase().as_str() {
        "gtf" => read_gtf_file(gff_file),
        "gff" => read_gff_file(gff_file),
        _ => panic!("Invalid file format provided. Please use 'gtf' or 'gff'."),
    };

    // call the new function to print exon information
    print_exon_info(&annotations);
    std::process::exit(0);

    // debugging the annotation import
    preview_annotations(&annotations);

    // Preview site fields from the input file
    let input = BufReader::new(File::open(input_file).unwrap());
    println!("Previewing site fields from the input file:");
    for (i, line) in input.lines().enumerate() {
        if i >= 10 {
            break;
        }
        let line = line.unwrap();
        let site_fields: Vec<&str> = line.split('\t').collect();
        println!("Site fields: {:?}", site_fields);
    }

    let input = BufReader::new(File::open(input_file).unwrap());
    let mut output = File::create(output_file).unwrap();

    // for line in input.lines() {
    //     let line = line.unwrap();
    //     let site_fields: Vec<&str> = line.split('\t').collect();
    //
    //     if let Some(genomic_coordinates) = convert_transcriptomic_to_genomic_coordinates(&site_fields, &annotations) {
    //         writeln!(output, "{}", genomic_coordinates).unwrap();
    //     }
    // }

}
