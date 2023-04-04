use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use crate::parse_annotation::{Transcript, read_gtf_file, read_gff_file};
use std::collections::HashMap;

pub fn convert_transcriptomic_to_genomic_coordinates(
    site_fields: &[&str], // input: tab-separated transcriptome position fields
    annotations: &HashMap<String, Transcript>, // input: parsed annotation object
) -> Option<String> { // return type: Option<String), either Some(String) or None
    // Check if there are at least 4 fields in site_fields
    if site_fields.len() < 4 {
        return None; // If not, return None
    }

    // Extract the transcript ID with version from the first field
    let transcript_id_with_version = site_fields[0];
    // Remove the version from transcript ID
    let transcript_id = transcript_id_with_version.split('.').next().unwrap();

    // Parse the transcriptomic position (second field) as a u64
    let position: u64 = site_fields[1].parse().unwrap();

    // Initialize the current_position variable to keep track of the cumulative exon lengths
    let mut current_position = 0;

    // Check if there is a transcript associated with the given transcript_id
    if let Some(transcript) = annotations.get(transcript_id) {
        let exons = if transcript.strand.as_deref() == Some("-") {
            transcript.exons.iter().rev().collect::<Vec<_>>()
        } else {
            transcript.exons.iter().collect::<Vec<_>>()
        };

        // Iterate through each exon in the transcript
        for exon_data in &exons {
            // Calculate the exon length
            let exon_length = exon_data.end - exon_data.start + 1;
            // Check if the current exon contains the transcriptomic position
            if current_position + exon_length >= position {
                // Calculate the genomic position based on the strand
                let genomic_position = if transcript.strand.as_deref() == Some("+") {
                    position - current_position + exon_data.start
                } else {
                    exon_data.end - (position - current_position + 2)
                };
                // Get the chromosome name
                let chrom = &transcript.chromosome;
                // Get the genomic strand
                let genomic_strand = &transcript.strand;
                // Join the additional columns (fields after the second one) using a tab character
                let additional_columns = site_fields[2..].join("\t");
                // Return the formatted output string
                return Some(format!(
                    "{}\t{}\t{}\t\t\t{}\t{}",
                    chrom, genomic_position, genomic_position + 1, genomic_strand.as_deref().unwrap_or(""), additional_columns
                ));
            }
            // Increment the current_position by the exon_length
            current_position += exon_length;
        }
    }

    // If no suitable transcript is found, print a warning and return None
    eprintln!("Warning: No associated transcripts found for site '{}'.", transcript_id);
    None
}



pub fn run_liftover(matches: &clap::ArgMatches) {
    let gff_file = matches.value_of("gff").unwrap();
    let input_file = matches.value_of("input").unwrap();
    let output_file = matches.value_of("output");
    let format = matches.value_of("format").unwrap_or("gff");

    let annotations = if format == "gtf" {
        read_gtf_file(gff_file)
    } else {
        read_gff_file(gff_file)
    };

    // Print the annotations in a table
    // preview_annotations(&annotations);
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

        // Update the header for the output file
        let output_header = format!(
            "chromosome\tstart\tend\tname\tscore\tstrand\t{}",
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

// pub fn preview_annotations(annotations: &HashMap<String, Transcript>) {
//     println!("Number of annotations: {}", annotations.len()); // Add this line
//     for (key, transcript) in annotations {
//         println!("Annotations start");
//         println!("Transcript ID: {}", key);
//         println!("{:#?}", transcript);
//         println!("Annotations end");
//     }
// }

pub fn print_exon_info(annotations: &HashMap<String, Transcript>) {
    for (key, transcript) in annotations {
        println!("Transcript ID: {}", key);
        for (i, exon_data) in transcript.exons.iter().enumerate() {
            println!("Exon {}: start={}, end={}", i + 1, exon_data.start, exon_data.end);
        }
    }
}
