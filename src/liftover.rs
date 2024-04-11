use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::error::Error;
// use crate::parse_annotation::{Transcript, read_gtf_file, read_gff_file};
use crate::parse_gtf::{Transcript, read_annotation_file};
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

    // Parse the transcriptomic position (second field) as u64
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

            if current_position + exon_length > position {

                // Calculate the genomic position based on the strand
                let genomic_position = if transcript.strand.as_deref() == Some("+") {
                    position - current_position + exon_data.start - 1
                } else {
                    exon_data.end - (position - current_position + 2) + 1
                };

                // Get the chromosome name
                let chrom = &transcript.chromosome;

                // Get the genomic strand
                let genomic_strand = &transcript.strand;

                // Keep all original columns in the output of the liftover 
                let additional_columns = site_fields.join("\t");

                // Return the formatted output string
                return Some(format!(
                    "{}\t{}\t{}\t\t\t{}\t{}\t{}",
                    chrom, genomic_position, genomic_position + 1, genomic_strand.as_deref().unwrap_or(""), site_fields[0], additional_columns
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

pub fn run_liftover(matches: &clap::ArgMatches, has_header: bool) -> Result<(), Box<dyn Error>> {

    // TODO: implement format matching for GFF3 file parsing 
    // let default_format = String::from("gtf");
    //let format = matches.get_one("format").unwrap_or(&default_format);

    let gtf_file: String = matches.get_one::<String>("gtf").unwrap().to_string();
    let input_file: String = matches.get_one::<String>("input").unwrap().to_string();
    let output_file: Option<String> = matches.get_one::<String>("output").map(|s: &String| s.to_string());
    
    // By default, read in the annotations as GTF file
    // TODO: implement GFF3 parsing 
    let annotations = read_annotation_file(&gtf_file, true)?;

    // Print the annotations in a table
    // eprintln!("Previewing transcript annotations\n");
    // preview_annotations(&annotations);

    let mut input_reader = BufReader::new(File::open(input_file.clone()).unwrap_or_else(|_| panic!("Cannot open input file: {}", input_file)));

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
            if let Err(_) = writeln!(output_writer, "{}", genomic_coordinates) {
                break;
            }
        }
        line.clear();
    }
    Ok(())
}

pub fn preview_annotations(annotations: &HashMap<String, Transcript>) {
    println!("Number of annotations: {}", annotations.len()); // Add this line
    for (key, transcript) in annotations {
        println!("Annotations start");
        println!("Transcript ID: {}", key);
        println!("{:#?}", transcript);
        println!("Annotations end");
    }
}

pub fn print_exon_info(annotations: &HashMap<String, Transcript>) {
    for (key, transcript) in annotations {
        println!("Transcript ID: {}", key);
        for (i, exon_data) in transcript.exons.iter().enumerate() {
            println!("Exon {}: start={}, end={}", i + 1, exon_data.start, exon_data.end);
        }
    }
}
