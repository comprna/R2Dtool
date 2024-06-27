use std::fs::File;
use std::io::{BufRead, BufReader, Write, BufWriter};
use std::error::Error;
use crate::parse_gtf::{Transcript, read_annotation_file};
use std::collections::HashMap;
use rayon::prelude::*;

fn convert_transcriptomic_to_genomic_coordinates(
    site_fields: &[&str],
    annotations: &HashMap<String, Transcript>,
    has_version: bool
) -> Option<String> {
    if site_fields.len() < 4 {
        return None;
    }

    let transcript_id_with_version = site_fields[0];
    let transcript_id = if has_version {
        transcript_id_with_version
    } else {
        transcript_id_with_version.split('.').next()?
    };

    let position: u64 = site_fields[1].parse().ok()?;
    let mut current_position = 0;

    let transcript = annotations.get(transcript_id)?;
    let mut exons = transcript.exons.clone();
    if transcript.strand.as_deref() == Some("-") {
        exons.sort_by(|a, b| b.start.cmp(&a.start));
    } else {
        exons.sort_by(|a, b| a.start.cmp(&b.start));
    }

    for exon_data in &exons {
        let exon_length = exon_data.end - exon_data.start + 1;
        if current_position + exon_length > position {
            let genomic_position = if transcript.strand.as_deref() == Some("+") {
                position - current_position + exon_data.start - 1
            } else {
                exon_data.end - (position - current_position) - 1
            };

            let chrom = &transcript.chromosome;
            let genomic_strand = &transcript.strand;
            let additional_columns = site_fields[1..].join("\t");

            return Some(format!(
                "{}\t{}\t{}\t\t\t{}\t{}\t{}",
                chrom, genomic_position, genomic_position + 1, genomic_strand.as_deref().unwrap_or(""), site_fields[0], additional_columns
            ));
        }
        current_position += exon_length;
    }

    eprintln!("Warning: No associated transcripts found for site '{}'.", transcript_id);
    None
}
pub fn run_liftover(matches: &clap::ArgMatches, has_header: bool, has_version: bool) -> Result<(), Box<dyn Error>> {

    let gtf_file: String = matches.get_one::<String>("gtf").unwrap().to_string();
    let input_file: String = matches.get_one::<String>("input").unwrap().to_string();
    let output_file: Option<String> = matches.get_one::<String>("output").map(|s: &String| s.to_string());
    
    let annotations = read_annotation_file(&gtf_file, true, has_version)?;

    let mut input_reader = BufReader::with_capacity(512 * 1024, File::open(input_file.clone())?);
    let mut output_writer = BufWriter::with_capacity(512 * 1024, match output_file {
        Some(file_name) => Box::new(File::create(file_name)?) as Box<dyn Write>,
        None => Box::new(std::io::stdout()) as Box<dyn Write>,
    });

    if has_header {
        let mut header = String::new();
        input_reader.read_line(&mut header)?;
        writeln!(output_writer, "chromosome\tstart\tend\tname\tscore\tstrand\t{}", header.trim())?;
    }

    // Process lines in parallel
    let results: Vec<_> = input_reader.lines()
        .par_bridge()
        .filter_map(|line| {
            let line = line.ok()?;
            let site_fields: Vec<&str> = line.trim().split('\t').collect();
            if has_header && site_fields[0] == "transcript" {
                return None;
            }
            convert_transcriptomic_to_genomic_coordinates(&site_fields, &annotations, has_version)
        })
        .collect();

    // Write results
    for result in results {
        writeln!(output_writer, "{}", result)?;
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
