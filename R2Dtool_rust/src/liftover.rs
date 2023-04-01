use multimap::MultiMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use crate::parse_annotation::{Exon, read_gtf_file, read_gff_file};

pub fn convert_transcriptomic_to_genomic_coordinates(
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
