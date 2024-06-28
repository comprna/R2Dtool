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

    // Process the rest of the lines in parallel
    let results: Vec<_> = input_reader.lines()
        .par_bridge()
        .filter_map(|line| {
            let line = line.ok()?;
            let site_fields: Vec<&str> = line.trim().split('\t').collect();
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

// Unit tests 
#[cfg(test)]
mod tests {
    use super::*;
    use tempfile;
    use crate::parse_gtf::Exon;

    fn create_test_transcript(chromosome: &str, strand: &str, exons: Vec<(u64, u64)>) -> Transcript {
        Transcript {
            transcript_id: "test_transcript".to_string(),
            gene_id: Some("test_gene".to_string()),
            gene_name: Some("test_gene".to_string()),
            utr5_len: None,
            cds_len: None,
            utr3_len: None,
            has_missing_features: false,
            exons: exons.into_iter().map(|(start, end)| Exon {
                seq_id: chromosome.to_string(),
                source: "test".to_string(),
                feature_type: "exon".to_string(),
                start,
                end,
                length: end - start + 1,
                score: None,
                strand: strand.to_string(),
                frame: None,
                attributes: HashMap::new(),
                feature: Some("exon".to_string()),
            }).collect(),
            biotype: Some("protein_coding".to_string()),
            splice_junction_positions: Vec::new(),
            strand: Some(strand.to_string()),
            chromosome: chromosome.to_string(),
            cds_starts: Vec::new(),
            cds_ends: Vec::new(),
            transcript_length: None,
        }
    }

    #[test]
    fn test_convert_transcriptomic_to_genomic_coordinates_positive_strand() {
        let mut annotations = HashMap::new();
        annotations.insert("transcript1".to_string(), create_test_transcript("chr1", "+", vec![(100, 200), (300, 400)]));
        
        let site_fields = vec!["transcript1", "50", "A", "T"];
        let result = convert_transcriptomic_to_genomic_coordinates(&site_fields, &annotations, false);
        assert_eq!(result, Some("chr1\t149\t150\t\t\t+\ttranscript1\t50\tA\tT".to_string()));
    }

    #[test]
    fn test_convert_transcriptomic_to_genomic_coordinates_negative_strand() {
        let mut annotations = HashMap::new();
        annotations.insert("transcript1".to_string(), create_test_transcript("chr1", "-", vec![(100, 200), (300, 400)]));
        
        let site_fields = vec!["transcript1", "50", "A", "T"];
        let result = convert_transcriptomic_to_genomic_coordinates(&site_fields, &annotations, false);
        assert_eq!(result, Some("chr1\t349\t350\t\t\t-\ttranscript1\t50\tA\tT".to_string()));
    }

    #[test]
    fn test_convert_transcriptomic_to_genomic_coordinates_with_version() {
        let mut annotations = HashMap::new();
        annotations.insert("transcript1.1".to_string(), create_test_transcript("chr1", "+", vec![(100, 200)]));
        
        let site_fields = vec!["transcript1.1", "50", "A", "T"];
        let result = convert_transcriptomic_to_genomic_coordinates(&site_fields, &annotations, true);
        assert_eq!(result, Some("chr1\t149\t150\t\t\t+\ttranscript1.1\t50\tA\tT".to_string()));
    }

    #[test]
    fn test_convert_transcriptomic_to_genomic_coordinates_invalid_transcript() {
        let annotations = HashMap::new();
        let site_fields = vec!["invalid_transcript", "50", "A", "T"];
        let result = convert_transcriptomic_to_genomic_coordinates(&site_fields, &annotations, false);
        assert_eq!(result, None);
    }

    #[test]
    fn test_convert_transcriptomic_to_genomic_coordinates_out_of_bounds() {
        let mut annotations = HashMap::new();
        annotations.insert("transcript1".to_string(), create_test_transcript("chr1", "+", vec![(100, 200)]));
        
        let site_fields = vec!["transcript1", "150", "A", "T"];
        let result = convert_transcriptomic_to_genomic_coordinates(&site_fields, &annotations, false);
        assert_eq!(result, None);
    }

    #[test]
    fn test_convert_transcriptomic_to_genomic_coordinates_real_data() {
        let mut annotations = HashMap::new();
        annotations.insert("ENST00000400109.2".to_string(), create_test_transcript("13", "-", vec![(19304593, 19305625)]));

        let site_fields = vec!["ENST00000400109.2", "87", "88", "a", "10", "+", "10", "0.00"];
        let result = convert_transcriptomic_to_genomic_coordinates(&site_fields, &annotations, true);
        
        assert_eq!(
            result,
            Some("13\t19305537\t19305538\t\t\t-\tENST00000400109.2\t87\t88\ta\t10\t+\t10\t0.00".to_string())
        );
    }

    #[test]
    fn test_run_liftover() {
        use clap::{Arg, Command};
    
        // Mock input data
        let input_data = "transcript\tstart\tend\tbase\tcoverage\tstrand\tN_valid_cov\tfraction_modified\nENST00000400109.2\t87\t88\ta\t10\t+\t10\t0.00\n";
        let gtf_data = "13\thavana\ttranscript\t19304593\t19305625\t.\t-\t.\tgene_id \"ENSG00000215349\"; gene_version \"2\"; transcript_id \"ENST00000400109\"; transcript_version \"2\"; gene_name \"MRPL3P1\"; gene_source \"havana\"; gene_biotype \"processed_pseudogene\"; transcript_name \"MRPL3P1-201\"; transcript_source \"havana\"; transcript_biotype \"processed_pseudogene\"; tag \"basic\"; tag \"Ensembl_canonical\"; transcript_support_level \"NA\";\n13\thavana\texon\t19304593\t19305625\t.\t-\t.\tgene_id \"ENSG00000215349\"; gene_version \"2\"; transcript_id \"ENST00000400109\"; transcript_version \"2\"; exon_number \"1\"; gene_name \"MRPL3P1\"; gene_source \"havana\"; gene_biotype \"processed_pseudogene\"; transcript_name \"MRPL3P1-201\"; transcript_source \"havana\"; transcript_biotype \"processed_pseudogene\"; exon_id \"ENSE00001541585\"; exon_version \"2\"; tag \"basic\"; tag \"Ensembl_canonical\"; transcript_support_level \"NA\";\n";
    
        // Create temporary files
        let temp_dir = tempfile::tempdir().unwrap();
        let input_path = temp_dir.path().join("input.txt");
        let gtf_path = temp_dir.path().join("test.gtf");
        let output_path = temp_dir.path().join("output.txt");
        println!("{:?}", temp_dir.path());
        std::fs::write(&input_path, input_data).unwrap();
        std::fs::write(&gtf_path, gtf_data).unwrap();
    
        // Create mock ArgMatches
        let matches = Command::new("test")
            .arg(Arg::new("gtf").short('g').long("gtf").required(true))
            .arg(Arg::new("input").short('i').long("input").required(true))
            .arg(Arg::new("output").short('o').long("output"))
            .get_matches_from(vec![
                "test",
                "-g", gtf_path.to_str().unwrap(),
                "-i", input_path.to_str().unwrap(),
                "-o", output_path.to_str().unwrap(),
            ]);
    
        // Run liftover
        let result = run_liftover(&matches, true, false);
        assert!(result.is_ok());
    
        // Check output
        let output_content = std::fs::read_to_string(output_path).unwrap();

        // Print the output content
        println!("Output content: {}", output_content);

        let expected_output = "chromosome\tstart\tend\tname\tscore\tstrand\ttranscript\tstart\tend\tbase\tcoverage\tstrand\tN_valid_cov\tfraction_modified\n13\t19305537\t19305538\t\t\t-\tENST00000400109.2\t87\t88\ta\t10\t+\t10\t0.00\n";
        assert_eq!(output_content, expected_output);
    }
}

