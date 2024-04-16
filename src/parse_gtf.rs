use std::collections::{HashMap, HashSet};
use std::collections::hash_map::Entry;
use multimap::MultiMap;
use bio::io::gff;
use log::{warn, error, debug};

// Define exon structure
#[derive(Debug, PartialEq, Clone, Default)]
pub struct Exon {
    pub seq_id: String,
    pub source: String,
    pub feature_type: String,
    pub start: u64,
    pub end: u64,
    pub length: u64,
    pub score: Option<f64>,
    pub strand: String,
    pub frame: Option<char>,
    pub attributes: HashMap<String, String>,
    pub feature: Option<String>,
}

// Define transcript structure
// each transcript maps to a gene 
// and contains one or more exons
#[derive(Debug, PartialEq, Clone)]
pub struct Transcript {
    pub transcript_id: String,
    pub gene_id: Option<String>,
    pub gene_name: Option<String>,
    pub utr5_len: Option<u64>,
    pub cds_len: Option<u64>,
    pub utr3_len: Option<u64>,
    pub has_missing_features: bool,
    pub exons: Vec<Exon>,
    pub biotype: Option<String>,
    pub splice_junction_positions: Vec<u64>,
    pub strand: Option<String>,
    pub chromosome: String,
    pub cds_starts: Vec<u64>,
    pub cds_ends: Vec<u64>,
    pub transcript_length: Option<u64>, 
}

// Set default traits for transcripts 
impl Default for Transcript {
    fn default() -> Self {
        Self {
            transcript_id: String::new(),
            gene_id: None, 
            gene_name: None, 
            utr5_len: None,
            cds_len: None,
            utr3_len: None,
            has_missing_features: false,
            exons: Vec::new(),
            biotype: None,
            splice_junction_positions: Vec::new(),
            strand: None,
            chromosome: String::new(),    
            cds_starts: Vec::new(),
            cds_ends: Vec::new(),
            transcript_length: None,
        }
    }
}

// Parse gff attributes
pub fn parse_gff_attributes(attributes: &MultiMap<String, String>) -> HashMap<String, String> {
    attributes.iter().map(|(k, v)| (k.clone(), v.clone())).collect()
}

// Read annotation file 
// TODO: Implement parsing for GFF3 files, when is_gtf = False 
pub fn read_annotation_file(file_path: &str, is_gtf: bool) -> Result<HashMap<String, Transcript>, Box<dyn std::error::Error>> {
    let mut transcripts: HashMap<String, Transcript> = HashMap::new();
    let mut ignored_features: HashMap<String, u32> = HashMap::new();

    let mut reader = if is_gtf {
        gff::Reader::from_file(file_path, gff::GffType::GTF2)?
    } else {
        gff::Reader::from_file(file_path, gff::GffType::GFF3)?
    };

    // Define possible keys corresponding to gene biotype 
    let biotype_keys = vec!["transcript_biotype", "transcript_type", "gene_type", "gene_biotype"];

    fn find_biotype(attributes: &HashMap<String, String>, keys: &[&str]) -> Option<String> {
        keys.iter().find_map(|key| attributes.get(*key)).cloned()
    }

    let mut encountered_transcripts: HashSet<String> = HashSet::new();

    for record in reader.records() {
        let record = match record {
            Ok(record) => record,
            Err(err) => {
                error!("Error reading record: {}", err);
                continue;
            }
        };

        debug!("Processing record: {:?}", record);

        if record.feature_type() == "#" {
            continue;
        }

        let attributes = parse_gff_attributes(record.attributes());
        debug!("Parsed attributes: {:?}", attributes);

        let transcript_id_attr = if is_gtf { "transcript_id" } else { "ID" };
        let gene_id_attr = if is_gtf { "gene_id" } else { "Parent" };

        let transcript_id_with_version = match attributes.get(transcript_id_attr) {
            Some(id) => id,
            None => {
                warn!("Missing '{}' in record: {:?}. Skipping...", transcript_id_attr, record);
                continue;
            }
        };
        let transcript_id = match transcript_id_with_version.split('.').next() {
            Some(id) => id.to_string(),
            None => {
                warn!("Invalid transcript ID format: {}. Skipping...", transcript_id_with_version);
                continue;
            }
        };

        debug!("Transcript ID: {}", transcript_id);

        if *record.start() > *record.end() {
            warn!("Invalid coordinates for feature: start > end. Skipping record: {:?}", record);
            continue;
        }

        let feature_length = (*record.end() - *record.start() + 1) as u64;

        match record.feature_type() {
            "exon" => {
                let exon = Exon {
                    seq_id: record.seqname().to_string(),
                    source: record.source().to_string(),
                    feature_type: record.feature_type().to_string(),
                    start: *record.start(),
                    end: *record.end(),
                    score: record.score().map(|s| s as f64),
                    strand: record.strand().map(|s| s.strand_symbol().to_string()).unwrap_or_else(|| ".".to_string()),
                    frame: record.frame().chars().next().filter(|&c| c == '0' || c == '1' || c == '2' || c == '.'),
                    attributes: attributes.clone(),
                    feature: Some(record.feature_type().to_string()),
                    length: feature_length,
                };

                debug!("Exon: {:?}", exon);

                match transcripts.entry(transcript_id.clone()) {
                    Entry::Occupied(mut entry) => {
                        let transcript = entry.get_mut();
                        if transcript.chromosome != exon.seq_id {
                            warn!("Transcript {} has exons on different chromosomes. Removing transcript...", transcript.transcript_id);
                            entry.remove();
                            continue;
                        }
                        transcript.exons.push(exon);
                        debug!("Added exon to existing transcript: {}", transcript_id);
                    }
                    Entry::Vacant(entry) => {
                        let biotype = find_biotype(&attributes, &biotype_keys);
                        let strand = record.strand().map(|s| s.strand_symbol().to_string()).unwrap_or_else(|| ".".to_string());
                        let chromosome = record.seqname().to_string();
                        let gene_id = attributes.get(gene_id_attr).cloned();
                        let gene_name = attributes.get("gene_name").cloned();
                        let transcript = Transcript {
                            transcript_id: transcript_id.clone(),
                            gene_id,
                            gene_name,
                            utr5_len: None,
                            cds_len: None,
                            utr3_len: None,
                            has_missing_features: false,
                            exons: vec![exon],
                            biotype,
                            splice_junction_positions: Vec::new(),
                            strand: Some(strand),
                            chromosome,
                            cds_starts: Vec::new(),
                            cds_ends: Vec::new(),
                            transcript_length: None,
                        };
                        entry.insert(transcript);
                        debug!("Created new transcript: {}", transcript_id);
                    }
                }
            }
            "CDS" => {
                if let Some(transcript) = transcripts.get_mut(&transcript_id) {
                    transcript.cds_len = Some(transcript.cds_len.unwrap_or(0) + feature_length);
                    transcript.cds_starts.push(*record.start());
                    transcript.cds_ends.push(*record.end());
                    debug!("Updated CDS length for transcript: {}", transcript_id);
                }
            }
            "five_prime_utr" | "5UTR" => {
                if let Some(transcript) = transcripts.get_mut(&transcript_id) {
                    transcript.utr5_len = Some(transcript.utr5_len.unwrap_or(0) + feature_length);
                    debug!("Updated 5' UTR length for transcript: {}", transcript_id);
                }
            }
            "three_prime_utr" | "3UTR" => {
                if let Some(transcript) = transcripts.get_mut(&transcript_id) {
                    transcript.utr3_len = Some(transcript.utr3_len.unwrap_or(0) + feature_length);
                    debug!("Updated 3' UTR length for transcript: {}", transcript_id);
                }
            }
            "transcript" | "mRNA" => {
                if !encountered_transcripts.insert(transcript_id.clone()) {
                    warn!("More than one transcript with the same name '{}' found. Please check the annotation file.", transcript_id);
                }
            }
            "UTR" => {
                
                // Construct a UTR feature similar to how you construct an Exon
                let utr_feature = Exon {  
                    seq_id: record.seqname().to_string(),
                    source: record.source().to_string(),
                    feature_type: record.feature_type().to_string(),
                    start: *record.start(),
                    end: *record.end(),
                    length: feature_length,
                    score: record.score().map(|s| s as f64),
                    strand: record.strand().map(|s| s.strand_symbol().to_string()).unwrap_or_else(|| ".".to_string()),
                    frame: None, 
                    attributes: attributes.clone(),
                    feature: Some("UTR".to_string()), // Explicitly marking the feature as UTR
                };
        
                debug!("UTR Feature: {:?}", utr_feature);
        
                // Add the UTR feature to the corresponding transcript
                match transcripts.entry(transcript_id.clone()) {
                    Entry::Occupied(mut entry) => {
                        let transcript = entry.get_mut();
                        transcript.exons.push(utr_feature); // Here we are adding UTRs to the exons vector; consider renaming or creating a separate vector for clarity
                        debug!("Added UTR to existing transcript: {}", transcript_id);
                    }
                    Entry::Vacant(entry) => {
                        // Similar logic as for exons, creating a new transcript if not already present
                        let biotype = find_biotype(&attributes, &biotype_keys);
                        let strand = record.strand().map(|s| s.strand_symbol().to_string()).unwrap_or_else(|| ".".to_string());
                        let chromosome = record.seqname().to_string();
                        let gene_id = attributes.get(gene_id_attr).cloned();
                        let gene_name = attributes.get("gene_name").cloned();
                        let transcript = Transcript {
                            transcript_id: transcript_id.clone(),
                            gene_id,
                            gene_name,
                            utr5_len: None,
                            cds_len: None,
                            utr3_len: None,
                            has_missing_features: false,
                            exons: vec![utr_feature], // Adding the UTR feature
                            biotype,
                            splice_junction_positions: Vec::new(),
                            strand: Some(strand),
                            chromosome,
                            cds_starts: Vec::new(),
                            cds_ends: Vec::new(),
                            transcript_length: None,
                        };
                        entry.insert(transcript);
                        debug!("Created new transcript with UTR: {}", transcript_id);
                    }
                }
            }
        
            other => {
                *ignored_features.entry(other.to_string()).or_insert(0) += 1;
            }
        }
    }

    debug!("Number of transcripts parsed: {}", transcripts.len());

    for transcript in transcripts.values_mut() {
        
        // First, ensure only exon features are considered for splice junctions
        let mut exon_features: Vec<&Exon> = transcript.exons.iter()
            .filter(|exon| exon.feature.as_ref().map_or(false, |ft| ft == "exon"))
            .collect();
    
        // Sort by start position to ensure the correct order for calculating splice positions
        exon_features.sort_by_key(|exon| exon.start);
    
        // Calculate splice junction positions based on sorted exons
        if exon_features.len() > 1 {
            for i in 0..(exon_features.len() - 1) {
                let splice_pos = exon_features[i].end + 1; // Position after the end of the current exon
                transcript.splice_junction_positions.push(splice_pos);
            }
        }
    
        // Calculate the total length of the transcript based on exon lengths
        let total_exon_length: u64 = exon_features
            .iter()
            .map(|exon| exon.length)
            .sum();
        transcript.transcript_length = Some(total_exon_length);
    
        check_transcript_length(transcript);
    }
    
    if !ignored_features.is_empty() {
        warn!("Ignored the following annotation features:");
        for (feature_type, count) in ignored_features {
            warn!("- {}: {} occurrences", feature_type, count);
        }
    }
    
    // Post-processing: Calculate UTR lengths based on their position relative to CDS
    // for transcript in transcripts.values_mut() {
    //     // Ensure features are sorted by their start position
    //     transcript.exons.sort_by_key(|exon| exon.start);
    //     println!("UTR detected 0");
    //     let cds_range = transcript.exons.iter()
    //         .filter(|exon| exon.feature_type == "CDS")
    //         .fold(None, |acc: Option<(u64, u64)>, exon| match acc {
    //             None => Some((exon.start, exon.end)),
    //             Some((start, end)) => Some((start.min(exon.start), end.max(exon.end))),
    //         });
        
    //     println!("UTR detected 1");
    //     if let Some((cds_start, cds_end)) = cds_range {
    //         println!("UTR detected 2");
    //         for feature in &transcript.exons {
    //             if feature.feature_type == "UTR" {
    //                 println!("UTR detected 3");
    //                 // Adjust logic based on strand
    //                 let strand = transcript.strand.as_deref().unwrap_or(".");
    //                 match strand {
    //                     "+" | "." => {
    //                         if feature.end < cds_start {
    //                             // 5'UTR for positive strand
    //                             transcript.utr5_len = Some(transcript.utr5_len.unwrap_or(0) + feature.length);
    //                         } else if feature.start > cds_end {
    //                             // 3'UTR for positive strand
    //                             transcript.utr3_len = Some(transcript.utr3_len.unwrap_or(0) + feature.length);
    //                         }
    //                     },
    //                     "-" => {
    //                         if feature.start > cds_end {
    //                             // 5'UTR for negative strand (logic is reversed)
    //                             transcript.utr5_len = Some(transcript.utr5_len.unwrap_or(0) + feature.length);
    //                         } else if feature.end < cds_start {
    //                             // 3'UTR for negative strand (logic is reversed)
    //                             transcript.utr3_len = Some(transcript.utr3_len.unwrap_or(0) + feature.length);
    //                         }
    //                     },
    //                     _ => {}
    //                 }
    //             }
    //         }
    //     }
    //     }

    for transcript in transcripts.values_mut() {
        //println!("Processing transcript: {}", transcript.transcript_id);
        
        // Sort the CDS start and end vectors to ensure they are in the correct order
        transcript.cds_starts.sort_unstable();
        transcript.cds_ends.sort_unstable();

        // Assume each transcript has at least one CDS; thus, the first CDS start and last CDS end define the CDS range
        if let (Some(&cds_start), Some(&cds_end)) = (transcript.cds_starts.first(), transcript.cds_ends.last()) {
            for feature in &transcript.exons {
                if feature.feature_type == "UTR" {
                    
                    let strand = transcript.strand.as_deref().unwrap_or(".");
                    match strand {
                        "+" | "." => {
                            if feature.end < cds_start {
                                // 5'UTR for positive strand
                                transcript.utr5_len = Some(transcript.utr5_len.unwrap_or(0) + feature.length);
                                // println!("5'UTR detected with length {}", feature.length);

                            } else if feature.start > cds_end {
                                // 3'UTR for positive strand
                                transcript.utr3_len = Some(transcript.utr3_len.unwrap_or(0) + feature.length);
                                // println!("3'UTR detected with length {}", feature.length);
                            }
                        },
                        "-" => {
                            if feature.start > cds_end {
                                // 5'UTR for negative strand (logic is reversed)
                                transcript.utr5_len = Some(transcript.utr5_len.unwrap_or(0) + feature.length);
                                // println!("5'UTR (neg strand) detected with length {}", feature.length);

                            } else if feature.end < cds_start {
                                // 3'UTR for negative strand (logic is reversed)
                                transcript.utr3_len = Some(transcript.utr3_len.unwrap_or(0) + feature.length);
                                // println!("3'UTR (neg strand) detected with length {}", feature.length);
                            }
                        },
                        _ => {}
                    }
                }
            }
        } else {
            // eprintln!("No CDS range calculated for transcript: {}", transcript.transcript_id);
        }
    }
    
    
    debug!("Finished parsing annotation file.");

    Ok(transcripts)
}

// Ensure that exon lengths correspond to CDS + UTR lengths for trancripts that contain CDS 
fn check_transcript_length(transcript: &mut Transcript) {
    // Proceed with the check only if there's annotated CDS or UTR length
    if transcript.cds_len.unwrap_or(0) > 0 || transcript.utr5_len.unwrap_or(0) > 0 || transcript.utr3_len.unwrap_or(0) > 0 {
        let total_exon_length: u64 = transcript.exons.iter().map(|exon| exon.length).sum();

        let transcript_length = transcript.utr5_len.unwrap_or(0)
            + transcript.cds_len.unwrap_or(0)
            + if transcript.cds_len.unwrap_or(0) > 0 { 3 } else { 0 } // Adjust for potential stop codon offset
            + transcript.utr3_len.unwrap_or(0);

        let length_difference = if total_exon_length > transcript_length {
            total_exon_length as i64 - transcript_length as i64
        } else {
            transcript_length as i64 - total_exon_length as i64
        };

        if total_exon_length != transcript_length {
            if length_difference == 3 {
                warn!(
                    "Transcript {} has a 3 nucleotide offset in exon lengths likely due to the stop codon. Total exon length: {}, Transcript length: {}",
                    transcript.transcript_id, total_exon_length, transcript_length
                );
            } else {
                warn!(
                    "Transcript {} has inconsistent exon lengths. Total exon length: {}, Transcript length: {}",
                    transcript.transcript_id, total_exon_length, transcript_length
                );
            }
            transcript.has_missing_features = true;
        }
    }
}


// Unit tests 
#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;
    use env_logger;

    fn init() {
        let _ = env_logger::builder().is_test(true).try_init();
    }

    #[test]
    fn test_read_gtf_file() {
        init();
        let gtf_file = Path::new("./test/GRCm39_subset.gtf");
        let transcripts = read_annotation_file(gtf_file.to_str().unwrap(), true).unwrap();

        assert!(!transcripts.is_empty());
        // TODO: Add assertions to check the expected behavior and results
    }

    #[test]
    fn test_specific_transcript_info() {
        init(); 

        
        let gtf_file_path = "./test/gencode_v38.gtf";
        let target_transcript_id = "ENST00000579823";

        // Read the GTF file
        let transcripts = read_annotation_file(gtf_file_path, true).expect("Failed to read GTF file");
        
        // Check if the target transcript ID is present and print its details
        if let Some(transcript) = transcripts.get(target_transcript_id) {
            println!("Transcript ID: {}", transcript.transcript_id);
            println!("Gene ID: {:?}", transcript.gene_id);
            println!("Gene Name: {:?}", transcript.gene_name);
            println!("Chromosome: {}", transcript.chromosome);
            println!("Strand: {:?}", transcript.strand);
            println!("Exons Count: {}", transcript.exons.len());
            
        } else {
            println!("Transcript ID {} not found.", target_transcript_id);
        }
    }

    #[test]
    fn test_print_all_transcript_info() {
        init(); // Initialize logging if necessary

        // file path 
        let gtf_file_path = "./test/gencode_v38.gtf";

        // Read the GTF file
        let transcripts = read_annotation_file(gtf_file_path, true).expect("Failed to read GTF file");

        // Iterate over all transcripts and print their details
        for (transcript_id, transcript) in transcripts.iter() {
            println!("Transcript ID: {}", transcript_id);
            println!("Gene ID: {:?}", transcript.gene_id);
            println!("Gene Name: {:?}", transcript.gene_name);
            println!("Chromosome: {}", transcript.chromosome);
            println!("Strand: {:?}", transcript.strand);
            println!("Exons Count: {}", transcript.exons.len());
            println!("Biotype: {:?}", transcript.biotype);
            println!("CDS Length: {:?}", transcript.cds_len);
            println!("5' UTR Length: {:?}", transcript.utr5_len);
            println!("3' UTR Length: {:?}", transcript.utr3_len);
            println!("Has Missing Features: {}", transcript.has_missing_features);
            println!("Splice Junction Positions: {:?}", transcript.splice_junction_positions);

            // Print exon details for each transcript
            for exon in &transcript.exons {
                println!("\tExon: {:?}", exon);
            }

            println!("---------------------------------------------------");
        }
    }

    // // Test the parsing of GFF3 files 
    // #[test]
    // fn test_read_gff_file() {
    //     init();
        
    //     let gff_file = Path::new("./test/GRCm39_subset.gff3");
    //     let transcripts = read_annotation_file(gff_file.to_str().unwrap(), false).unwrap();
    //     assert!(!transcripts.is_empty());
    //     // TODO: Add assertions to check the expected behavior and results
    // }
}