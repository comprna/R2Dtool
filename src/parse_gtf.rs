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
                        };
                        entry.insert(transcript);
                        debug!("Created new transcript: {}", transcript_id);
                    }
                }
            }
            "CDS" => {
                if let Some(transcript) = transcripts.get_mut(&transcript_id) {
                    transcript.cds_len = Some(transcript.cds_len.unwrap_or(0) + feature_length);
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
            other => {
                *ignored_features.entry(other.to_string()).or_insert(0) += 1;
            }
        }
    }

    debug!("Number of transcripts parsed: {}", transcripts.len());

    for transcript in transcripts.values_mut() {
        transcript.exons.sort_by_key(|exon| exon.start);

        if transcript.exons.len() > 1 {
            for i in 0..(transcript.exons.len() - 1) {
                let splice_pos = transcript.exons[i].end + 1;
                transcript.splice_junction_positions.push(splice_pos);
            }
        }

        check_transcript_length(transcript);
    }

    if !ignored_features.is_empty() {
        warn!("Ignored the following annotation features:");
        for (feature_type, count) in ignored_features {
            warn!("- {}: {} occurrences", feature_type, count);
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