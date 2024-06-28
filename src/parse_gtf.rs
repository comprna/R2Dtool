use std::collections::{HashMap, HashSet};
use multimap::MultiMap;
use bio::io::gff;
use rayon::prelude::*;
use std::io::{BufRead, BufReader};
use std::fs::File;
use std::io::Cursor;
use std::sync::Arc;
use log::warn;

// exon struct 
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

// transcript struct:
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

// default traits for transcripts 
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

pub fn parse_gff_attributes(attributes: &MultiMap<String, String>) -> HashMap<String, String> {
    attributes.iter().map(|(k, v)| (k.clone(), v.clone())).collect()
}
  
fn find_biotype(attributes: &HashMap<String, String>, keys: &[&str]) -> Option<String> {
    keys.iter().find_map(|key| attributes.get(*key)).cloned()
}

pub fn read_annotation_file(file_path: &str, is_gtf: bool, has_version: bool) -> Result<HashMap<String, Transcript>, Box<dyn std::error::Error>> {

    let file = File::open(file_path)?;
    let reader = BufReader::with_capacity(512 * 1024, file); // buffer size of 512KB
    let biotype_keys = Arc::new(vec!["transcript_biotype", "transcript_type", "gene_type", "gene_biotype"]);

    // process the file in chunks
    let chunk_size = 20000; // 20k lines per chunk
    let chunks: Vec<Vec<String>> = reader.lines()
        .collect::<Result<Vec<_>, _>>()?
        .chunks(chunk_size)
        .map(|chunk| chunk.to_vec())
        .collect();

    // process chunks in parallel 
    let results: Vec<(HashMap<String, Transcript>, HashMap<String, u32>, HashSet<String>)> = chunks.par_iter()
        .map(|chunk| process_chunk(chunk, is_gtf, has_version, &biotype_keys))
        .collect();

    // merge results 
    let mut transcripts = HashMap::new();
    let mut ignored_features = HashMap::new();
    let mut skipped_par_genes = HashSet::new();

    for (chunk_transcripts, chunk_ignored, chunk_skipped) in results {
        for (id, transcript) in chunk_transcripts {
            transcripts.entry(id)
                .and_modify(|existing: &mut Transcript| {
                    existing.exons.extend(transcript.exons.clone());
                    existing.cds_starts.extend(transcript.cds_starts.clone());
                    existing.cds_ends.extend(transcript.cds_ends.clone());
                    existing.cds_len = Some(existing.cds_len.unwrap_or(0) + transcript.cds_len.unwrap_or(0));
                    existing.utr5_len = Some(existing.utr5_len.unwrap_or(0) + transcript.utr5_len.unwrap_or(0));
                    existing.utr3_len = Some(existing.utr3_len.unwrap_or(0) + transcript.utr3_len.unwrap_or(0));
                })
                .or_insert(transcript);
        }
        for (feature, count) in chunk_ignored {
            *ignored_features.entry(feature).or_insert(0) += count;
        }
        skipped_par_genes.extend(chunk_skipped);
    }

    // final processing
    transcripts.par_iter_mut().for_each(|(_, transcript)| {
        let mut exon_features: Vec<&Exon> = transcript.exons.iter()
            .filter(|exon| exon.feature.as_ref().map_or(false, |ft| ft == "exon"))
            .collect();
        exon_features.sort_by_key(|exon| exon.start);

        if exon_features.len() > 1 {
            transcript.splice_junction_positions = exon_features.windows(2)
                .map(|window| window[0].end + 1)
                .collect();
        }

        transcript.transcript_length = Some(exon_features.iter().map(|exon| exon.length).sum());

        check_transcript_length(transcript);

        transcript.cds_starts.sort_unstable();
        transcript.cds_ends.sort_unstable();

        if let (Some(&cds_start), Some(&cds_end)) = (transcript.cds_starts.first(), transcript.cds_ends.last()) {
            for feature in &transcript.exons {
                if feature.feature_type == "UTR" {
                    let strand = transcript.strand.as_deref().unwrap_or(".");
                    match strand {
                        "+" | "." => {
                            if feature.end < cds_start {
                                transcript.utr5_len = Some(transcript.utr5_len.unwrap_or(0) + feature.length);
                            } else if feature.start > cds_end {
                                transcript.utr3_len = Some(transcript.utr3_len.unwrap_or(0) + feature.length);
                            }
                        },
                        "-" => {
                            if feature.start > cds_end {
                                transcript.utr5_len = Some(transcript.utr5_len.unwrap_or(0) + feature.length);
                            } else if feature.end < cds_start {
                                transcript.utr3_len = Some(transcript.utr3_len.unwrap_or(0) + feature.length);
                            }
                        },
                        _ => {}
                    }
                }
            }
        }
    });

    // report ignored feature types
    if !ignored_features.is_empty() {
        warn!("Ignored the following annotation features:");
        for (feature_type, count) in ignored_features {
            warn!("- {}: {} occurrences", feature_type, count);
        }
    }

    // report skipped PAR genes 
    if !skipped_par_genes.is_empty() {
        eprintln!("Skipped {} unique genes with '_PAR_' in their identifiers while parsing GTF", skipped_par_genes.len());
    }
    
    Ok(transcripts)
}

// process each chunk 
fn process_chunk(
    chunk: &[String],
    is_gtf: bool,
    has_version: bool,
    biotype_keys: &Arc<Vec<&str>>
) -> (HashMap<String, Transcript>, HashMap<String, u32>, HashSet<String>) {

    let mut transcripts = HashMap::new();
    let mut ignored_features = HashMap::new();
    let mut skipped_par_genes = HashSet::new();

    // convert chunk into a single string
    let chunk_data = chunk.join("\n");
    let cursor = Cursor::new(chunk_data);

    // create GFF reader from cursor
    let gff_type = if is_gtf { gff::GffType::GTF2 } else { gff::GffType::GFF3 };
    let mut reader = gff::Reader::new(cursor, gff_type);

    for record_result in reader.records() {
        let record = match record_result {
            Ok(record) => record,
            Err(_) => continue,  
        };

        let attributes = parse_gff_attributes(record.attributes());
        let transcript_id_attr = if is_gtf { "transcript_id" } else { "ID" };
        let gene_id_attr = if is_gtf { "gene_id" } else { "Parent" };

        // skip PARs and record 
        if let Some(gene_id) = attributes.get(gene_id_attr) {
            if gene_id.contains("_PAR_") {
                skipped_par_genes.insert(gene_id.clone()); 
                continue;
            }
        }

        let transcript_id_with_version = match attributes.get(transcript_id_attr) {
            Some(id) => id,
            None => continue,
        };

        let transcript_id = if has_version {
            transcript_id_with_version.to_string()
        } else {
            match transcript_id_with_version.split('.').next() {
                Some(id) => id.to_string(),
                None => continue,
            }
        };

        if *record.start() > *record.end() {
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
                    length: feature_length,
                    score: record.score().map(|s| s as f64),
                    strand: record.strand().map(|s| s.strand_symbol().to_string()).unwrap_or_else(|| ".".to_string()),
                    frame: record.frame().chars().next().filter(|&c| c == '0' || c == '1' || c == '2' || c == '.'),
                    attributes: attributes.clone(),
                    feature: Some(record.feature_type().to_string()),
                };

                transcripts.entry(transcript_id.clone())
                    .and_modify(|existing: &mut Transcript| {
                        if existing.chromosome != exon.seq_id {
                            existing.has_missing_features = true;
                        }
                        existing.exons.push(exon.clone());
                    })
                    .or_insert_with(|| {
                        let biotype = find_biotype(&attributes, &biotype_keys);
                        let strand = record.strand().map(|s| s.strand_symbol().to_string()).unwrap_or_else(|| ".".to_string());
                        let chromosome = record.seqname().to_string();
                        let gene_id = attributes.get(gene_id_attr).cloned();
                        let gene_name = attributes.get("gene_name").cloned();
                        Transcript {
                            transcript_id: transcript_id.clone(),
                            gene_id,
                            gene_name,
                            utr5_len: None,
                            cds_len: None,
                            utr3_len: None,
                            has_missing_features: false,
                            exons: vec![exon.clone()],
                            biotype,
                            splice_junction_positions: Vec::new(),
                            strand: Some(strand),
                            chromosome,
                            cds_starts: Vec::new(),
                            cds_ends: Vec::new(),
                            transcript_length: None,
                        }
                    });
            }
            "CDS" => {
                if let Some(transcript) = transcripts.get_mut(&transcript_id) {
                    transcript.cds_len = Some(transcript.cds_len.unwrap_or(0) + feature_length);
                    transcript.cds_starts.push(*record.start());
                    transcript.cds_ends.push(*record.end());
                }
            }
            "five_prime_utr" | "5UTR" => {
                if let Some(transcript) = transcripts.get_mut(&transcript_id) {
                    transcript.utr5_len = Some(transcript.utr5_len.unwrap_or(0) + feature_length);
                }
            }
            "three_prime_utr" | "3UTR" => {
                if let Some(transcript) = transcripts.get_mut(&transcript_id) {
                    transcript.utr3_len = Some(transcript.utr3_len.unwrap_or(0) + feature_length);
                }
            }
            "transcript" | "mRNA" => {
                // TODO: parse transcript/mRNA features 
            }
            "UTR" => {
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
                    feature: Some("UTR".to_string()),
                };

                transcripts.entry(transcript_id.clone())
                    .and_modify(|existing: &mut Transcript| {
                        existing.exons.push(utr_feature.clone());
                    })
                    .or_insert_with(|| {
                        let biotype = find_biotype(&attributes, &biotype_keys);
                        let strand = record.strand().map(|s| s.strand_symbol().to_string()).unwrap_or_else(|| ".".to_string());
                        let chromosome = record.seqname().to_string();
                        let gene_id = attributes.get(gene_id_attr).cloned();
                        let gene_name = attributes.get("gene_name").cloned();
                        Transcript {
                            transcript_id: transcript_id.clone(),
                            gene_id,
                            gene_name,
                            utr5_len: None,
                            cds_len: None,
                            utr3_len: None,
                            has_missing_features: false,
                            exons: vec![utr_feature.clone()],
                            biotype,
                            splice_junction_positions: Vec::new(),
                            strand: Some(strand),
                            chromosome,
                            cds_starts: Vec::new(),
                            cds_ends: Vec::new(),
                            transcript_length: None,
                        }
                    });
            }
            other => {
                *ignored_features.entry(other.to_string()).or_insert(0) += 1;
            }
        }
    }

    (transcripts, ignored_features, skipped_par_genes)
}

// ensure UTR and CDS lengths are consistent with exon lengths
fn check_transcript_length(transcript: &mut Transcript) {
    if transcript.cds_len.unwrap_or(0) > 0 || transcript.utr5_len.unwrap_or(0) > 0 || transcript.utr3_len.unwrap_or(0) > 0 {
        let total_exon_length: u64 = transcript.exons.iter().map(|exon| exon.length).sum();

        let transcript_length = transcript.utr5_len.unwrap_or(0)
            + transcript.cds_len.unwrap_or(0)
            + if transcript.cds_len.unwrap_or(0) > 0 { 3 } else { 0 }
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
        let gtf_file = Path::new("./test/GRCh38.110_subset.gtf");
        let transcripts = read_annotation_file(gtf_file.to_str().unwrap(), true, false).unwrap();

        assert!(!transcripts.is_empty());
        // TODO: Add assertions to check the expected behavior and results
    }

    #[test]
    fn test_specific_transcript_info() {
        init(); 

        
        let gtf_file_path = "./test/GRCh38.110_subset.gtf";
        let target_transcript_id = "ENST00000400109";

        // Read the GTF file
        let transcripts = read_annotation_file(gtf_file_path, true, false).expect("Failed to read GTF file");
        
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
        let gtf_file_path = "./test/GRCh38.110_subset.gtf";

        // Read the GTF file
        let transcripts = read_annotation_file(gtf_file_path, true, false).expect("Failed to read GTF file");

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
}