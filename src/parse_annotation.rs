use std::collections::{HashMap, HashSet};
use multimap::MultiMap;
use bio::io::gff;

// Define exon structure
#[derive(Debug, PartialEq, Clone)] // Add Clone trait
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
// Each transcript maps to a agene 
// And contains one or more exons
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

// Default traits for transcript 
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

// Read the gtf file; create a HashMap<String, Transcript>, each key (a String) will map to a single Transcript object.
// each transcript has up to several exon objects stored in the exons field, which is a vector of Exon structs.
pub fn read_gtf_file(gtf_file: &str) -> HashMap<String, Transcript> {
    
    
    let mut transcripts: HashMap<String, Transcript> = HashMap::new();

    // Count unused features 
    let mut ignored_features: HashMap<String, u32> = HashMap::new();


    let mut reader = gff::Reader::from_file(gtf_file, gff::GffType::GTF2)
        .expect("Unable to open the gtf file");

    let biotype_keys = vec!["transcript_biotype", "transcript_type", "gene_type"];

    fn find_biotype(attributes: &HashMap<String, String>, keys: &[&str]) -> Option<String> {
        for key in keys {
            if let Some(value) = attributes.get(*key) {
                return Some(value.to_string());
            }
        }
        None
    }

    let mut encountered_transcripts: HashSet<String> = HashSet::new();

    for record in reader.records() {
        let record = match record {
            Ok(record) => record,
            Err(err) => {
                eprintln!("Error reading record: {}", err);
                continue;
            }
        };

        if record.feature_type() == "#" {
            continue;
        }

        let attributes = parse_gff_attributes(record.attributes());

        let transcript_id_with_version = match attributes.get("transcript_id") {
            Some(id) => id,
            None => {
                eprintln!(
                    "Error: Missing 'transcript_id' in record: {:?}. Skipping...",
                    record
                );
                continue;
            }
        };
        let transcript_id = transcript_id_with_version.split('.').next().unwrap().to_string();

        let transcript = transcripts.entry(transcript_id.clone()).or_insert_with(|| { // missing brace here
            let biotype = find_biotype(&attributes, &biotype_keys);
            let strand = record.strand().expect("Missing strand information").strand_symbol().to_string();
            let chromosome = record.seqname().to_string();
            let gene_id = attributes.get("gene_id").map(|s| s.to_string());
            let gene_name = attributes.get("gene_name").map(|s| s.to_string());
            Transcript {
                transcript_id: transcript_id.clone(),
                gene_id,
                gene_name,
                utr5_len: None,
                cds_len: None,
                utr3_len: None,
                has_missing_features: false,
                exons: Vec::new(),
                biotype,
                splice_junction_positions: Vec::new(),
                strand: Some(strand),
                chromosome,
            }
        }); 

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
                    strand: record.strand().expect("Missing strand information").strand_symbol().to_string(),
                    frame: record.frame().chars().next(),
                    attributes: attributes.clone(),
                    feature: Some(record.feature_type().to_string()),
                    length: feature_length,
                };
                transcript.exons.push(exon);
            }
            "CDS" => {
                transcript.cds_len = Some(transcript.cds_len.unwrap_or(0) + feature_length);
            }
            "five_prime_utr" | "5UTR" => {
                transcript.utr5_len = Some(transcript.utr5_len.unwrap_or(0) + feature_length);
            }
            "three_prime_utr" | "3UTR" => {
                transcript.utr3_len = Some(transcript.utr3_len.unwrap_or(0) + feature_length);
            }
            "transcript" => {
                if !encountered_transcripts.insert(transcript_id.clone()) {
                    eprintln!(
                        "Warning: More than one transcript with the same name '{}' found. Please check the GTF file.",
                        transcript_id
                    );
                }
            }
            other => {
                *ignored_features.entry(other.to_string()).or_insert(0) += 1;
            }
            }
        
    }

    transcripts.retain(|_, transcript| {
        let exon_seq_id = transcript.exons.get(0).map(|exon| &exon.seq_id);
        if let Some(exon_seq_id) = exon_seq_id {
            transcript.exons.iter().all(|exon| exon.seq_id == *exon_seq_id)
        } else {
            false
        }
    });

    for transcript in transcripts.values_mut() {
        transcript.exons.sort_by_key(|exon| exon.start);

        if transcript.exons.len() > 1 {
            for i in 0..(transcript.exons.len() - 1) {
                let splice_pos = transcript.exons[i].end + 1;
                transcript.splice_junction_positions.push(splice_pos);
            }
        }
    }

    if !ignored_features.is_empty() {
        eprintln!("Ignored the following annotation features:");
        for (feature_type, count) in ignored_features {
            eprintln!("- {}: {} occurrences", feature_type, count);
        }
    }

    transcripts
}

// read gff file
pub fn read_gff_file(gff_file: &str) -> HashMap<String, Transcript> {

    let mut transcripts: HashMap<String, Transcript> = HashMap::new();

    // Count unused features 
    // let mut ignored_features: HashMap<String, u32> = HashMap::new();


    let mut reader = gff::Reader::from_file(gff_file, gff::GffType::GFF3).unwrap();

    let biotype_keys = vec!["gene_biotype", "gene_type", "transcript_biotype"];

    fn find_biotype(attributes: &HashMap<String, String>, keys: &[&str]) -> Option<String> {
        for key in keys {
            if let Some(value) = attributes.get(*key) {
                return Some(value.to_string());
            }
        }
        None
    }

    for record in reader.records() {

        // println!("GFF Record found: {:?}", record);

        let record = record.unwrap();

        // Handle exon, CDS, 5' UTR, and 3' UTR features with both naming conventions
        if ["exon", "CDS", "five_prime_utr", "5UTR", "three_prime_utr", "3UTR"].contains(&record.feature_type()) {
            let feature_type = match record.feature_type() {
                "five_prime_utr" | "5UTR" => "5UTR",
                "three_prime_utr" | "3UTR" => "3UTR",
                _ => record.feature_type(),
            };

            let attributes = parse_gff_attributes(record.attributes());

            let exon = Exon {
                seq_id: record.seqname().to_string(),
                source: record.source().to_string(),
                feature_type: feature_type.to_string(),
                start: *record.start(),
                end: *record.end(),
                score: record.score().map(|s| s as f64),
                strand: record.strand().expect("Missing strand information").strand_symbol().to_string(),
                frame: record.frame().chars().next(),
                attributes: attributes.clone(),
                feature: Some(feature_type.to_string()),
                length: (*record.end() - *record.start() + 1) as u64,
            };

            if let Some(transcript_id_with_version) = exon.attributes.get("transcript_id") {
                let transcript_id = transcript_id_with_version.split('.').next().unwrap(); // Remove version from transcript ID
                let transcript_id = transcript_id.to_string();

                // Safety check to ensure that all exons are on the same chromosome (PAR_Y case)
                transcripts.retain(|_, transcript| {
                    if transcript.chromosome == exon.seq_id {
                        transcript.exons.push(exon.clone());
                        true
                    } else {
                        eprintln!(
                            "Warning: Transcript {} has exons on different chromosomes. Removing transcript from object.",
                            transcript.transcript_id
                        );
                        false
                    }
                });

                let transcript = transcripts.entry(transcript_id.clone()).or_insert_with(|| {
                    let biotype = find_biotype(&exon.attributes, &biotype_keys);
                    let strand = record.strand().expect("Missing strand information").strand_symbol().to_string(); // Get strand for transcript
                    let chromosome = record.seqname().to_string();
                    let gene_id = attributes.get("gene_id").map(|s| s.to_string());
                    let gene_name = attributes.get("gene_name").map(|s| s.to_string());
                    Transcript {
                        transcript_id: transcript_id.clone(),
                        gene_id,
                        gene_name,
                        utr5_len: None,
                        cds_len: None,
                        utr3_len: None,
                        has_missing_features: false,
                        exons: Vec::new(),
                        biotype,
                        splice_junction_positions: Vec::new(),
                        strand: Some(strand),
                        chromosome,
                    }
                });

                transcript.exons.push(exon);
            }
        }
    }

    // Update splice_junction_positions for each transcript
    for transcript in transcripts.values_mut() {
        transcript.exons.sort_by_key(|exon| exon.start);

        if transcript.exons.len() > 1 {
            for i in 0..(transcript.exons.len() - 1) {
                let splice_pos = transcript.exons[i].end + 1;
                transcript.splice_junction_positions.push(splice_pos);
            }
        }
    }

    println!("Size of the transcripts hashmap: {}", transcripts.len()); 
    transcripts

}

// parse gff attributes
pub fn parse_gff_attributes(attributes: &MultiMap<String, String>) -> HashMap<String, String> {
    let mut attr_map = HashMap::new();
    for (key, value) in attributes.iter() {
        attr_map.insert(key.clone(), value.clone());
    }
    attr_map
}

// safety code

// check transcript length
// fn check_transcript_length(transcript: &Transcript) {
//     let mut total_exon_length: u64 = 0;
//     for exon in &transcript.exons {
//         total_exon_length += exon.length;
//     }
//
//     let transcript_length = transcript.utr5_len.unwrap_or(0)
//     + transcript.cds_len.unwrap_or(0)
//     + transcript.utr3_len.unwrap_or(0) + 3;
//
//     if total_exon_length != transcript_length {
//         eprintln!(
//             "Warning: Transcript {} has inconsistent exon lengths. Total exon length: {}, Transcript length: {}",
//             transcript.transcript_id, total_exon_length, transcript_length
//         );
//     }
// }
