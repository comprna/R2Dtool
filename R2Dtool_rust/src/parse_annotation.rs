use std::collections::HashMap;
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

// Define transcript structure, which contains exons
#[derive(Debug, PartialEq, Clone)]
pub struct Transcript {
    pub transcript_id: String,
    pub utr5_len: Option<u64>,
    pub cds_len: Option<u64>,
    pub utr3_len: Option<u64>,
    pub has_missing_features: bool,
    pub exons: Vec<Exon>,
    pub biotype: Option<String>,
    pub splice_junction_positions: Vec<u64>,
    pub strand: Option<String>,
    pub chromosome: String, // Add this line
}


// Default trait for transript structure
impl Default for Transcript {
    fn default() -> Self {
        Self {
            transcript_id: String::new(),
            utr5_len: None,
            cds_len: None,
            utr3_len: None,
            has_missing_features: false,
            exons: Vec::new(),
            biotype: None,
            splice_junction_positions: Vec::new(),
            strand: None,
            chromosome: String::new(), // Add this line
        }
    }
}

// Read the gtf file
// create a HashMap<String, Transcript>, each key (a String) will map to a single Transcript object.
// each transcript has up to several exon objects stored in the exons field, which is a vector of Exon structs.
pub fn read_gtf_file(gtf_file: &str) -> HashMap<String, Transcript> {

    // create a empty hashmap, transcripts
    let mut transcripts: HashMap<String, Transcript> = HashMap::new();

    // initialise the reader
    let mut reader = gff::Reader::from_file(gtf_file, gff::GffType::GTF2).unwrap();

    // define possible keys of the biotype; includes  Ensembl, gencode, flybase, pombase, (others?)
    let biotype_keys = vec!["gene_biotype", "gene_type", "transcript_biotype"];

    // slice the biotype, if it exists
    fn find_biotype(attributes: &HashMap<String, String>, keys: &[&str]) -> Option<String> {
        for key in keys {
            if let Some(value) = attributes.get(*key) {
                return Some(value.to_string());
            }
        }
        None
    }

    // loop over the records of the gtf file
    for record in reader.records() {

        // print lines to Debug
        // println!("Record found: {:?}", record);

        // unwrap the record and store in 'record' object
        let record = record.unwrap();

        // skip header records if present
        if record.feature_type() == "#" {
            continue;
        }

        // parse attributes
        let attributes = parse_gff_attributes(record.attributes());
        if let Some(transcript_id_with_version) = attributes.get("transcript_id") {
            let transcript_id = transcript_id_with_version.split('.').next().unwrap(); // Remove version from transcript ID
            let transcript_id = transcript_id.to_string();

            let transcript = transcripts.entry(transcript_id.clone()).or_insert_with(|| {
                let biotype = find_biotype(&attributes, &biotype_keys);
                let strand = record.strand().expect("Missing strand information").strand_symbol().to_string(); // Get strand for transcript
                let chromosome = record.seqname().to_string(); // Define chromosome here <-- ADD THIS LINE
                Transcript {
                    transcript_id: transcript_id.clone(),
                    utr5_len: None,
                    cds_len: None,
                    utr3_len: None,
                    has_missing_features: false,
                    exons: Vec::new(),
                    biotype,
                    splice_junction_positions: Vec::new(),
                    strand: Some(strand),
                    chromosome, // Add this line
                }
            });

            // calculate feature length
            let feature_length = (*record.end() - *record.start() + 1) as u64;

            // store exon properties, if an exon is found
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
                _ => (),
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

    println!("Size of the transcripts hashmap: {}", transcripts.len()); // Add this line
    transcripts
}

// read gff file
pub fn read_gff_file(gff_file: &str) -> HashMap<String, Transcript> {
    let mut transcripts: HashMap<String, Transcript> = HashMap::new();
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
                attributes: attributes,
                feature: Some(feature_type.to_string()),
                length: (*record.end() - *record.start() + 1) as u64,
            };
            if let Some(transcript_id_with_version) = exon.attributes.get("transcript_id") {
                let transcript_id = transcript_id_with_version.split('.').next().unwrap(); // Remove version from transcript ID
                let transcript_id = transcript_id.to_string();

                let transcript = transcripts.entry(transcript_id.clone()).or_insert_with(|| {
                    let biotype = find_biotype(&exon.attributes, &biotype_keys);
                    let strand = record.strand().expect("Missing strand information").strand_symbol().to_string(); // Get strand for transcript
                    let chromosome = record.seqname().to_string();
                    Transcript {
                        transcript_id: transcript_id.clone(),
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

    println!("Size of the transcripts hashmap: {}", transcripts.len()); // Add this line
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

// pub fn multimap_to_hashmap(multimap: &MultiMap<String, String>) -> HashMap<String, String> {
//     let mut hashmap = HashMap::new();
//     for (key, value) in multimap.iter() {
//         hashmap.insert(key.clone(), value.clone());
//     }
//     hashmap
// }
//



// pub fn parse_annotations(annotations: &MultiMap<String, Exon>) -> HashMap<String, Transcript> {
//     let mut transcripts: HashMap<String, Transcript> = HashMap::new();
//
//     for (transcript_id, exons) in annotations.iter_all() {
//         let transcript_id_no_version = transcript_id.split('.').next().unwrap().to_string();
//         let exons: Vec<Exon> = exons.clone().into_iter().collect();
//
//         if let Some(transcript) = transcripts.get_mut(&transcript_id_no_version) {
//             transcript.exons = exons;
//         } else {
//             let transcript = Transcript {
//                 transcript_id: transcript_id.clone(),
//                 utr5_len: None,
//                 cds_len: None,
//                 utr3_len: None,
//                 has_missing_features: false,
//                 exons,
//                 biotype: None, // Add this line
//                 splice_junction_positions: Vec::new(), // Add this line
//
//             };
//             transcripts.insert(transcript_id_no_version.clone(), transcript);
//         }
//     }
//
//     // println!("Number of exons per transcript:");
//     // for (transcript_id, exons) in annotations {
//     //     println!("{}: {} exons", transcript_id, exons.len());
//     // }
//
//     transcripts
// }

// pub fn parse_annotations(annotations: &MultiMap<String, Exon>) -> HashMap<String, Transcript> {
//     let mut transcripts: HashMap<String, Transcript> = HashMap::new();
//
//     for (transcript_id, exon) in annotations.iter() {
//         let transcript = transcripts.entry(transcript_id.to_string()).or_insert_with(|| Transcript {
//             transcript_id: transcript_id.clone(),
//             utr5_len: None,
//             cds_len: None,
//             utr3_len: None,
//             has_missing_features: false,
//             exons: Vec::new(),
//         });
//
//         transcript.exons.push(exon.to_owned());
//     }
//
//     // Check if any of the features are missing
//     for (_, transcript) in transcripts.iter_mut() {
//         if transcript.utr5_len.is_none() || transcript.cds_len.is_none() || transcript.utr3_len.is_none() {
//             transcript.has_missing_features = true;
//         }
//     }
//
//     // Print the number of exons per transcript
//     println!("Number of exons per transcript:");
//     for (transcript_id, exons) in annotations {
//         println!("{}: {} exons", transcript_id, exons.len());
//     }
//
//
//     transcripts
// }



// pub fn preview_annotations(annotations: &MultiMap<String, Exon>) {
//     println!("Previewing the first transcript:");
//
//     for (key, exons) in annotations.iter_all() {
//         println!("Transcript ID: {}", key);
//         for exon in exons {
//             println!("  Exon: {:?}", exon);
//         }
//         // break; // Add this line to exit the loop after the first transcript
//     }
// }

// pub fn print_exon_info(annotations: &MultiMap<String, Exon>) {
//     // Print table headers
//     println!(
//         "{:<10} {:<10} {:<15} {:<15} {:<7}",
//         "Start", "End", "Transcript ID", "Gene ID", "Strand"
//     );
//     println!("{}", "-".repeat(57));
//
//     let na_gene_id = String::from("N/A");
//
//     for (transcript_id, exon) in annotations.iter() {
//         let gene_id_str = exon.attributes.get("gene_id").unwrap_or(&na_gene_id);
//
//         println!(
//             "{:<10} {:<10} {:<15} {:<15} {:<7}",
//             exon.start, exon.end, transcript_id, gene_id_str, exon.strand
//         );
//     }
// }
