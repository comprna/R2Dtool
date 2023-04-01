use std::collections::HashMap;
use multimap::MultiMap;
use bio::io::gff;

// define exon structure
#[derive(Debug, PartialEq)]
pub struct Exon {
    pub seq_id: String,
    pub source: String,
    pub feature_type: String,
    pub start: u64,
    pub end: u64,
    pub score: Option<f64>,
    pub strand: String,
    pub frame: Option<char>,
    pub attributes: HashMap<String, String>,
}

// read the gtf file
pub fn read_gtf_file(gtf_file: &str) -> MultiMap<String,Exon> {
    let mut annotations = MultiMap::new();
    let mut reader = gff::Reader::from_file(gtf_file, gff::GffType::GTF2).unwrap();
    for record in reader.records() {
        let record = record.unwrap();

        if record.feature_type() == "#" {
            continue;
        }

        if record.feature_type() == "exon" {
            let attributes = parse_gff_attributes(record.attributes());

            let exon = Exon {
                seq_id: record.seqname().to_string(),
                source: record.source().to_string(),
                feature_type: record.feature_type().to_string(),
                start: *record.start(),
                end: *record.end(),
                score: record.score().map(|s| s as f64),
                strand: record.strand().expect("Missing strand information").strand_symbol().to_string(),
                frame: record.frame().chars().next(),
                attributes: attributes,
            };
            if let Some(transcript_id) = exon.attributes.get("transcript_id") {
                annotations.insert(transcript_id.to_string(), exon);
            }
        }
    }

    annotations
}

pub fn read_gff_file(gff_file: &str) -> MultiMap<String, Exon> {
    let mut annotations = MultiMap::new();
    let mut reader = gff::Reader::from_file(gff_file, gff::GffType::GFF3).unwrap();

    for record in reader.records() {
        let record = record.unwrap();
        let attributes = record.attributes().clone();

        if record.feature_type() == "exon" {
            let exon = Exon {
                seq_id: record.seqname().to_string(),
                source: record.source().to_string(),
                feature_type: record.feature_type().to_string(),
                start: *record.start(),
                end: *record.end(),
                score: record.score().map(|s| s as f64),
                strand: record.strand().expect("Missing strand information").strand_symbol().to_string(),
                frame: record.frame().chars().next(),
                attributes: multimap_to_hashmap(&attributes),
            };
            if let Some(transcript_id) = exon.attributes.get("Parent") {
                annotations.insert(transcript_id.to_string(), exon);
            }
        }
    }

    annotations
}

pub fn multimap_to_hashmap(multimap: &MultiMap<String, String>) -> HashMap<String, String> {
    let mut hashmap = HashMap::new();
    for (key, value) in multimap.iter() {
        hashmap.insert(key.clone(), value.clone());
    }
    hashmap
}

// parse gff attributes
pub fn parse_gff_attributes(attributes: &MultiMap<String, String>) -> HashMap<String, String> {
    let mut attr_map = HashMap::new();
    for (key, value) in attributes.iter() {
        attr_map.insert(key.clone(), value.clone());
    }
    attr_map
}

// pub fn preview_annotations(annotations: &MultiMap<String, Exon>) {
//     println!("Previewing the first transcript:");
//
//     for (key, exons) in annotations.iter_all() {
//         println!("Transcript ID: {}", key);
//         // for exon in exons {
//         //     println!("  Exon: {:?}", exon);
//         // }
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
