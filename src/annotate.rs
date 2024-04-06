use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::error::Error;
use std::collections::HashMap;
use crate::parse_gtf::{Transcript, read_annotation_file};


#[derive(Debug, Clone)]
pub struct SpliceSite {
    pub transcript_id: String,
    pub tx_coord: u64,
}

// Additional functions to calculate missing features
fn calculate_tx_len(utr5_len: u64, cds_len: u64, utr3_len: u64) -> u64 {
    utr5_len + cds_len + utr3_len + 3
}

fn calculate_cds_start(utr5_len: u64) -> u64 {
    utr5_len
}

fn calculate_cds_end(utr5_len: u64, cds_len: u64) -> u64 {
    utr5_len + cds_len
}

fn calculate_tx_end(utr5_len: u64, cds_len: u64, utr3_len: u64) -> u64 {
    utr5_len + cds_len + utr3_len
}

fn calculate_meta_coordinates(tx_coord: u64, utr5_len: u64, cds_len: u64, utr3_len: u64) -> (f64, i64, i64) {
    let cds_start = utr5_len;
    let cds_end = utr5_len + cds_len;
    let _tx_end = cds_end + utr3_len;

    let rel_pos = if tx_coord < cds_start {
        tx_coord as f64 / cds_start as f64
    } else if tx_coord < cds_end {
        1.0 + (tx_coord - utr5_len) as f64 / cds_len as f64
    } else {
        2.0 + (tx_coord - utr5_len - cds_len) as f64 / utr3_len as f64
    };

    let abs_cds_start = tx_coord as i64 - cds_start as i64;
    let abs_cds_end = tx_coord as i64 - cds_end as i64;

    (rel_pos, abs_cds_start, abs_cds_end)
}

fn generate_splice_sites(transcripts: &HashMap<String, Transcript>) -> HashMap<String, Vec<SpliceSite>> {
    let mut splice_sites_map: HashMap<String, Vec<SpliceSite>> = HashMap::new();

    for (transcript_id, transcript) in transcripts {
        let transcript_id_no_version = transcript_id.split('.').next().unwrap().to_string();
        let mut splice_sites: Vec<SpliceSite> = Vec::new();

        if !transcript.exons.is_empty() {
            let ref exons = transcript.exons;
            let mut cumsum_width = 0;

            for (i, exon) in exons.iter().enumerate() {
                let exon_length = exon.end - exon.start + 1;
                cumsum_width += exon_length;

                if i < exons.len() - 1 {
                    splice_sites.push(SpliceSite {
                        transcript_id: transcript_id_no_version.clone(),
                        tx_coord: cumsum_width as u64,
                    });
                }
            }
        }
        splice_sites_map.insert(transcript_id_no_version, splice_sites);
    }

    splice_sites_map
}

fn splice_site_distances(tx_coord: u64, splice_sites: &[SpliceSite]) -> (Option<i64>, Option<i64>) {
    let mut upstream_distance: Option<i64> = None;
    let mut downstream_distance: Option<i64> = None;

    for splice_site in splice_sites {
        let site_distance = tx_coord as i64 - splice_site.tx_coord as i64;

        if site_distance >= 0 {
            if downstream_distance.is_none() || site_distance < downstream_distance.unwrap() {
                downstream_distance = Some(site_distance);
            }
        } else {
            if upstream_distance.is_none() || site_distance.abs() < upstream_distance.unwrap() {
                upstream_distance = Some(site_distance.abs());
            }
        }
    }

    (upstream_distance, downstream_distance)
}

pub fn run_annotate(matches: &clap::ArgMatches, has_header: bool) -> Result<(), Box<dyn Error>> {
   
    eprintln!("Running the annotate functionality...");
   
    let gtf_file: String = matches.get_one::<String>("gtf").unwrap().to_string();
    let input_file: String = matches.get_one::<String>("input").unwrap().to_string();
    let output_file: Option<String> = matches.get_one::<String>("output").map(|s: &String| s.to_string());

    // By default, read in the annotations as GTF file
    // TODO: implement GFF3 parsing         
    // let default_format = String::from("gtf");
    // let format = matches.get_one("format").unwrap_or(&default_format);
    let annotations = read_annotation_file(&gtf_file, true)?;
    
    // Print the annotations in a table
    // eprintln!("Previewing transcript annotations\n");
    // preview_annotations(&annotations);

    let transcripts = annotations;

    let mut input_reader = BufReader::new(File::open(input_file.clone()).unwrap_or_else(|_| panic!("Cannot open input file: {}", input_file)));

    let mut output_writer: Box<dyn Write> = match output_file {
        Some(file) => Box::new(File::create(file).expect("Cannot create output file")),
        None => Box::new(std::io::stdout()),
    };

    let splice_sites = generate_splice_sites(&transcripts);

    let mut header = String::new();
    if has_header {
        input_reader.read_line(&mut header).unwrap();
        let header_fields: Vec<&str> = header.trim().split('\t').collect();

        let output_header = format!(
    "{}\tgene_id\tgene_name\ttranscript_biotype\ttx_len\tcds_start\tcds_end\ttx_end\trel_pos\tabs_cds_start\tabs_cds_end\tup_junc_dist\tdown_junc_dist",
    header_fields.join("\t")
);
        writeln!(output_writer, "{}", output_header).unwrap();
    }

    for line in input_reader.lines() {
        let line = line.unwrap();
        let fields: Vec<&str> = line.split('\t').collect();
        let transcript_id_with_version = fields[0];
        let transcript_id = transcript_id_with_version.split('.').next().unwrap();
        let tx_coord: u64 = fields[1].parse().unwrap();
    
        if let Some(transcript) = transcripts.get(transcript_id) {
            
            // Initialize all fields to "NA"
            let mut tx_len = "NA".to_string();
            let mut cds_start = "NA".to_string();
            let mut cds_end = "NA".to_string();
            let mut tx_end = "NA".to_string();
            let mut rel_pos = "NA".to_string();
            let mut abs_cds_start = "NA".to_string();
            let mut abs_cds_end = "NA".to_string();
            let mut up_junc_dist = "NA".to_string();
            let mut down_junc_dist = "NA".to_string(); 
        
            // Always populate gene_id, gene_name, and biotype if available
            let gene_id = transcript.gene_id.clone().unwrap_or_else(|| "NA".to_string());
            let gene_name = transcript.gene_name.clone().unwrap_or_else(|| "NA".to_string());
            let biotype = transcript.biotype.clone().unwrap_or_else(|| "NA".to_string());
    
            // Calculate transcript length if possible
            if let (Some(utr5_len), Some(cds_len), Some(utr3_len)) = (transcript.utr5_len, transcript.cds_len, transcript.utr3_len) {
                tx_len = calculate_tx_len(utr5_len, cds_len, utr3_len).to_string();
                cds_start = calculate_cds_start(utr5_len).to_string();
                cds_end = calculate_cds_end(utr5_len, cds_len).to_string();
                tx_end = calculate_tx_end(utr5_len, cds_len, utr3_len).to_string();
                let calculated_values = calculate_meta_coordinates(tx_coord, utr5_len, cds_len, utr3_len);
                rel_pos = format!("{:.5}", calculated_values.0);
                abs_cds_start = calculated_values.1.to_string();
                abs_cds_end = calculated_values.2.to_string();
            }
    
            // Handle splice sites if available
            if let Some(splice_sites) = splice_sites.get(transcript_id) {
                let calculated_distances = splice_site_distances(tx_coord, splice_sites);
                up_junc_dist = calculated_distances.0.map_or("NA".to_string(), |x| x.to_string());
                down_junc_dist = calculated_distances.1.map_or("NA".to_string(), |x| x.to_string());
            }
    
            // Construct and write the output line
            let output_line = format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                line,
                gene_id,
                gene_name,
                biotype,
                tx_len,
                cds_start,
                cds_end,
                tx_end,
                rel_pos,
                abs_cds_start,
                abs_cds_end,
                up_junc_dist,
                down_junc_dist
            );
            if let Err(_) = writeln!(output_writer, "{}", output_line) {
                break;
            }
        } else {
            // Handle the case where no transcript data is found
            if let Err(_) = writeln!(output_writer, "{}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA", line) {
                break;
            }
        }
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
