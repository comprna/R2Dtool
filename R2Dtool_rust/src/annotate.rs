use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use crate::parse_annotation::{Transcript, read_gtf_file, read_gff_file};

// New struct for SpliceSite
#[derive(Debug, Clone)]
pub struct SpliceSite {
    pub transcript_id: String,
    pub tx_coord: u64,
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

// New function to generate splice sites
fn generate_splice_sites(transcripts: &HashMap<String, Transcript>) -> HashMap<String, Vec<SpliceSite>> {
    let mut splice_sites_map: HashMap<String, Vec<SpliceSite>> = HashMap::new();

    for (transcript_id, transcript) in transcripts {
        // println!("Processing transcript: {:?}", transcript_id);
        let transcript_id_no_version = transcript_id.split('.').next().unwrap().to_string();
        let mut splice_sites: Vec<SpliceSite> = Vec::new();

        if !transcript.exons.is_empty() {
            let ref exons = transcript.exons;
            // println!("Exons for transcript {}: {:?}", transcript_id, exons); // Add this line to print exons
            let mut cumsum_width = 0;

            for (i, exon) in exons.iter().enumerate() {
                let exon_length = exon.end - exon.start + 1;
                // println!("Exon {}: start {}, end {}, length {}", i, exon.start, exon.end, exon_length);

                cumsum_width += exon_length;

                if i < exons.len() - 1 {
                    splice_sites.push(SpliceSite {
                        transcript_id: transcript_id_no_version.clone(),
                        tx_coord: cumsum_width as u64,
                    });
                }
                // println!("Splice sites after exon {}: {:?}", i, splice_sites);
            }
        }
        // println!("Splice sites for transcript {}: {:?}", transcript_id_no_version, splice_sites);
        splice_sites_map.insert(transcript_id_no_version, splice_sites);
    }

    splice_sites_map
}




// New function to find splice site distances
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


// Updated run_annotate function
pub fn run_annotate(matches: &clap::ArgMatches) {
    println!("Running the annotate functionality...");

    let gff_file = matches.value_of("gff").unwrap();
    let input_file = matches.value_of("input").unwrap();
    let output_file = matches.value_of("output");
    let format = matches.value_of("format").unwrap_or("gff");

    let annotations = if format == "gtf" {
        read_gtf_file(gff_file)
    } else {
        read_gff_file(gff_file)
    };

    let transcripts = annotations;


    let mut input_reader = BufReader::new(File::open(input_file).expect("Cannot open input file"));
    let has_header = matches.is_present("header");

    let mut output_writer: Box<dyn Write> = match output_file {
        Some(file) => Box::new(File::create(file).expect("Cannot create output file")),
        None => Box::new(std::io::stdout()),
    };

    // preview transcripts hashmap
    // println!("Transcripts: {:?}", transcripts);

    // Generate splice sites
    let splice_sites = generate_splice_sites(&transcripts);

    // preview splice sites and exit
    // println!("{:?}", splice_sites);
    // std::process::exit(0);


    let mut header = String::new();
    if has_header {
        input_reader.read_line(&mut header).unwrap();
        let header_fields: Vec<&str> = header.trim().split('\t').collect();

        let output_header = format!(
            "{}\trel_pos\tabs_cds_start\tabs_cds_end\tupstream_ss_distance\tdownstream_ss_distance",
            header_fields.join("\t")
        );
        writeln!(output_writer, "{}", output_header).unwrap();
    }

    for line in input_reader.lines() {
        let line = line.unwrap();
        let fields: Vec<&str> = line.split('\t').collect();
        let transcript_id_with_version = fields[0];
        let transcript_id = transcript_id_with_version.split('.').next().unwrap(); // Add this line to remove version from transcript ID
        let tx_coord: u64 = fields[1].parse().unwrap();
        // println!("Processing transcript ID: {:?}", transcript_id); // Add this print statement

        if let Some(transcript) = transcripts.get(transcript_id) {
            // println!("Found transcript: {:?}", transcript); // Add this print statement
            if let (Some(utr5_len), Some(cds_len), Some(utr3_len)) = (transcript.utr5_len, transcript.cds_len, transcript.utr3_len) {
                let (rel_pos, abs_cds_start, abs_cds_end) = calculate_meta_coordinates(tx_coord, utr5_len, cds_len, utr3_len);

                // Calculate splice site distances
                if let Some(splice_sites) = splice_sites.get(transcript_id) {
                    // println!("Found splice sites: {:?}", splice_sites); // Add this print statement
                    let (upstream_ss_distance, downstream_ss_distance) = splice_site_distances(tx_coord, splice_sites);

                    let output_line = format!(
                        "{}\t{}\t{}\t{}\t{}\t{}",
                        line,
                        rel_pos,
                        abs_cds_start,
                        abs_cds_end,
                        upstream_ss_distance.unwrap_or(-1),
                        downstream_ss_distance.unwrap_or(-1)
                    );
                    writeln!(output_writer, "{}", output_line).unwrap();
                } else {
                    writeln!(output_writer, "{}\tNA\tNA\tNA\tNA\tNA", line).unwrap();
                }
            } else {
                writeln!(output_writer, "{}\tNA\tNA\tNA\tNA\tNA", line).unwrap();
            }
        } else {
            writeln!(output_writer, "{}\tNA\tNA\tNA\tNA\tNA", line).unwrap();
        }
    }
}
