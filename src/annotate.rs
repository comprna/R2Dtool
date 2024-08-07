use std::fs;
use std::fs::{File,OpenOptions};
use std::io::{BufRead, BufReader, Write};
use std::error::Error;
use std::collections::HashMap;
use crate::parse_gtf::{Transcript, read_annotation_file, Exon};
use rayon::prelude::*;
use std::sync::{Arc, Mutex};


#[derive(Debug, Clone)]
pub struct SpliceSite {
    pub transcript_id: String,
    pub tx_coord: u64,
}

fn calculate_cds_start(utr5_len: u64) -> u64 {
    utr5_len
}

fn calculate_cds_end(utr5_len: u64, cds_len: u64) -> u64 {
    utr5_len + cds_len
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


fn generate_splice_sites(transcripts: &HashMap<String, Transcript>) -> Arc<Mutex<HashMap<String, Vec<SpliceSite>>>> {
    let splice_sites_map = Arc::new(Mutex::new(HashMap::new()));

    let file_path = "splice_sites_map.txt";
    let file = OpenOptions::new().write(true).create(true).truncate(true).open(file_path).expect("Unable to open file");
    let file = Arc::new(Mutex::new(file));

    transcripts.par_iter().for_each(|(transcript_id, transcript)| {
        let transcript_id_no_version = transcript_id.split('.').next().unwrap().to_string();
        let mut splice_sites: Vec<SpliceSite> = Vec::new();
        let mut cumulative_length: u64 = 0;

        if !transcript.exons.is_empty() {
            let mut exons: Vec<&Exon> = transcript.exons.iter()
                .filter(|exon| exon.feature.as_ref().map_or(false, |ft| ft == "exon"))
                .collect();

            if transcript.strand.as_deref() == Some("+") {
                exons.sort_by_key(|exon| exon.start);
            } else if transcript.strand.as_deref() == Some("-") {
                exons.sort_by_key(|exon| std::cmp::Reverse(exon.start));
            }

            for (i, exon) in exons.iter().enumerate() {

                // Update cumulative length to include exon length
                cumulative_length += exon.length;
                if i < exons.len() - 1 {  // Ignore the last exon 
                    splice_sites.push(SpliceSite {
                        transcript_id: transcript_id_no_version.clone(),
                        tx_coord: cumulative_length,  // Position within the transcript
                    });
                }
                // Debug output 
                // println!("Transcript ID: {}, Exon Length: {}, Cumulative Length: {}", transcript_id_no_version, exon.length, cumulative_length);
            }
        }

        // Writing to file (thread-safe)
        let mut file = file.lock().unwrap();
        writeln!(file, "Transcript ID: {}", transcript_id_no_version).expect("Unable to write to file");
        for splice_site in &splice_sites {
            writeln!(file, "Splice Site: {}", splice_site.tx_coord).expect("Unable to write to file");
        }
        writeln!(file).expect("Unable to write to file"); // Write a newline for better readability

        let mut map = splice_sites_map.lock().unwrap();
        map.insert(transcript_id_no_version, splice_sites);
    });

    splice_sites_map
}

fn splice_site_distances(tx_coord: u64, splice_sites: &[SpliceSite]) -> (Option<i64>, Option<i64>) {
    
    // Initialize upstream and downstream distances as None
    let mut upstream_distance: Option<i64> = None;
    let mut downstream_distance: Option<i64> = None;

    // Iterate over all splice sites
    for splice_site in splice_sites {
        
        // Calculate the distance from the current transcript coordinate to the splice site
        let site_distance = splice_site.tx_coord as i64 - tx_coord as i64;

        // Determine if the site is upstream or downstream
        // If site_distance is negative, the splice site is upstream
        if site_distance < 0 {
            let positive_distance = site_distance.abs();
            if upstream_distance.is_none() || positive_distance < upstream_distance.unwrap() {
                upstream_distance = Some(positive_distance);
            }
        }
        // Similar logic for downstream 
        else if site_distance > 0 {
            if downstream_distance.is_none() || site_distance < downstream_distance.unwrap() {
                downstream_distance = Some(site_distance);
            }
        }
    }

    // Return the smallest positive distances for upstream and downstream junctions
    (upstream_distance, downstream_distance)
}


pub fn run_annotate(matches: &clap::ArgMatches, has_header: bool, has_version: bool) -> Result<(), Box<dyn Error>> {
   
    // eprintln!("Running the annotate functionality...");
   
    let gtf_file: String = matches.get_one::<String>("gtf").unwrap().to_string();
    let input_file: String = matches.get_one::<String>("input").unwrap().to_string();
    let output_file: Option<String> = matches.get_one::<String>("output").map(|s: &String| s.to_string());

    // By default, read in the annotations as GTF file
    // TODO: implement GFF3 parsing         
    // let default_format = String::from("gtf");
    // let format = matches.get_one("format").unwrap_or(&default_format);
    let annotations = read_annotation_file(&gtf_file, true, has_version)?;
    
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
    "{}\tgene_id\tgene_name\ttranscript_biotype\ttx_len\tcds_start\tcds_end\ttx_end\ttranscript_metacoordinate\tabs_cds_start\tabs_cds_end\tup_junc_dist\tdown_junc_dist",
    header_fields.join("\t")
);
        writeln!(output_writer, "{}", output_header).unwrap();
    }

    for line in input_reader.lines() {
        let line = line.unwrap();
        let fields: Vec<&str> = line.split('\t').collect();
        let transcript_id_with_version = fields[0];


        let transcript_id = if has_version {
            transcript_id_with_version
        } else {
            transcript_id_with_version.split('.').next().unwrap()
        };
    

        let tx_coord: u64 = fields[1].parse().unwrap();
    
        if let Some(transcript) = transcripts.get(transcript_id) {
            
            // Initialize all fields to "NA"
            let tx_len;
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
    
            tx_len = transcript.transcript_length.map_or("NA".to_string(), |len| len.to_string());
            if let (Some(utr5_len), Some(cds_len), Some(utr3_len)) = (transcript.utr5_len, transcript.cds_len, transcript.utr3_len) {
                cds_start = calculate_cds_start(utr5_len).to_string();
                cds_end = calculate_cds_end(utr5_len, cds_len).to_string();
                tx_end = tx_len.clone();
                let calculated_values = calculate_meta_coordinates(tx_coord, utr5_len, cds_len, utr3_len);
                rel_pos = format!("{:.5}", calculated_values.0);
                abs_cds_start = calculated_values.1.to_string();
                abs_cds_end = calculated_values.2.to_string();
            }

            // Lock the mutex to safely access the shared hashmap
            let map_lock = splice_sites.lock().expect("Failed to lock the mutex");

            // Handle splice sites if available
            if let Some(splice_sites) = map_lock.get(transcript_id) {
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
    let file_path = "splice_sites_map.txt";
    if let Err(e) = fs::remove_file(file_path) {
        eprintln!("Failed to delete file '{}': {:?}", file_path, e);
    }

    Ok(())
}


pub fn preview_annotations(annotations: &HashMap<String, Transcript>) {
    eprintln!("Number of annotations: {}", annotations.len()); 
    for (key, transcript) in annotations {
        eprintln!("Annotations start");
        eprintln!("Transcript ID: {}", key);
        eprintln!("{:#?}", transcript);
        eprintln!("Annotations end");
    }
}

// Unit tests 
#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::NamedTempFile;

    #[test]
    fn test_calculate_cds_end() {
        assert_eq!(calculate_cds_end(100, 200), 300);
        assert_eq!(calculate_cds_end(0, 100), 100);
    }

    #[test]
    fn test_calculate_meta_coordinates() {
        // Test UTR5 region
        let (rel_pos, abs_cds_start, abs_cds_end) = calculate_meta_coordinates(50, 100, 200, 100);
        assert_eq!(rel_pos, 0.5);
        assert_eq!(abs_cds_start, -50);
        assert_eq!(abs_cds_end, -250);

        // Test CDS region
        let (rel_pos, abs_cds_start, abs_cds_end) = calculate_meta_coordinates(150, 100, 200, 100);
        assert_eq!(rel_pos, 1.25);
        assert_eq!(abs_cds_start, 50);
        assert_eq!(abs_cds_end, -150);

        // Test UTR3 region
        let (rel_pos, abs_cds_start, abs_cds_end) = calculate_meta_coordinates(350, 100, 200, 100);
        assert_eq!(rel_pos, 2.5);
        assert_eq!(abs_cds_start, 250);
        assert_eq!(abs_cds_end, 50);
    }

    #[test]
    fn test_splice_site_distances() {
        let splice_sites = vec![
            SpliceSite { transcript_id: "test".to_string(), tx_coord: 50 },
            SpliceSite { transcript_id: "test".to_string(), tx_coord: 100 },
            SpliceSite { transcript_id: "test".to_string(), tx_coord: 150 },
        ];

        // Test coordinate before all splice sites
        let (up, down) = splice_site_distances(25, &splice_sites);
        assert_eq!(up, None);
        assert_eq!(down, Some(25));

        // Test coordinate between splice sites
        let (up, down) = splice_site_distances(75, &splice_sites);
        assert_eq!(up, Some(25));
        assert_eq!(down, Some(25));

        // Test coordinate after all splice sites
        let (up, down) = splice_site_distances(200, &splice_sites);
        assert_eq!(up, Some(50));
        assert_eq!(down, None);
    }

    #[test]
    fn test_run_annotate() {
        use std::path::Path;
    
        // Use actual files from the ./test/ directory
        let input_file = Path::new("./test/m6A_isoform_sites_GRCh38_subset.bed");
        let gtf_file = Path::new("./test/GRCh38.110_subset.gtf");
        let output_file = NamedTempFile::new().unwrap();
    
        // Create mock ArgMatches
        let matches = clap::Command::new("test")
            .arg(clap::Arg::new("gtf").short('g').long("gtf").required(true))
            .arg(clap::Arg::new("input").short('i').long("input").required(true))
            .arg(clap::Arg::new("output").short('o').long("output").required(true))
            .get_matches_from(vec![
                "test",
                "-g", gtf_file.to_str().unwrap(),
                "-i", input_file.to_str().unwrap(),
                "-o", output_file.path().to_str().unwrap(),
            ]);
    
        // Run annotate
        let result = run_annotate(&matches, true, false);
        assert!(result.is_ok());
    
        // Read and check output
        let output = std::fs::read_to_string(output_file.path()).unwrap();
        let output_lines: Vec<&str> = output.lines().take(2).collect();
        assert_eq!(output_lines.len(), 2); // Header + 1 data line
    
        // Check header
        let expected_header = "transcript\tstart\tend\tbase\tcoverage\tstrand\tN_valid_cov\tfraction_modified\tgene_id\tgene_name\ttranscript_biotype\ttx_len\tcds_start\tcds_end\ttx_end\ttranscript_metacoordinate\tabs_cds_start\tabs_cds_end\tup_junc_dist\tdown_junc_dist";
        assert_eq!(output_lines[0], expected_header);
    
        // Check data line
        let data_fields: Vec<&str> = output_lines[1].split('\t').collect();
        assert_eq!(data_fields.len(), 20); // Ensure we have the correct number of fields
    
        // Check specific fields
        assert_eq!(data_fields[0], "ENST00000381989.4");
        assert_eq!(data_fields[1], "2682");
        assert_eq!(data_fields[2], "2683");
        assert_eq!(data_fields[3], "a");
        assert_eq!(data_fields[4], "10");
        assert_eq!(data_fields[5], "+");
        assert_eq!(data_fields[6], "10");
        assert_eq!(data_fields[7], "0.00");
        assert_eq!(data_fields[8], "ENSG00000102699");
        assert_eq!(data_fields[9], "PARP4");
        assert_eq!(data_fields[10], "protein_coding");
    
        // Check numeric fields 
        assert!(data_fields[11].parse::<u64>().is_ok()); // tx_len
        assert!(data_fields[12].parse::<u64>().is_ok()); // cds_start
        assert!(data_fields[13].parse::<u64>().is_ok()); // cds_end
        assert!(data_fields[14].parse::<u64>().is_ok()); // tx_end
        
        let metacoordinate: f64 = data_fields[15].parse().unwrap();
        assert!((metacoordinate - 1.50425).abs() < 0.00001);
    
        assert!(data_fields[16].parse::<i64>().is_ok()); // abs_cds_start
        assert!(data_fields[17].parse::<i64>().is_ok()); // abs_cds_end
        assert!(data_fields[18].parse::<u64>().is_ok()); // up_junc_dist
        assert!(data_fields[19].parse::<u64>().is_ok()); // down_junc_dist
    
        println!("Actual output: {}", output_lines[1]);
    }
}