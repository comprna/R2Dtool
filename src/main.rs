#![allow(non_snake_case)]
use clap::{Arg,Command};
use std::process::Command as ProcessCommand;
use std::env;

// modules
pub mod parse_annotation;
pub mod annotate;
pub mod liftover;
pub mod parse_gtf;
use std::path::Path;

#[cfg(test)]
mod tests;

const RELATIVE_SCRIPT_PATH: &str = "../../scripts/";

pub fn generate_r_plots(script_name: &str, args: &[String], script_dir: Option<&Path>) -> Result<(), String> {
    let script_path = if let Some(dir) = script_dir {
        dir.join(script_name)
    } else {
        let current_exe = env::current_exe().map_err(|e| format!("Failed to get current executable path: {}", e))?;
        let script_dir = current_exe.parent()
            .ok_or("Failed to get parent directory of executable")?
            .join(RELATIVE_SCRIPT_PATH);
        script_dir.join(script_name)
    };

    println!("Executing R script: {:?}", script_path);

    let output = ProcessCommand::new("Rscript")
        .arg(&script_path)
        .args(args)
        .output()
        .map_err(|e| format!("Failed to execute R script: {}", e))?;

    // Print stdout
    println!("R script stdout:");
    println!("{}", String::from_utf8_lossy(&output.stdout));

    // Print stderr
    eprintln!("R script stderr:");
    eprintln!("{}", String::from_utf8_lossy(&output.stderr));

    if output.status.success() {
        println!("R script ran successfully.");
        Ok(())
    } else {
        Err(format!("R script exited with non-zero status: {:?}", output.status))
    }
}

fn main() {
    let matches = Command::new("R2Dtool")
        .version("1.0.0")
        .author("AJ Sethi, compRNA <aditya.sethi@anu.edu.au; CompRNA@ANU365.onmicrosoft.com>")
        .about("R2Dtool")
        .arg_required_else_help(true) 
        .subcommand(
            Command::new("liftover")
                .about("Converts transcriptomic to genomic coordinates")
                .arg(
                    Arg::new("gtf")
                    .short('g')
                    .long("gtf")
                    .value_name("GTF_FILE")
                    .help("Path to GTF gene structure annotation")
                    .required(true)
                )
                .arg(
                    Arg::new("input")
                    .short('i')
                    .long("input")
                    .value_name("INPUT_FILE")
                    .help("Path to input file with transcriptomic coordinates")
                    .required(true)
                )
                .arg(
                    Arg::new("header")
                    .short('H')
                    .long("header")
                    .help("Indicates that the input file has a header in line 1")
                    .action(clap::ArgAction::SetTrue)
                )
                .arg(
                    Arg::new("transcript-version")
                    .short('t')
                    .long("transcript-version")
                    .help("Retain transcript version information (. delimited) in col 1")
                    .action(clap::ArgAction::SetTrue)
                )
                .arg(
                    Arg::new("output")
                    .short('o')
                    .long("output")
                    .value_name("OUTPUT_FILE")
                    .help("Path to output file")
                )
                // .arg(
                //     Arg::new("format")
                //     .short('f')
                //     .long("format")
                //     .value_name("FORMAT")
                //     .help("Specify the gene strutcture annotation format: gtf or gff (default: gff)")
                // )
        )
        .subcommand(
            Command::new("annotate")
                .about("Annotates transcriptomic sites with genomic cooridnates")
                .arg(
                    Arg::new("gtf")
                    .short('g')
                    .long("gtf")
                    .value_name("GTF_FILE")
                    .help("Path to GTF gene structure annotation")
                    .required(true)
                )
                .arg(
                    Arg::new("input")
                    .short('i')
                    .long("input")
                    .value_name("INPUT_FILE")
                    .help("Path to input file with transcriptomic coordinates")
                    .required(true)
                )
                .arg(
                    Arg::new("header")
                    .short('H')
                    .long("header")
                    .help("Indicates that the input file has a header in line 1")
                    .action(clap::ArgAction::SetTrue)
                )
                .arg(
                    Arg::new("transcript-version")
                    .short('t')
                    .long("transcript-version")
                    .help("Retain transcript version information (. delimited) in col 1")
                    .action(clap::ArgAction::SetTrue)
                )
                .arg(
                    Arg::new("output")
                    .short('o')
                    .long("output")
                    .value_name("OUTPUT_FILE")
                    .help("Path to output file")
                )
                // .arg(
                //     Arg::new("format")
                //     .short('f')
                //     .long("format")
                //     .value_name("FORMAT")
                //     .help("Specify the gene strutcture annotation format: gtf or gff (default: gff)")
                // )
        )
        .subcommand(
            Command::new("plotMetaTranscript")
                .about("Plot meta transcript distribution")
                .arg(Arg::new("annotated_sites")
                    .short('i')
                    .long("input")
                    .value_name("FILE")
                    .help("Path to annotated file generated by R2Dtool annotate")
                    .required(true))
                .arg(Arg::new("output_path")
                    .short('o')
                    .long("output")
                    .value_name("FILE")
                    .help("Path to output plot (including extension, e.g. .svg or .png)")
                    .required(true))
                .arg(Arg::new("filter_field")
                    .short('f')
                    .long("filter-field")
                    .value_name("FIELD")
                    .help("Column used to select significant sites")
                    .required(true))
                .arg(Arg::new("cutoff")
                    .short('u')
                    .long("cutoff")
                    .value_name("VALUE")
                    .help("Integer value to use for determining significant sites")
                    .required(true))
                .arg(Arg::new("cutoff_type")
                    .short('t')
                    .long("cutoff-type")
                    .value_name("TYPE")
                    .help("Select sites that have a significance lower or higher than the cutoff value (lower/upper)")
                    .required(true))
                .arg(Arg::new("confidence_method")
                    .short('c')
                    .long("confidence-method")
                    .value_name("METHOD")
                    .help("Strategy for displaying confidence intervals: loess (default) or binomial")
                    .default_value("loess"))
                .arg(Arg::new("save_table")
                    .short('s')
                    .long("save-table")
                    .value_name("FILE")
                    .help("Save the aggregated metatranscript data as a tab-separated file"))
                .arg(Arg::new("show_labels")
                    .short('l')
                    .long("show-labels")
                    .help("Display transcript region labels (5' UTR, CDS, 3'UTR) on the plot")
                    .action(clap::ArgAction::SetTrue))
        )
        .subcommand(
            Command::new("plotMetaJunction")
                .about("Plot meta junction distribution")
                .arg(Arg::new("annotated_sites")
                    .short('i')
                    .long("input")
                    .value_name("FILE")
                    .help("Path to annotated file generated by R2Dtool annotate")
                    .required(true))
                .arg(Arg::new("output_path")
                    .short('o')
                    .long("output")
                    .value_name("FILE")
                    .help("Path to output plot (including extension, e.g. .svg or .png)")
                    .required(true))
                .arg(Arg::new("filter_field")
                    .short('f')
                    .long("filter-field")
                    .value_name("FIELD")
                    .help("Column used to select significant sites")
                    .required(true))
                .arg(Arg::new("cutoff")
                    .short('u')
                    .long("cutoff")
                    .value_name("VALUE")
                    .help("Integer value to use for determining significant sites")
                    .required(true))
                .arg(Arg::new("cutoff_type")
                    .short('t')
                    .long("cutoff-type")
                    .value_name("TYPE")
                    .help("Select sites that have a significance lower or higher than the cutoff value (lower/upper)")
                    .required(true))
                .arg(Arg::new("confidence_method")
                    .short('c')
                    .long("confidence-method")
                    .value_name("METHOD")
                    .help("Strategy for displaying confidence intervals: loess (default) or binomial")
                    .default_value("loess"))
                .arg(Arg::new("save_table")
                    .short('s')
                    .long("save-table")
                    .value_name("FILE")
                    .help("Save the aggregated metajunction data as a tab-separated file"))
        )
        .subcommand(
            Command::new("plotMetaCodon")
                .about("Plot meta codon distribution")
                .arg(
                    Arg::new("start_codon")
                        .short('s')
                        .long("start")
                        .help("Plot distribution around start codon")
                        .action(clap::ArgAction::SetTrue)
                )
                .arg(
                    Arg::new("stop_codon")
                        .short('e')
                        .long("stop")
                        .help("Plot distribution around stop codon")
                        .action(clap::ArgAction::SetTrue)
                )
                .arg(
                    Arg::new("input_file")
                        .short('i')
                        .long("input")
                        .value_name("FILE")
                        .help("Path to the annotated transcriptomic sites file")
                        .required(true)
                )
                .arg(
                    Arg::new("output_file")
                        .short('o')
                        .long("output")
                        .value_name("FILE")
                        .help("Path where the plot will be saved (include file extension, e.g., .png or .svg)")
                        .required(true)
                )
                .arg(
                    Arg::new("field_name")
                        .short('f')
                        .long("filter-field")
                        .value_name("FIELD")
                        .help("The name of the column in input_file used to filter significant sites")
                        .required(true)
                )
                .arg(
                    Arg::new("cutoff_value")
                        .short('u')
                        .long("cutoff")
                        .value_name("VALUE")
                        .help("Numeric value defining the threshold for significance")
                        .required(true)
                )
                .arg(
                    Arg::new("cutoff_type")
                        .short('t')
                        .long("cutoff-type")
                        .value_name("TYPE")
                        .help("Specifies the comparison direction, either 'lower' or 'upper', to determine significance")
                        .required(true)
                        .value_parser(["lower", "upper"])
                )
                .arg(
                    Arg::new("confidence_method")
                        .short('c')
                        .long("confidence-method")
                        .value_name("METHOD")
                        .help("Strategy for displaying confidence intervals: loess (default) or binomial")
                        .default_value("loess")
                        .value_parser(["loess", "binom"])
                )
                .arg(
                    Arg::new("save_table")
                        .short('a')
                        .long("save-table")
                        .value_name("FILE")
                        .help("Save the aggregated metacodon data as a tab-separated file")
                )
                .group(
                    clap::ArgGroup::new("codon_type")
                        .required(true)
                        .args(&["start_codon", "stop_codon"]),
                )
        )
        .get_matches();

    // Liftover 
    if let Some(liftover_matches) = matches.subcommand_matches("liftover") {
        let has_header = liftover_matches.get_flag("header");
        let has_version = liftover_matches.get_flag("transcript-version");
        
        eprintln!("Running liftover...");

        if let Err(e) = liftover::run_liftover(liftover_matches, has_header, has_version) {
            eprintln!("Error running liftover: {}", e);
        }
    }

    // Annotate
    if let Some(annotate_matches) = matches.subcommand_matches("annotate") {
        let has_header = annotate_matches.get_flag("header");
        let has_version = annotate_matches.get_flag("transcript-version");

        eprintln!("Running annotate...");
        
        if let Err(e) = annotate::run_annotate(annotate_matches, has_header, has_version) {
            eprintln!("Error running annotate: {}", e);
        }
    }

    // PlotMetaTranscript
    if let Some(matches) = matches.subcommand_matches("plotMetaTranscript") {
        let mut args: Vec<String> = vec![
            matches.get_one::<String>("annotated_sites").unwrap().clone(),
            matches.get_one::<String>("output_path").unwrap().clone(),
            matches.get_one::<String>("filter_field").unwrap().clone(),
            matches.get_one::<String>("cutoff").unwrap().clone(),
            matches.get_one::<String>("cutoff_type").unwrap().clone(),
        ];

        if let Some(method) = matches.get_one::<String>("confidence_method") {
            args.extend_from_slice(&["-c".to_string(), method.clone()]);
        }
        if matches.contains_id("save_table") {
            args.extend_from_slice(&["-o".to_string(), matches.get_one::<String>("save_table").unwrap().clone()]);
        }
        if matches.get_flag("show_labels") {
            args.push("-l".to_string());
        }

        println!("Arguments being passed to R2_plotMetaTranscript.R:");
        println!("Rscript R2_plotMetaTranscript.R {}", args.join(" "));

        match generate_r_plots("R2_plotMetaTranscript.R", &args) {
            Ok(_) => println!("PlotMetaTranscript generated successfully."),
            Err(e) => {
                eprintln!("Error: {}", e);
                eprintln!("PlotMetaTranscript generation failed. Please check the error messages above.");
                std::process::exit(1);
            }
        }
    }

    // plotMetaJunction 
    if let Some(matches) = matches.subcommand_matches("plotMetaJunction") {
        let mut args: Vec<String> = vec![
            matches.get_one::<String>("annotated_sites").unwrap().clone(),
            matches.get_one::<String>("output_path").unwrap().clone(),
            matches.get_one::<String>("filter_field").unwrap().clone(),
            matches.get_one::<String>("cutoff").unwrap().clone(),
            matches.get_one::<String>("cutoff_type").unwrap().clone(),
        ];

        if let Some(method) = matches.get_one::<String>("confidence_method") {
            args.extend_from_slice(&["-c".to_string(), method.clone()]);
        }
        if let Some(table) = matches.get_one::<String>("save_table") {
            args.extend_from_slice(&["-o".to_string(), table.clone()]);
        }

        println!("Arguments being passed to R2_plotMetaJunction.R:");
        println!("Rscript R2_plotMetaJunction.R {}", args.join(" "));

        match generate_r_plots("R2_plotMetaJunction.R", &args) {
            Ok(_) => println!("PlotMetaJunction generated successfully."),
            Err(e) => {
                eprintln!("Error: {}", e);
                eprintln!("PlotMetaJunction generation failed. Please check the error messages above.");
                std::process::exit(1);
            }
        }
    }
    
// plotMetaCodon 
if let Some(matches) = matches.subcommand_matches("plotMetaCodon") {
        let codon_flag = if matches.get_flag("start_codon") { "-s" } else { "-e" };

        let mut args: Vec<String> = vec![
            matches.get_one::<String>("input_file").unwrap().clone(),
            matches.get_one::<String>("output_file").unwrap().clone(),
            matches.get_one::<String>("field_name").unwrap().clone(),
            matches.get_one::<String>("cutoff_value").unwrap().clone(),
            matches.get_one::<String>("cutoff_type").unwrap().clone(),
            codon_flag.to_string(),
        ];

        if let Some(method) = matches.get_one::<String>("confidence_method") {
            args.extend_from_slice(&["-c".to_string(), method.clone()]);
        }

        if let Some(table) = matches.get_one::<String>("save_table") {
            args.extend_from_slice(&["-o".to_string(), table.clone()]);
        }

        println!("Arguments being passed to R2_plotMetaCodon.R:");
        println!("Rscript R2_plotMetaCodon.R {}", args.join(" "));
    
        match generate_r_plots("R2_plotMetaCodon.R", &args) {
            Ok(_) => println!("PlotMetaCodon generated successfully."),
            Err(e) => {
                eprintln!("Error: {}", e);
                eprintln!("PlotMetaCodon generation failed. Please check the error messages above.");
                std::process::exit(1);
            }
        }
    }
}
