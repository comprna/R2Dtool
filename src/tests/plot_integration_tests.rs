use std::path::{PathBuf, Path};
use tempfile::NamedTempFile;
use crate::{annotate, generate_r_plots};
use std::env;

fn get_project_root() -> PathBuf {
    Path::new(&env!("CARGO_MANIFEST_DIR")).to_path_buf()
}

fn run_annotate(input: &Path, gtf: &Path, output: &Path) -> Result<(), Box<dyn std::error::Error>> {
    let matches = clap::Command::new("test")
        .arg(clap::Arg::new("gtf").short('g').long("gtf").required(true))
        .arg(clap::Arg::new("input").short('i').long("input").required(true))
        .arg(clap::Arg::new("output").short('o').long("output").required(true))
        .get_matches_from(vec![
            "test",
            "-g", gtf.to_str().unwrap(),
            "-i", input.to_str().unwrap(),
            "-o", output.to_str().unwrap(),
        ]);

    annotate::run_annotate(&matches, true, false)
}

#[test]
fn test_plotting_functions() -> Result<(), Box<dyn std::error::Error>> {
    let project_root = get_project_root();
    let script_dir = project_root.join("scripts");
    let test_dir = project_root.join("test");

    // Run annotate and save to a temp file
    let annotate_output = NamedTempFile::new()?;
    run_annotate(
        &test_dir.join("m6A_isoform_sites_GRCh38_subset.bed"),
        &test_dir.join("GRCh38.110_subset.gtf"),
        annotate_output.path()
    )?;

    // Test plotMetaTranscript with different configurations
    let meta_transcript_plots = [
        ("metatranscript_m6A_upper.png", "upper", true),
        ("metatranscript_m6A_upper_no_labels.png", "upper", false),
        ("metatranscript_m6A_lower.png", "lower", true),
    ];

    for (_, cutoff_type, show_labels) in &meta_transcript_plots {
        let output_path = NamedTempFile::new()?;
        let mut args = vec![
            annotate_output.path().to_str().unwrap().to_string(),
            output_path.path().to_str().unwrap().to_string(),
            "fraction_modified".to_string(),
            "10".to_string(),
            cutoff_type.to_string(),
        ];
        if *show_labels {
            args.push("-l".to_string());
        }
        generate_r_plots("R2_plotMetaTranscript.R", &args, Some(&script_dir))?;
        assert!(output_path.path().exists());
    }

    // Test plotMetaJunction with different configurations
    let meta_junction_plots = [
        ("metajunction_m6A_upper.png", "upper"),
        ("metajunction_m6A_lower.png", "lower"),
    ];

    for (_, cutoff_type) in &meta_junction_plots {
        let output_path = NamedTempFile::new()?;
        let args = vec![
            annotate_output.path().to_str().unwrap().to_string(),
            output_path.path().to_str().unwrap().to_string(),
            "fraction_modified".to_string(),
            "10".to_string(),
            cutoff_type.to_string(),
        ];
        generate_r_plots("R2_plotMetaJunction.R", &args, Some(&script_dir))?;
        assert!(output_path.path().exists());
    }

    // Test plotMetaCodon for start and stop codons
    let meta_codon_plots = [
        ("metacodon_start_upper.png", "upper", "-s"),
        ("metacodon_start_lower.png", "lower", "-s"),
        ("metacodon_stop_upper.png", "upper", "-e"),
        ("metacodon_stop_lower.png", "lower", "-e"),
    ];

    for (_, cutoff_type, codon_type) in &meta_codon_plots {
        let output_path = NamedTempFile::new()?;
        let args = vec![
            annotate_output.path().to_str().unwrap().to_string(),
            output_path.path().to_str().unwrap().to_string(),
            "fraction_modified".to_string(),
            "10".to_string(),
            cutoff_type.to_string(),
            codon_type.to_string(),
        ];
        generate_r_plots("R2_plotMetaCodon.R", &args, Some(&script_dir))?;
        assert!(output_path.path().exists());
    }

    Ok(())
}