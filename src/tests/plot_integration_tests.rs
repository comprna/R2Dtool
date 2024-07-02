use std::path::{PathBuf, Path};
use tempfile::NamedTempFile;
use crate::{annotate, generate_r_plots};
use std::env;
use std::fs;

fn get_project_root() -> PathBuf {
    Path::new(&env!("CARGO_MANIFEST_DIR")).to_path_buf()
}

fn create_output_directory() -> std::io::Result<PathBuf> {
    let output_dir = get_project_root().join("test_outputs");
    fs::create_dir_all(&output_dir)?;
    Ok(output_dir)
}

fn generate_output_filename(plot_type: &str, cutoff: &str, direction: &str, labels: bool, ci_method: &str, codon_type: Option<&str>) -> String {
    let mut filename = format!("{}_cutoff{}_{}_{}", plot_type, cutoff, direction, ci_method);
    if labels {
        filename.push_str("_with_labels");
    }
    if let Some(codon) = codon_type {
        filename.push_str(&format!("_{}_codon", codon));
    }
    filename.push_str(".png");
    filename
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

fn setup_test() -> Result<(PathBuf, PathBuf, NamedTempFile), Box<dyn std::error::Error>> {
    let project_root = get_project_root();
    let script_dir = project_root.join("scripts");
    let test_dir = project_root.join("test");

    let annotate_output = NamedTempFile::new()?;
    run_annotate(
        &test_dir.join("m6A_isoform_sites_GRCh38_subset.bed"),
        &test_dir.join("GRCh38.110_subset.gtf"),
        annotate_output.path()
    )?;

    Ok((script_dir, test_dir, annotate_output))
}

fn run_plot_test(script_name: &str, args: Vec<String>, script_dir: &Path) -> Result<(), Box<dyn std::error::Error>> {
    generate_r_plots(script_name, &args, Some(script_dir))?;
    let output_path = Path::new(&args[1]);
    assert!(output_path.exists());
    Ok(())
}

// metatranscript tests 
#[test]
fn test_plot_meta_transcript_cutoff10_upper_loess_with_labels() -> Result<(), Box<dyn std::error::Error>> {
    let (script_dir, _, annotate_output) = setup_test()?;
    let output_dir = create_output_directory()?;
    let filename = generate_output_filename("metatranscript", "10", "upper", true, "loess", None);
    let output_path = output_dir.join(filename);
    let args = vec![
        annotate_output.path().to_str().unwrap().to_string(),
        output_path.to_str().unwrap().to_string(),
        "fraction_modified".to_string(),
        "10".to_string(),
        "upper".to_string(),
        "-l".to_string(),
    ];
    run_plot_test("R2_plotMetaTranscript.R", args, &script_dir)
}

#[test]
fn test_plot_meta_transcript_cutoff0_lower_binom_without_labels() -> Result<(), Box<dyn std::error::Error>> {
    let (script_dir, _, annotate_output) = setup_test()?;
    let output_dir = create_output_directory()?;
    let filename = generate_output_filename("metatranscript", "0", "lower", false, "binom", None);
    let output_path = output_dir.join(filename);
    let args = vec![
        annotate_output.path().to_str().unwrap().to_string(),
        output_path.to_str().unwrap().to_string(),
        "fraction_modified".to_string(),
        "0".to_string(),
        "lower".to_string(),
        "-c".to_string(),
        "binom".to_string(),
    ];
    run_plot_test("R2_plotMetaTranscript.R", args, &script_dir)
}

#[test]
fn test_plot_meta_transcript_cutoff100_upper_loess_without_labels() -> Result<(), Box<dyn std::error::Error>> {
    let (script_dir, _, annotate_output) = setup_test()?;
    let output_dir = create_output_directory()?;
    let filename = generate_output_filename("metatranscript", "100", "upper", false, "loess", None);
    let output_path = output_dir.join(filename);
    let args = vec![
        annotate_output.path().to_str().unwrap().to_string(),
        output_path.to_str().unwrap().to_string(),
        "fraction_modified".to_string(),
        "100".to_string(),
        "upper".to_string(),
    ];
    run_plot_test("R2_plotMetaTranscript.R", args, &script_dir)
}

// metajunction tests 
#[test]
fn test_plot_meta_junction_cutoff10_upper_loess() -> Result<(), Box<dyn std::error::Error>> {
    let (script_dir, _, annotate_output) = setup_test()?;
    let output_dir = create_output_directory()?;
    let filename = generate_output_filename("metajunction", "10", "upper", false, "loess", None);
    let output_path = output_dir.join(filename);
    let args = vec![
        annotate_output.path().to_str().unwrap().to_string(),
        output_path.to_str().unwrap().to_string(),
        "fraction_modified".to_string(),
        "10".to_string(),
        "upper".to_string(),
    ];
    run_plot_test("R2_plotMetaJunction.R", args, &script_dir)
}

#[test]
fn test_plot_meta_junction_cutoff0_lower_binom() -> Result<(), Box<dyn std::error::Error>> {
    let (script_dir, _, annotate_output) = setup_test()?;
    let output_dir = create_output_directory()?;
    let filename = generate_output_filename("metajunction", "0", "lower", false, "binom", None);
    let output_path = output_dir.join(filename);
    let args = vec![
        annotate_output.path().to_str().unwrap().to_string(),
        output_path.to_str().unwrap().to_string(),
        "fraction_modified".to_string(),
        "0".to_string(),
        "lower".to_string(),
        "-c".to_string(),
        "binom".to_string(),
    ];
    run_plot_test("R2_plotMetaJunction.R", args, &script_dir)
}

// metacodon tests
#[test]
fn test_plot_meta_codon_start_cutoff10_upper_loess() -> Result<(), Box<dyn std::error::Error>> {
    let (script_dir, _, annotate_output) = setup_test()?;
    let output_dir = create_output_directory()?;
    let filename = generate_output_filename("metacodon", "10", "upper", false, "loess", Some("start"));
    let output_path = output_dir.join(filename);
    let args = vec![
        annotate_output.path().to_str().unwrap().to_string(),
        output_path.to_str().unwrap().to_string(),
        "fraction_modified".to_string(),
        "10".to_string(),
        "upper".to_string(),
        "-s".to_string(),
    ];
    run_plot_test("R2_plotMetaCodon.R", args, &script_dir)
}

#[test]
fn test_plot_meta_codon_stop_cutoff0_lower_binom() -> Result<(), Box<dyn std::error::Error>> {
    let (script_dir, _, annotate_output) = setup_test()?;
    let output_dir = create_output_directory()?;
    let filename = generate_output_filename("metacodon", "0", "lower", false, "binom", Some("stop"));
    let output_path = output_dir.join(filename);
    let args = vec![
        annotate_output.path().to_str().unwrap().to_string(),
        output_path.to_str().unwrap().to_string(),
        "fraction_modified".to_string(),
        "0".to_string(),
        "lower".to_string(),
        "-e".to_string(),
        "-c".to_string(),
        "binom".to_string(),
    ];
    run_plot_test("R2_plotMetaCodon.R", args, &script_dir)
}

#[test]
fn test_plot_meta_codon_start_cutoff100_upper_loess() -> Result<(), Box<dyn std::error::Error>> {
    let (script_dir, _, annotate_output) = setup_test()?;
    let output_dir = create_output_directory()?;
    let filename = generate_output_filename("metacodon", "100", "upper", false, "loess", Some("start"));
    let output_path = output_dir.join(filename);
    let args = vec![
        annotate_output.path().to_str().unwrap().to_string(),
        output_path.to_str().unwrap().to_string(),
        "fraction_modified".to_string(),
        "100".to_string(),
        "upper".to_string(),
        "-s".to_string(),
    ];
    run_plot_test("R2_plotMetaCodon.R", args, &script_dir)
}

#[test]
fn test_plot_meta_codon_start_cutoff100_lower_loess() -> Result<(), Box<dyn std::error::Error>> {
    let (script_dir, _, annotate_output) = setup_test()?;
    let output_dir = create_output_directory()?;
    let filename = generate_output_filename("metacodon", "100", "lower", false, "loess", Some("start"));
    let output_path = output_dir.join(filename);
    let args = vec![
        annotate_output.path().to_str().unwrap().to_string(),
        output_path.to_str().unwrap().to_string(),
        "fraction_modified".to_string(),
        "100".to_string(),
        "lower".to_string(),
        "-s".to_string(),
    ];
    run_plot_test("R2_plotMetaCodon.R", args, &script_dir)

// helper function for saving plot tables 
fn run_plot_test_with_table(script_name: &str, args: Vec<String>, script_dir: &Path) -> Result<(), Box<dyn std::error::Error>> {
    let mut table_args = args.clone();
    let table_path = Path::new(&args[1]).with_extension("tsv");
    table_args.extend_from_slice(&["-o".to_string(), table_path.to_str().unwrap().to_string()]);
    generate_r_plots(script_name, &table_args, Some(script_dir))?;
    assert!(Path::new(&args[1]).exists());
    assert!(table_path.exists());
    Ok(())
}

#[test]
fn test_plot_meta_transcript_with_table() -> Result<(), Box<dyn std::error::Error>> {
    let (script_dir, _, annotate_output) = setup_test()?;
    let output_dir = create_output_directory()?;
    let filename = generate_output_filename("metatranscript", "10", "upper", false, "loess", None);
    let output_path = output_dir.join(filename);
    let args = vec![
        annotate_output.path().to_str().unwrap().to_string(),
        output_path.to_str().unwrap().to_string(),
        "fraction_modified".to_string(),
        "10".to_string(),
        "upper".to_string(),
    ];
    run_plot_test_with_table("R2_plotMetaTranscript.R", args, &script_dir)
}

#[test]
fn test_plot_meta_junction_with_table() -> Result<(), Box<dyn std::error::Error>> {
    let (script_dir, _, annotate_output) = setup_test()?;
    let output_dir = create_output_directory()?;
    let filename = generate_output_filename("metajunction", "10", "upper", false, "loess", None);
    let output_path = output_dir.join(filename);
    let args = vec![
        annotate_output.path().to_str().unwrap().to_string(),
        output_path.to_str().unwrap().to_string(),
        "fraction_modified".to_string(),
        "10".to_string(),
        "upper".to_string(),
    ];
    run_plot_test_with_table("R2_plotMetaJunction.R", args, &script_dir)
}

#[test]
fn test_plot_meta_codon_with_table() -> Result<(), Box<dyn std::error::Error>> {
    let (script_dir, _, annotate_output) = setup_test()?;
    let output_dir = create_output_directory()?;
    let filename = generate_output_filename("metacodon", "10", "upper", false, "loess", Some("start"));
    let output_path = output_dir.join(filename);
    let args = vec![
        annotate_output.path().to_str().unwrap().to_string(),
        output_path.to_str().unwrap().to_string(),
        "fraction_modified".to_string(),
        "10".to_string(),
        "upper".to_string(),
        "-s".to_string(),
    ];
    run_plot_test_with_table("R2_plotMetaCodon.R", args, &script_dir)
}