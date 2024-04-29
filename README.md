# R2Dtool

Genomics utilities for handling, integrating, and viualising isoform-mapped RNA feature data.    


  - **Integrate** transcriptome-mapped data: R2Dtool performs liftover of transcriptome-mapped RNA features to their corresponding genomic coordinates

  - **Annotate** transcriptome-mapped sites: R2Dtool annotates transcript-specific metatranscript coordinates and absolute and relative distances to annotated transcript landmarks, in an isoform-specific manner.

  - **Visualise** isoform-aware RNA feature distributions: R2Dtool introduces isoform-aware **metatranscript** plots and **metajunction** plots to study the positonal distribution of RNA features around annotated RNA landmarks.


> [!NOTE]
> [Link to R2Dtool preprint on biorXiv](https://doi.org/10.1101/2022.09.23.509222)    
> Sethi, A. J., Mateos, P. A., Hayashi, R., Shirokikh, N. & Eyras, E. R2Dtool: Positional Interpretation of RNA-Centric Information in the Context of Transcriptomic and Genomic Features. biorXiv 2022.09.23.509222 (2022)
> doi:10.1101/2022.09.23.509222.

# Contents 

   - [Usage](#usage)
       - [Handling isoform-mapped RNA sites](#handling-isoform-mapped-rna-sites)
       - [Visualising isoform-aware RNA feature metadistributions](#visualising-isoform-aware-rna-feature-metadistributions)
   - [Installation](#installation)
   - [Input and output data structure](#input-and-output-data-structure)
   - [Utilities: Convert CHEUI model II output to a BED-like input](#utilities)
   - [General notes](#general-notes)

     
# Usage

## Handling isoform-mapped RNA sites 

**Liftover** transcriptome-mapped features to genomic coordinates:

```
Usage: r2d liftover -i <input> -g <gtf>

Arguments:
    -i, --input <input>: Path to tab-separated transcriptome sites in BED format. 
    -g, --gtf <annotation>: Path to gene structure annotation in GTF format.

Options:
    -H, --header: Indicates the input file has a header, which will be preserved in the output [Default: False]
    -o, --output <OUTPUT>: Path to output file [Default: STDOUT]
    -t, --transscript-version: Indicates that '.'-delimited transcript version information is present in col1 and should be considered during liftover [default: False].
```
- Liftover prepends 6 columns to the input file, containing the genome coordinates of the transcript features in BED format
- All data in the original input are preserved in the output and shifted by 6 columns 

**Annotate** transcriptome-mapped sites with isoform-specific distances to transcript landmarks:

```
Usage: r2d annotate -i <input> -g <gtf>

Arguments:
    -i, --input <input>: Path to tab-separated transcriptome sites in BED format. 
    -g, --gtf <annotation>: Path to gene structure annotation in GTF format.

Options:
    -H, --header: Indicates the input file has a header, which will be preserved in the output [Default: False]
    -o, --output <OUTPUT>: Path to output file [Default: STDOUT]
    -t, --transscript-version: Indicates that '.'-delimited transcript version is present in col1 and should be considered during annotation [default: False].

```

Annotation adds the following information to the epitranscriptomic sites as additional coluumns, relying on the gene structure GTF to generate these data.

```
transcript_biotype | gene_name | gene_id | tx_len | cds_len | utr5_len | utr3_len | cds_start | cds_end | tx_end  | transcript_metacoordinate | abs_cds_start | abs_cds_end | up_junc_dist | down_junc_dist
```

- ```tx_len```, ```cds-len```, ```utr5-len``` and ```utr3-len```represents the lengths of the entire transcript, the coding sequence, and the 5' and 3' untranslated regions (UTRs)
- ```cds_start``` and ```cds_end``` represent the positions of the coding sequence start and end compared to the transcript.
- ```transcript_metacoordinate``` represents the scaled metatrascript position of the given RNA feature, between 0 and 3, where 0 represents transcript start-site, 1 represents CDS start, 2 represent CDS end, and 3 represents the 3' transcript end. 
- ```abs_cds_start``` and ```abs_cds_end``` represent the absolute distance (in nt) of a given feature from the cds start and end
- ```up_junc_dist``` and ```down_junc_dist``` repreesnt the absolute distance (in nt) of a given site from the nearest upstream and downstream splice-junction contained in a given transcript

> [!NOTE]
> - ```liftover``` and ```annotate``` requires columns 1-6 to contain feature coordinates against a transcriptome reference
> - ```annotate``` can be perfomed before, but __not__ after, ```liftover```
> - ```annotate``` requires protein-coding gene models in order to calculate feature metacoordinates and distances to CDS starts and ends 
> - The GTF file used must exactly match the transcriptome to which features are mapped 

## Visualising isoform-aware RNA feature metadistributions

**Plot** the **metatranscript** distribution of RNA features:

```

Rscript R2_plotMetaTranscript.R [annotated transcriptomic sites] [save path for plot, including .png/.svg extension] [name of filter field] [cutoff for significant sites] [lower/upper cutoff for significant sites]

Mandatory positional arguments:
- path to annotated file generated by R2Dtool annotate 
- path to output plot (including extension, e.g. .svg or .png)
- filter field; column used to select significant sites
- cutoff; integer value to use for determining significant
- "lower"/"upper": Select sites that have a significanace _lower_ or _higher_ than the cutoff value in the significance field 

Optional arguments: 
- The flag "-c [method]" can be used to change the strategy used for displaying confidence intervals between loess (default) or binoial confidence intervals (-c "binom")
- The flag "-o [/path/to/table.tsv]" can be used to save the aggregated metatranscript data as a tab-separated file 
- The flag "-l" displays transcript region labels (5' UTR, CDS, 3'UTR) on the plot (default = FALSE)
```



**Plot** the **metacodon** distribution of RNA features:

```
Rscript R2_plotMetaCodon.R [start/stop codon] [options] <input_file> <output_file> <field_name> <cutoff_value> <cutoff_type>

Mandatory positional arguments: 
- -s or -e to indicate start or stop codon 
- input_file: Path to the annotated transcriptomic sites file.
- output_file: Path where the plot will be saved (include file extension, e.g., .png or .svg).
- field_name: The name of the column in input_file used to filter significant sites.
- cutoff_value: Numeric value defining the threshold for significance.
- cutoff_type: Specifies the comparison direction, either "lower" or "upper", to determine significance.

Optional arguments: 
- The flag "-c [method]" can be used to change the strategy used for displaying confidence intervals between loess (default) or binoial confidence intervals (-c "binom")
- The flag "-o [/path/to/table.tsv]" can be used to save the aggregated metacodon data as a tab-separated file 
```
**Plot** the **metajunction** distribution of RNA features:

```
Rscript R2_plotMetaJunction.R [annotated transcriptomic sites] [save path for plot, including .png/.svg extension] [name of filter field] [cutoff for significant sites] [lower/upper cutoff for significant sites]

Mandatory positional arguments: 
- path to annotated file generated by R2Dtool annotate 
- path to output plot (including extension, e.g. .svg or .png)
- filter field; column used to select significant sites
- cutoff; integer value to use for determining significant
- "lower"/"upper": Select sites that have a significnace _lower_ or _higher_ than the cutoff value in the significance field 


Optional arguments: 
- The flag "-c [method]" can be used to change the strategy used for displaying confidence intervals between loess (default) or binoial confidence intervals (-c "binom")
- The flag "-o [/path/to/table.tsv]" can be used to save the aggregated metajunction data as a tab-separated file 
```


# Installation and dependencies 

#### Dependencies 

> [!NOTE]
> R2Dtool's core ```liftover``` and ```annotation``` functionality are built into an executable program built on Rust, and dependencies are managed during compilation by Cargo. 
> R2Dtool's ```metatranscript``` and ```metajunction``` plotting functinality are implemented as executable R scripts, and depend on the following packages:

```
r==4.2.0
tidyverse==1.3.1
genomicFeatures==1.48.0
rtracklayer==1.56.0
binom==binom_1.1
```
#### Compiling R2Dtool from source: 

```
# requires rust 
git clone git@github.com:comprna/R2Dtool.git && cd R2Dtool
cargo build --release

# add R2DTOOL to PATH
export PATH="$PATH:$(pwd)/target/release"
```

#### Testing R2Dtool installation 

```
cd R2Dtool && export PATH="$PATH:$(pwd)/target/release/"
mkdir ./test/outputs/ 2>/dev/null

# liftover transcriptomic sites to genomic coordinates
r2d liftover -H -g ./test/GRCm39_subset.gtf -i ./test/out_CHEUI_modelII.bed > ./test/outputs/liftover.bed

# annotate bed-like transcriptomic sites with metatranscript coordinates, distance to splice junctions, transcript structure and transcript biotype
r2d annotate -H -g ./test/GRCm39_subset.gtf -i ./test/out_CHEUI_modelII.bed > ./test/outputs/annotate.bed

# make metatranscript plot using annotated sites 
Rscript /home/150/as7425/R2Dtool/scripts/R2_plotMetaTranscript.R "./test/outputs/annotate.bed" "./test/outputs/metagene.png" "probability" "0.9999" "upper"

# make metajunction plot using annotated sites 
Rscript /home/150/as7425/R2Dtool/scripts/R2_plotMetaJunction.R "./test/outputs/annotate.bed" "./test/outputs/metagene.png" "probability" "0.9999" "upper"
```

> Test data was generated using [cheui](https://github.com/comprna/CHEUI) and converted to bed-like coordinates using [R2Dtool utilities](https://github.com/comprna/R2Dtool/blob/main/scripts/cheui_to_bed.sh).



# Input and output data structure

R2Dtool is designed to work with tab-delimited, plain text BED-like files (optionally) with a column names in row 1. Following UCSC's recommendations, input and output BED files should be tab-delimited plain-text:

- column 1 should represent transcript,
- column 2 and 3 should represent the coordinate of interest in zero-based, half open coordinates (e.g. the third nucleotide of a transcript is represented as (2,3),
- column 5, representing strand, is assumed as being (+)
- columns 4, and any column greater than 5 are preserved during annotation or liftover, and may be used to hold additional information, e.g. predictions of RNA modification stoichiometry, probability, etc.


# Gene annotation requirements 

R2Dtool is designed to work with GTF version 2 annotations

- Compatible feature types (col3) include 'transcript'/'mRNA', 'exon', 'CDS', 'five_prime_utr'/'5UTR'/'UTR' and 'three_prime_utr'/'3UTR'/'UTR'
- 'exon', 'transcript', and 'UTR' feature types are __required__ for liftover and annotation 
- Column 9 should optinally contain 'transcript_id', 'gene_id', 'gene_name', and some variety of biotype attribute (e.g. 'transcript_biotype', 'transcript_type', 'gene_type', or 'gene_biotype'). This data is used for R2Dtool's annotation function. 
- Transcripts that are duplicated between autosomal loci and PAR regions are assigned only to autosomal regions (_PAR_ gene entries are skipped by R2Dtool)

Critically, the GTF annotation provided __must__ contain identical gene structures used to generate the transcriptome, including identical transcript names in the FASTA header. One option is to use genomes, transcriptomes and gene structures from the same genome release. Another option is for users to generate their own transcriptome using a genome and gene structure file, e.g. using [gffread](https://ccb.jhu.edu/software/stringtie/gff.shtml#gffread_ex).  

# Utilities

### Convert CHEUI model II output to a BED-like input

- This script transposes CHEUI coordinates by +3 (BED interval start) and +4 (BED interval end) to represent a single nucleotides, and rearranges columns to become bed-like.

```
bash cheui_to_bed.sh [cheui model II output file] [cheui_to_bed output file]
```





```
       .=.                              
      '==c|
      [)-+|    <(RNA to DNA)
      //'_|        
 snd /]==;\                                                                                                                                                               
```









