# R2Dtool

R2Dtool is a set of genomics utilities for handling, integrating, and viualising isoform-mapped RNA feature data.

  - **Integrate** transcriptome-mapped data: R2Dtool performs liftover of transcriptome-mapped RNA features to their corresponding genomic coordinates

  - **Annotate** transcriptome-mapped sites: R2Dtool annotates transcript-specific metatranscript coordinates and absolute and relative distances to annotated transcript landmarks, in an isoform-specific manner.

  - **Visualise** isoform-aware RNA feature distributions: R2Dtool introduces isoform-aware **metatranscript** plots and **metajunction** plots to study the positonal distribution of RNA features around annotated RNA landmarks.


> [!NOTE]
> [Link to R2Dtool preprint on biorXiv](https://doi.org/10.1101/2022.09.23.509222)    
> Sethi, A. J., Mateos, P. A., Hayashi, R., Shirokikh, N. & Eyras, E. R2Dtool: Positional Interpretation of RNA-Centric Information in the Context of Transcriptomic and Genomic Features. biorXiv 2022.09.23.509222 (2022)
> doi:10.1101/2022.09.23.509222.

# Contents 

   - [Quick-start guide](#quick-start-guide)
       - [Requirements](#requirements)
       - [Compiling from source](#compiling-from-source)
       - [R2Dtool tutorial](#r2dtool-tutorial)

   - [General usage](#usage)
       - [Handling isoform-mapped RNA sites](#handling-isoform-mapped-rna-sites)
       - [Visualising isoform-aware RNA feature metadistributions](#visualising-isoform-aware-rna-feature-metadistributions)

   - [Input and output data structure](#input-and-output-data-structure)
   - [Gene annotation requirements](#gene-annotation-rquirements)

# Quick-start guide 

#### Requirements

R2Dtool has been tested on Mac OS X 14.1, and Red Hat Enterprise Linux 4.18.0, and should be broadly compatible with any system that support Rustc and R. 

- R2Dtool requires ```rustc``` and ```cargo``` for compilation. 
- Additionally, ```R```, along with ```tidyverse``` and ```binom``` R packages, are requred to generate metaplots. 

#### Compiling from source: 

```
# clone R2Dtool and compile using Cargo 
git clone git@github.com:comprna/R2Dtool.git && cd R2Dtool
cargo build --release
```

#### R2Dtool tutorial

To test the R2Dtool installation, we provide a toy dataset of isoform-mapped m6A predictions in HeLa cells, taken from [Figure 1E and 1F of the manuscript.](https://doi.org/10.1101/2022.09.23.509222) 

The input dataset, ```m6A_isoform_sites_GRCh38_subset.bed```, contains the positions and stoichiometry of m6A sites mapped to individual transcript isoforms, in a tab-separated ```BED3+``` format:  
The input file has headers, so we will pass the '-H' flag to R2Dtool to specify that headers are present. 

```
$ head -n 5 ./test/m6A_isoform_sites_GRCh38_subset.bed

transcript      start   end     base    coverage        strand  N_valid_cov     fraction_modified
ENST00000381989.4       2682    2683    a       10      +       10      0.00
ENST00000381989.4       2744    2745    a       10      +       10      0.00
ENST00000381989.4       2928    2929    a       10      +       10      0.00
ENST00000381989.4       2985    2986    a       10      +       10      0.00
```

R2Dtool allows us to liftover these sites to genomic coordinates, calculate the distances of m6A sites to annotated gene features, and visualise the isoform-resolved distributions of RNA feature around genomic landmarks. 


In order to compare the positions of our m6A sites to annotated genomic features, or to view their genomic positions in the genome browser, we can use ```r2d liftover```. 
Liftover adds the genomic coordinates of transcriptome positions to the input bed file. 

```
# add R2Dtool to PATH
$ cd R2Dtool && export PATH="$PATH:$(pwd)/target/release/"

# liftover transcriptomic sites to genomic coordinates
$ r2d liftover -H -g ./test/GRCh38.110_subset.gtf -i ./test/m6A_isoform_sites_GRCh38_subset.bed > ./test/liftover.bed

$ head -n 5 ./test/liftover.bed
chromosome      start   end     name    score   strand  transcript      start   end     base    coverage        strand  N_valid_cov     fraction_modified
13      24455165        24455166                        -       ENST00000381989.4       2682    2683    a       10      +       10      0.00
13      24455103        24455104                        -       ENST00000381989.4       2744    2745    a       10      +       10      0.00
13      24452564        24452565                        -       ENST00000381989.4       2928    2929    a       10      +       10      0.00
13      24452507        24452508                        -       ENST00000381989.4       2985    2986    a       10      +       10      0.00
```
The genomic cooridnates have been added to columns 1-6. We can now open liftover.bed in the genome browser, or compare the positions of our m6A sites to features previously annotated in genomic coordinates, e.g. using [Bedtools](https://bedtools.readthedocs.io/en/latest/). 

We can also use R2Dtool to annotate the position of each m6A site, assigning a metatranscript region to the site, and calculating distances to local transcript landmarks. 
This can be accomplished using R2Dtool's annotate function: 

```
# annotate bed-like transcriptomic sites with metatranscript coordinates, distance to splice junctions, transcript structure and transcript biotype
$ r2d annotate -H -g ./test/GRCh38.110_subset.gtf -i ./test/m6A_isoform_sites_GRCh38_subset.bed > ./test/annotate.bed

$ head -n 5  ./test/annotate.bed 
transcript      start   end     base    coverage        strand  N_valid_cov     fraction_modified       gene_id gene_name       transcript_biotype      tx_len  cds_start       cds_end tx_end  transcript_metacoordinate       abs_cds_start   abs_cds_end     up_junc_dist    down_junc_dist
ENST00000381989.4       2682    2683    a       10      +       10      0.00    ENSG00000102699 PARP4   protein_coding  5437    74      5246    5437    1.50425 2608    -2564   46      150
ENST00000381989.4       2744    2745    a       10      +       10      0.00    ENSG00000102699 PARP4   protein_coding  5437    74      5246    5437    1.51624 2670    -2502   108     88
ENST00000381989.4       2928    2929    a       10      +       10      0.00    ENSG00000102699 PARP4   protein_coding  5437    74      5246    5437    1.55182 2854    -2318   28      160
ENST00000381989.4       2985    2986    a       10      +       10      0.00    ENSG00000102699 PARP4   protein_coding  5437    74      5246    5437    1.56284 2911    -2261   85      103
```

If we want to understand the distribution of these features across the length of transcripts, we can make metatranscript plots. 

In this plot, we will show the proportion of sites that have more than 10% m6A methylation, across a metatrancscript model: 

```
# make metatranscript plot using annotated sites 
$ r2d plotMetaTranscript -i "./test/annotate.bed" -o "./test/metatranscript_m6A.png" -f "fraction_modified" -u "10" -t "upper" -l 
# `-f "fraction_modified" -u "10" -t "upper"` means we only consider sites with greater than "10" in the "fraction_modified" column to be significantly m6A methylated 
# `-l` flag adds labels for metatrancript regions 
```

We can also study the distribution of m6A sites around exon-exon junctions in an isoform specific manner, using the plotMetaJunction function: 

```
# make metajunction plot using annotated sites 
$ r2d plotMetaJunction -i "./test/annotate.bed" -o "./test/metajunction_m6A.png" -f "fraction_modified" -u "10" -t "upper" 
# similarly to before, we only consider sites with more than 10% methylatin stoichiometry in the fraction modified column to be significantly methylated 
```

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

[More information on ```r2d liftover`` can be found on the R2Dtool wiki](https://github.com/comprna/R2Dtool/wiki/Further-information-on-r2d-liftover)

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

[More information on ```r2d annotate``` can be found on the R2Dtool wiki](https://github.com/comprna/R2Dtool/wiki/Further-information-on-r2d-annotate)

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

[More information on ```r2d``` plot functions can be found on the R2Dtool wiki](github.com/comprna/R2Dtool/wiki/Visualising-RNA-feature-distributions-with-R2Dtool)

# Input and output data structure

R2Dtool is designed to work with tab-delimited, plain text BED-like files (optionally) with a column names in row 1. Following UCSC's recommendations, input and output BED files should be tab-delimited plain-text:

- column 1 should represent transcript,
- column 2 and 3 should represent the coordinate of interest in zero-based, half open coordinates (e.g. the third nucleotide of a transcript is represented as (2,3),
- column 5, representing strand, is assumed as being (+)
- columns 4, and any column greater than 5 are preserved during annotation or liftover, and may be used to hold additional information, e.g. predictions of RNA modification stoichiometry, probability, etc.

An example of appropriate input data is availale in the ```./test/m6A_isoform_sites_GRCh38_subset.bed``` directory. 

Notes: 
* In general, most upstream methods for producing isoform-resolved RNA feature maps include header headers in col1, but many downstream tools, such as genome browsers and bedtools do not expecte headers. 
* R2Dtool ```liftover``` and ```annotate``` can take input files without headers (default), or with headers (-H flag)
* We recommend using input files with headers for R2Dtools analysis
* ```r2d annotate``` adds headers to the output file, and these are required for plotting. 

# Gene annotation requirements 

R2Dtool is designed to work with GTF version 2 annotations

- Compatible feature types (col3) include 'transcript'/'mRNA', 'exon', 'CDS', 'five_prime_utr'/'5UTR'/'UTR' and 'three_prime_utr'/'3UTR'/'UTR'
- 'exon', 'transcript', and 'UTR' feature types are __required__ for liftover and annotation 
- Column 9 should optinally contain 'transcript_id', 'gene_id', 'gene_name', and some variety of biotype attribute (e.g. 'transcript_biotype', 'transcript_type', 'gene_type', or 'gene_biotype'). This data is used for R2Dtool's annotation function. 
- Transcripts that are duplicated between autosomal loci and PAR regions are assigned only to autosomal regions (_PAR_ gene entries are skipped by R2Dtool)

Critically, the GTF annotation provided __must__ contain identical gene structures used to generate the transcriptome, including identical transcript names in the FASTA header. One option is to use genomes, transcriptomes and gene structures from the same genome release. Another option is for users to generate their own transcriptome using a genome and gene structure file, e.g. using [gffread](https://ccb.jhu.edu/software/stringtie/gff.shtml#gffread_ex).  

[More information about selecting the right gene annotation is available through the R2Dtool Wiki](https://github.com/comprna/R2Dtool/wiki/R2Dtool-GTF-requirements)



```
       .=.                              
      '==c|
      [)-+|    <(RNA to DNA)
      //'_|        
 snd /]==;\                                                                                                                                                               
```









