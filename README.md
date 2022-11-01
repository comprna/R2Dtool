# R2Dtool

A utility for long-read isoform-centric epitranscriptomics that:
  - Annotates (epi)transcriptomic positions with transcript-specific metatranscript coordinates and proximity to adjacent splice-junctions, and
  - Transposes transcriptomic coordinates to their underlying genomic coordinates to enable the comparison of epitranscriptomic sites between overlapping transcript isoforms, and also to enable visualization of epitranscriptomic sites on a genome browser.

```
                                          ┌──────── liftover ───────►   Genome-centric transcriptomic sites
                                          │
Isoform-centric transcriptomic sites  ────┼──────── annotate ───────►   Annotated isoform-centric transcriptomic sites
                                          │
       .=.                                └───── lift + annotate ───►   Annotated genome-centric transcriptomic sites
      '==c|
      [)-+|    <(RNA to DNA)
      //'_|        
 snd /]==;\                                                                                                                                                                     
```

   - [Dependencies and installation](#dependencies-and-installation)
   - [Usage](#usage)
   - [Input and output data structure](#Input-and-output-data-structure)
   - [Utilities: Convert CHEUI model II output to a BED-like input](#utilities)
   - [General notes](#general-notes)

------------------------------------------

## Dependencies and installation

R2Dtool scripts require no installation beyond the availability of R and the listed packages:

Dependencies:
```
r==4.2.0
tidyverse==1.3.1
genomicFeatures==1.48.0
rtracklayer==1.56.0
```

Downloading and testing R2Dtool:

```
# download R2Dtool from github

git clone git@github.com:comprna/R2Dtool.git
cd R2Dtool

# annotate bed-like transcriptomic sites with metatranscript coordinates, distance to splice junctions, transcript structure and transcript biotype

Rscript ./scripts/R2_annotate.R ./test/out_CHEUI_modelII.bed ./test/GRCm39_subset.gtf ./test/out_CHEUI_modelII_annotated.bed

# liftover annotated transcriptomic sites to genomic coordinates

Rscript ./scripts/R2_lift.R ./test/out_CHEUI_modelII_annotated.bed ./test/GRCm39_subset.gtf ./test/out_CMII_annotated_lifted.bed
```

Note: Test data was generated using [cheui](https://github.com/comprna/CHEUI) and converted to bed-like coordinates using [R2Dtool utilities](https://github.com/comprna/R2Dtool/blob/main/scripts/cheui_to_bed.sh).

## Usage

### Annotating transcriptomic sites (metatranscript coordinates, splice junction distances, gene structure)

Annotation adds the following information to the epitranscriptomic sites as additional coluumns, relying on the gene structure GTF to generate these data.

```
transcript_biotype | gene_name | gene_id | tx_len | cds_len | utr5_len | utr3_len | cds_start | cds_end | tx_end  | rel_pos | abs_cds_start | abs_cds_end | up_junc_dist | down_junc_dist
```

- *tx_len* represents the transcript length, and *cds-*, *utr5-* and *utr3-len* represent the length of the coding sequence, 5' UTR and 3' UTR respectively.
- *cds_start* and *cds_end* represent the positions of the coding sequence start and end compared to the transcript.
- *rel_pos* represents the scaled metatrascript position of the epitranscriptomic site, between 0 and 3, where 0 represents TSS, 1 represents CDS start, 2 represent CDS end, and 3 represents p(A) site.
- *abs_cds_start* and *abs_cds_end* represent the absolute distance (in nt) of a given site from the cds start and end
- *up_junc_dist* and *down_junc_dist* repreesnt the absolute distance (in nt) of a given site from the nearest upstream and downstream splice-junction (on that specific transcript).

```
Rscript R2_annotate.R [bed-like transcriptomic sites] [gtf annotation] [annotated bed-like output in transcriptomic coordinates]
```

### Liftover transcriptomic sites to genomic coordinates

Liftover converts a bed-like file of transcriptomic sites from transcriptomic to genomic coordinates.

```
Rscript R2_lift.R [BED-like transcriptomic sites] [GTF annotation] [annotated bed-like output in genomic coordinates]
```

Note: Liftover can be complete independantly of annotation, or following annotation.

### Make a metatranscript plot

```
Rscript R2_plotMetaTranscript.R [annotated transcriptomic sites] [save path for plot, including .png/.svg extension] [name of filter field] [cutoff for significant sites] [lower/upper cutoff for significant sites]
```

Notes:
- *filter field* corresponds to the column name by which significant sites are selected
- *cutoff* is the value of the filter field at which significant sites are chosen
- *lower/upper* allows the selection of significant sites which have a value either *lower* or *higher* than the cutoff.
- It is recommended that users use the plotMetaTranscript code as a starting point to generate publicaion quality metaTranscript plot. The code is annotated and easily edited.

### Make a metajunction plot

```
Rscript R2_plotMetaJunction.R [annotated transcriptomic sites] [save path for plot, including .png/.svg extension] [name of filter field] [cutoff for significant sites] [lower/upper cutoff for significant sites]
```

Notes:
- *filter field* corresponds to the column name by which significant sites are selected
- *cutoff* is the value of the filter field at which significant sites are chosen
- *lower/upper* allows the selection of significant sites which have a value either *lower* or *higher* than the cutoff.
- It is recommended that users use the plotMetaTranscript code as a starting point to generate publicaion quality metaTranscript plot. The code is annotated and easily edited.

## Input and output data structure

R2Dtool is designed to work with tab-delimited, plain text BED-like files with a header of column names in row 1. Following UCSC's recommendations, input and output BED files should be tab-delimited plain-text:

- column 1 should represent transcript,
- column 2 and 3 should represent the coordinate of interest in zero-based, half open coordinates (e.g. the third nucleotide of a transcript is represented as (2,3),
- column 5, representing strand, is assumed as being (+)
- columns 4, and any column greater than 5 are preserved during annotation or liftover, and may be used to hold additional information, e.g. predictions of RNA modification stoichiometry, probability, etc.

The GTF annotation provided must contain identical gene structures used to generate the transcriptome, including identical transcript names in the FASTA header. One option is to use genomes, transcriptomes and gene structures from the same genome release. Another option is for users to generate their own transcriptome using a genome and gene structure file, e.g. using gffread.  

## Utilities
### Convert CHEUI model II output to a BED-like input

- This script transposes CHEUI coordinates by +3 (BED interval start) and +4 (BED interval end) to represent a single nucleotides, and rearranges columns to become bed-like.

```
bash cheui_to_bed.sh [cheui model II output file] [cheui_to_bed output file]
```

## General notes
- Transcriptome to genome liftover relies only on transcript_id, site coordinates, and strand (columns 1:3,6)
- The annotation used must be in GTF format and correspond correctly to the transcriptome used
- Other columns aren't considered during annotation or liftover, and can be used to store additional information of interest
- A header containing column names is expected for annotation and liftover
- Issues with transcript_id and transcript_version? Try commenting out lines 32:33 and 40:41 in lift.R
- Using a GENCODE annotation? Try commenting out line 25 in lift.R and line 27 in annotate.R
- GENCODE annotations use 'transcript_type' whereas Ensembl annotations use 'transcript_biotype'?
- Because the output files will have a header, this may need to be removed prior to opening the file on some genome browsers.
