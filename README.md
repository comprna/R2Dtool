# R2DTool

A utility for long-read isoform-centric epitranscriptomics that:
  - Annotates transcriptomic position with transcript-specific metatranscript coordinates and proximity to adjacent splice-junctions, and 
  - Transposes transcriptomic coordinates to genomic coordinates to enable the comparison of epitranscriptomic sites between overlapping transcript isoforms, and also to enable visualization of epitranscriptomic sites on a genome browser.

```
                                          ┌──────── liftover ───────►   Genome-centric transcriptomic sites
                                          │
Isoform-centric transcriptomic sites  ────┼──────── annotate ───────►   Annotated isoform-centric transcriptomic sites
                                          │
                                          └───── lift + annotate ───►   Annotated genome-centric transcriptomic sites
```

   - [Dependencies and installation](#dependencies-and-installation)
   - [Usage](#usage)
   - [Input and output data structure](#Input-and-output-data-structure)
   - [Utilities: Convert CHEUI model II output to a BED-like input](#utilities--convert-cheui-model-ii-output-to-a-bed-like-input)
   - [General notes](#general-notes)

------------------------------------------

## Dependencies and installation

R2DTool scripts require no installation beyond the availability of R and the listed packages:

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

Rscript ./scripts/annotate.R ./test/out_CHEUI_modelII.bed ./test/GRCm39_subset.gtf ./test/out_CHEUI_modelII_annotated.bed

# liftover annotated transcriptomic sites to genomic coordinates 

Rscript ./scripts/lift.R ./test/out_CHEUI_modelII_annotated.bed ./test/GRCm39_subset.gtf ./test/out_CMII_annotated_lifted.bed
```

## Usage 

### Annotating transcriptomic sites with metatranscript coordinates, splice junction distances, and gene structure information 

Input file format: bed-like file of transcriptomic sites (see input file structure for more information)
Output file format: Identical to input file, with n identical columns, representing ?what?

```
Rscript R2_annotate.R [bed-like transcriptomic sites] [gtf annotation] [output file]
```

### Liftover transcriptomic sites to genomic coordinates

Input file format: BED-like file of transcriptomic sites (col 1 = transcript)     
Output file format: BED-like file of identical sites in genomic coordinates (col 1 = chromosome/scaffold)      

```
bash txliftover.sh [BED-like transcriptomic sites] [GTF annotation] [BED-like output file]
```

Note: Liftover can be complete indpendantly of annotation.
Note: columns (1-3,6) represent the genomic positions of


## Input and output data structure

R2Dtool is designed to work with tab-delimited, plain text BED-like files with a header of column names in row 1. Following UCSC's recommendations, input and output BED files should be tab-delimited plain-text: 

- column 1 should represent transcript, 
- column 2 and 3 should represent the coordinate of interest in zero-based, half open coordinates (e.g. the third nucleotide of a transcript is represented as (2,3),
- column 5, representing strand, is assumed as being (+)
- columns 4, and any column greater than 5 are preserved during annotation or liftover, and may be used to hold additional information, e.g. predictions of RNA modification stoichiometry, probability, etc. 

The GTF annotation provided must contain identical gene structures used to generate the transcriptome reference from the genomic reference. One option is to use genomes, transcriptomes and gene structures from the same genome release. Another is to generate a transcriptome using a genome and gene structure file, e.g. using  


Note: 
- Because the output files will have a header, this may need to be removed prior to opening the file on some genome browsers. 


### Utilities: Convert CHEUI model II output to a BED-like input
- This script transposes CHEUI coordinates by +3 (BED interval start) and +4 (BED interval end) to represent a single nuecleotide
```
bash cheui_to_bed.sh [cheui model II output file] [cheui_to_bed output file]
```

### General notes
- Transcriptome to genome liftover relies only on transcript_id, site coordinates, and strand (columns 1:3,6)
- The annotation used must be in GTF format and correspond correctly to the transcriptome used
- Other columns aren't considered during annotation or liftover, and can be used to store additional information of interest
- A header containing column names is expected for annotation and liftover
- Issues with transcript_id and transcript_version? Try commenting out lines 32:33 and 40:41 in lift.R
- Using a GENCODE annotation? Try commenting out line 25 in lift.R and line 27 in annotate.R 
- GENCODE annotations use 'transcript_type' whereas Ensembl annotations use 'transcript_biotype'?

```
             ___
          ,-'___'-.
        ,'  [(_)]  '.
       |_]||[][O]o[][|
     _ |_____________| _
    | []   _______   [] |
    | []   _______   [] |
   [| ||      _      || |]
    |_|| =   [=]     ||_|
    | || =   [|]     || |
    | ||      _      || |
    | ||||   (+)    (|| |
    | ||_____________|| |
    |_| \___________/ |_|
    / \      | |      / \
   /___\    /___\SSt /___\

```
Artwork taken from  https://ascii.co.uk/art/starwars1
