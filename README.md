# txannotate

A utility for manipulating transcriptomic coordinates for epitranscriptomic applications. 
Anotate transcriptomic sites with metatranscript coordinates, calculate distance of a transcriptomic site from the nearest upstream and downstream splice sites, and lift-over transcriptomic sites to to genomic coordinates to enable visualisation of epitranscriptomic sites on a genome browser.

### Installation and dependencies  

txannotate scripts require no installation beyond the availability of R and the listed packages:

Dependencies:
```
r==4.2.0
tidyverse==1.3.1
genomicFeatures==1.48.0
rtracklayer==1.56.0
```

Installation:
```
git clone git@github.com:comprna/txannotate.git
```

Testing the installation:
```
cd txannotate

# convert CHEUI model II output to BED

bash cheui_to_bed.sh ./test/CHEUI_modelII_subset.txt ./test/out_CHEUI_modelII.bed

# annotate bed-like transcriptomic sites with metatranscript coordinates, distance to splice junctions, transcript structure and biotype 

Rscript annotate.R ./test/out_CHEUI_modelII.bed ./test/GRCm39_subset.gtf ./test/out_CHEUI_modelII_annotated.bed

# liftover annotated transcriptomic sites

Rscript lift.R ./test/out_CHEUI_modelII_annotated.bed ./test/GRCm39_subset.gtf ./test/out_CMII_annotated_lifted.bed
```

### Annotate transcriptomic sites with metatranscript coordinates and transcript information

Input file format: bed-like file of transcriptomic sites (col 1 = transcript).   
Output file format: Input, with n additional columns corresponding to ...    

Note: txannotate.sh still needs to be adapted from cheui-liftover.sh
```
bash txannotate.sh [bed-like transcriptomic sites] [output file]
```

### Liftover transcriptomic sites to genomic coordinates

Input file format: bed-like file of transcriptomic sites (col 1 = transcript)     
Output file format: bed-like file of identical sites in genomic coordinates (col 1 = chromosome/scaffold)      

```
bash txliftover.sh [bed-like transcriptomic sites] [annotation] [output file]
```

### Utilities: Convert CHEUI model II output to a bed-like input
- This script transposes CHEUI coordinates by +3 (bed interval start) and +4 (bed interval end) to represent a single nuecleotide
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
