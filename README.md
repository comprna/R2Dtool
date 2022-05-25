# txannotate
Convert bed-like sites from transcriptomic to genomic coordinates using a genome-appropriate annotation

### Installation and dependencies  

txannotate scripts require no installation beyond the availability of R and the listed packages:

Installation:
```
git clone git@github.com:comprna/txannotate.git
```

Dependencies:
```
r==4.2.0
tidyverse==1.3.1
genomicFeatures==1.48.0
rtracklayer==1.56.0
```

### Annotate transcriptomic sites with metatranscript coordinates and transcript information

Input file format: bed-like file of transcriptomic sites (col 1 = transcript)
Output file format: Input, with n additional columns corresponding to ...

```
bash txannotate.sh [bed-like transcriptomic sites] [output file]
```

### Liftover transcriptomic sites to genomic coordinates

Input file format: bed-like file of transcriptomic sites (col 1 = transcript)
Output file format: bed-like file of identical sites in genomic coordinates (col 1 = chromosome/scaffold)

```
bash txliftover.sh [bed-like transcriptomic sites] [transcriptome-specific GTF annotation] [output file]
```

### Utilities: Convert CHEUI model II output to a bed-like input
```
bash cheui_to_bed.sh [cheui model II output file] [cheui_to_bed output file]
```

### General notes
- Transcriptome to genome liftover relies only on transcript_id, site coordinates, and strand (columns 1:3,6)
- Other columns aren't considered during annotation or liftover, and can be used to store additional information of interest
- A header containing column names is expected for annotation and liftover
