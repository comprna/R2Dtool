# txannotate
Convert bed-like sites from transcriptomic to genomic coordinates using a genome-appropriate annotation

### Dependencies  
```
r==4.2.0
tidyverse==1.3.1
genomicFeatures==1.48.0
rtracklayer==1.56.0
```

### Annotate bed-like transcriptomic sites with metatranscript coordinates  transcript information

tbd


### Liftover bed-like transcriptomic sites from transcriptomic coordinates to genomic coordinates

tbd

### Utilities: Convert CHEUI model II output to a bed-like input
```
bash cheui_to_bed.sh [cheui model II output file] [cheui_to_bed output file]
```

### General notes
Transcriptome to genome liftover relies only on transcript_id, site coordinates, and strand (columns 1:3,6)
