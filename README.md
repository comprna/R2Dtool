# txannotate
Convert bed-like sites from transcriptomic to genomic coordinates using a genome-appropriate annotation

### dependencies  
```
r==4.2.0
r::tidyverse==1.3.1
r::genomicFeatures==1.48.0
r::rtracklayer==1.56.0
```

### convert CHEUI model II output to a bed-like input
```
bash cheui_to_bed.sh [cheui model II output file] [cheui_to_bed output file]
```

### notes
Transcriptome to genome liftover relies only on transcript_id, site coordinates, and strand (columns 1:3,6)
