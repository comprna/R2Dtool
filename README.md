# txannotate
Annotate transcriptomic sites with metatranscript coordinates and convert to genomic coordinates using a genome annotation

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

bash txliftover.sh ./test/test_mouse_GRCm39_cheui-model-ii.txt \
./test/subset_GRCm39.104_annotation.gtf \
./test/test_mouse_GRCm39_cheui-model-ii.txt \
./test/liftover_output.txt
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
```
bash cheui_to_bed.sh [cheui model II output file] [cheui_to_bed output file]
```

### General notes
- Transcriptome to genome liftover relies only on transcript_id, site coordinates, and strand (columns 1:3,6)
- The annotation used must be in GTF format and correspond correctly to the transcriptome used
- Other columns aren't considered during annotation or liftover, and can be used to store additional information of interest
- A header containing column names is expected for annotation and liftover
- Issues with transcript_id and transcript_version? Try commenting out lines 32:33 and 40:41 in lift.R
