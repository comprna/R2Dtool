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
cat <(cat ${cheui_model_ii_output} | head -n 1) <(cat ${cheui_model_ii_output} |\
tail -n +2 |\
tr "_" "\t" |\
awk '{ print $1"\t"$2+3"\t"$2+4"\t"$2+3";"$3";"$4";"$5";"$6"\t.\t+"}') > ${bedlike_cheui_model_ii_output}
```

### notes
Transcriptome to genome liftover relies only on transcript_id, site coordinates, and strand (columns 1:3,6)
