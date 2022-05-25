# txannotate
Convert bed-like sites from transcriptomic to genomic coordinates using a genome-appropriate annotation


### convert CHEUI model II output to a bed-like input
```
cat <(cat ${cheui_model_ii_output} | head -n 1) <(cat ${cheui_model_ii_output} |\
tail -n +2 |\
tr "_" "\t" |\
awk '{ print $1"\t"$2+3"\t"$2+4"\t"$2+3";"$3";"$4";"$5";"$6"\t.\t+"}') > ${bedlike_cheui_model_ii_output}
```
