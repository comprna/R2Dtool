#!/bin/bash

# Written by AJ Sethi on 2022-09-08
# Last modified on 2024-04-09 using R2Dtool v2.0.0 

##################################################

# R2Dtool rust implementation 
cd ~/R2Dtool && cargo build --release 
export PATH="${PATH}:/home/150/as7425/R2Dtool/target/release/"

# Input data corresponds to HEK293 m6A predictions, calculated against the gencode V38 transcriptome 
export methylation_calls="/g/data/xc17/pm1122/xpore_hek293/results/HEK293T-WT_CHEUI_predictions_site_level_WT.txt"
export annotation="/g/data/lf10/as7425/genomes/human_genome/gencode/gencode_v38/gencode.v38.annotation.gtf"
export wd="/g/data/lf10/as7425/R2DTool_demo/"; mkdir -p ${wd}

# convert cheui data to bed-like
bash ~/R2Dtool/scripts/cheui_to_bed.sh ${methylation_calls} "${wd}/methylation_calls.bed"

# annotate the methylation calls against gencode v38 
# temporarty 
cd ~/R2Dtool && rm -rf target && cargo build --release 
annotation="/home/150/as7425/R2Dtool/test/gencode_v38.gtf"
cd $wd; rm splice_sites_map.txt 2>/dev/null
time r2d annotate -i "${wd}/methylation_calls.bed" -g ${annotation} -H > "${wd}/methylation_calls_annotated.bed"
cat splice_sites_map.txt | grep "ENST00000000233" -A 20 -B 20 

# compare to Rscript version
module load R 
# time Rscript ~/R2Dtool/scripts/R2_annotate.R "${wd}/methylation_calls.bed" "${annotation}" "${wd}/methylation_calls_annotated_R.bed"
cd $wd
for i in *; do head $i > ~/toLocal/${i##*/}; done


# liftover the annotated called to genomic coordinates 
time r2d liftover -i "${wd}/methylation_calls_annotated.bed" -g ${annotation} -H > "${wd}/methylation_calls_annotated_lifted.bed"

# plotMetaTranscript
module load R 
time Rscript ~/R2Dtool/scripts/R2_plotMetaTranscript.R "${wd}/methylation_calls_annotated_lifted.bed" "${wd}/metatranscript_out.png" "probability" "0.9999" "upper"

# plotMetaJunction
time Rscript ~/R2Dtool/scripts/R2_plotMetaJunction.R "${wd}/methylation_calls_annotated_lifted.bed" "${wd}/metajunction_out.png" "probability" "0.9999" "upper"


# annotate methylation calls using Gencode v38


time Rscript ~/R2Dtool/scripts/R2_annotate.R "${wd}/methylationCalls.bed" "${annotation}" "${wd}/methylationCalls_annotated.bed"
# takes about 9 minutes

# lift methylation calls to genomic coordinates using GENCODE v38
time Rscript ~/R2Dtool/scripts/R2_lift.R "${wd}/methylationCalls_annotated.bed" "${annotation}" "${wd}/methylationCalls_annotated_lifted.bed"
# takes about 5 minutes

# filter for significant sites
cat <(head -n 1 "${wd}/methylationCalls_annotated_lifted.bed") <(awk '($11>0.9999)' "${wd}/methylationCalls_annotated_lifted.bed") > "${wd}/methylationCalls_annotated_lifted_significant.bed"

# compress
cat "${wd}/methylationCalls_annotated_lifted_significant.bed" | gzip -c > "${wd}/methylationCalls_annotated_lifted_significant.bed.gz"

