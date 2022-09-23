#!/bin/bash

# Written by AJ Sethi on 2022-09-08
# Last modified on ...

##################################################

# add files to environment
export methylationCalls="/g/data/xc17/pm1122/xpore_hek293/results/HEK293T-WT_CHEUI_predictions_site_level_WT.txt"
export annotation="/g/data/lf10/as7425/genomes/human_genome/gencode/gencode_v38/gencode.v38.annotation.gtf"
export wd="/scratch/lf10/as7425/R2_example"; mkdir -p ${wd}

# convert cheui data to bed-like
bash ~/R2Dtool/scripts/cheui_to_bed.sh ${methylationCalls} "${wd}/methylationCalls.bed"

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

# plotMetaTranscript
time Rscript ~/R2Dtool/scripts/R2_plotMetaTranscript.R "~/localGadiData/2022-09-21_R2Dtool/methylationCalls_annotated_lifted.bed" "~/localGadiData/2022-09-21_R2Dtool/out.png" "probability" "0.9999" "upper"

# plotMetaJunction
time Rscript ~/R2Dtool/scripts/R2_plotMetaJunction.R "~/localGadiData/2022-09-21_R2Dtool/methylationCalls_annotated_lifted.bed" "~/localGadiData/2022-09-21_R2Dtool/out.png" "probability" "0.9999" "upper"
