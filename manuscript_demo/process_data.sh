#!/bin/bash 

# using minimap v2.1 and samtools v1.19
export PATH="$PATH:/g/data/lf10/as7425/apps/minimap2/"
module load samtools

# map hela DRS reads to GRCh38 cDNA transcriptome
export wd="/g/data/lf10/as7425/R2DTool_demo"; mkdir -p ${wd} 2>/dev/null
export reads="/g/data/lf10/as7425/2023_mrna-biogenesis-maps/analysis/2024-03-28_ASR012_HeLa-total-RNA004-rep1/ASR012_HeLa-total-RNA004-rep1_dorado_polyA_m6A_all_reads.fastq"
export genome="/g/data/lf10/as7425/genomes/human_genome/ensembl_release_110/Homo_sapiens.GRCh38.cdna.all.fa"

minimap2 -ax map-ont -y -k 14 -N 20 ${genome} ${reads} -t 104 | samtools view -bh > ${wd}/aligned_reads.bam

# filter for plus strand primary alignments with nonzero mapq
samtools view -b -F 2308 -@ 104 -q 1 ${wd}/aligned_reads.bam | samtools sort > ${wd}/primary.bam
samtools view -@ 104 -c ${wd}/primary.bam
samtools index ${wd}/primary.bam

# pileup the modification calls using modkit 
export PATH="$PATH:/g/data/lf10/as7425/apps/modkit/target/release"
modkit pileup -t 104 ${wd}/primary.bam "${wd}/pileup.bed" --log-filepath "${wd}/pileup.bed.log"

# convert to table with header for R2Dtool 
# filter for coverage >=10
printf "transcript\tstart\tend\tmod\tcov\tstrand\tstart\tend\tcolor\tX1\tX2\tcov\tfrac_mod\tn_mod\tn_canonical\tn_other_mod\tn_delete\tn_fail\tn_diff\tn_no_call\n" > ${wd}/R2D_input.bed
awk 'BEGIN {FS=OFS="\t"} {gsub(",", OFS, $0); $9="."; print}' "${wd}/pileup.bed" | tr " " "\t" | awk '($12 > 9)' >> ${wd}/R2D_input.bed

# verify that all sites mapped to 'A' nucleotide as expected
export wd="/g/data/lf10/as7425/R2DTool_demo"; mkdir -p ${wd} 2>/dev/null
export genome="/g/data/lf10/as7425/genomes/human_genome/ensembl_release_110/Homo_sapiens.GRCh38.cdna.all.fa"
tail -n +2 ${wd}/R2D_input.bed > ${wd}/R2D_input_noheader.bed
bedtools getfasta -s -fi ${genome} -bed ${wd}/R2D_input_noheader.bed | tail -n +2 | awk 'NR%2==1' | sort | uniq -c | sort -nr > ${wd}/input_sequence_context.txt
cat ${wd}/input_sequence_context.txt

# 336160 A
#     782 T
#     333 G
#     300 C


##################################################

# run R2Dtool 

export wd="/g/data/lf10/as7425/R2DTool_demo";
export sites="${wd}/R2D_input.bed"
export anno="/g/data/lf10/as7425/genomes/human_genome/ensembl_release_110/Homo_sapiens.GRCh38.110.chr_patch_hapl_scaff.gtf"

# R2Dtool rust implementation 
cd ~/R2Dtool && rm -rf target; cargo build --release 
export PATH="${PATH}:/home/150/as7425/R2Dtool/target/release/"

time r2d annotate -i "${sites}" -g ${anno} -H > "${wd}/methylation_calls_annotated.bed"
time r2d liftover -i "${sites}" -g ${anno} -H > "${wd}/methylation_calls_annotated_liftover.bed"

##################################################

# check sequence context in the liftover 

export genome="/g/data/lf10/as7425/genomes/human_genome/ensembl_release_110/Homo_sapiens.GRCh38.dna_sm.toplevel.fa"
tail -n +2 "${wd}/methylation_calls_annotated_liftover.bed" | cut -f1-6 | awk -F "\t" '($6 == "-")' > "${wd}/methylation_calls_annotated_liftover_noheader.bed" 
bedtools getfasta -s -fi ${genome} -bed "${wd}/methylation_calls_annotated_liftover_noheader.bed" | tail -n +2 | awk 'NR%2==1' | sort | uniq -c | sort -nr > ${wd}/liftover_sequence_context.txt
cat ${wd}/liftover_sequence_context.txt



##################################################

# make R2Dtool plots 

module load R 
time Rscript ~/R2Dtool/scripts/R2_plotMetaTranscript.R "${wd}/methylation_calls_annotated.bed" "${wd}/metatranscript_out_Rust.svg" "frac_mod" "9" "upper" -c loess -l -o "${wd}/transcript_data_Rust.tsv"
time Rscript ~/R2Dtool/scripts/R2_plotMetaJunction.R "${wd}/methylation_calls_annotated.bed" "${wd}/metajunction_out_binom_rust.svg" "frac_mod" "9" "upper" -c loess -o "$wd/table"


##################################################
