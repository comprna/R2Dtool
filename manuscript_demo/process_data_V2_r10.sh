#!/bin/bash 

# minimap v2.1
export PATH="$PATH:/g/data/lf10/as7425/apps/minimap2/"
module load samtools

# start with r10 promethion of hela transcriptome
export wd="/g/data/lf10/as7425/R2DTool_demo"; mkdir -p ${wd} 2>/dev/null
export reads="/g/data/lf10/as7425/2023_mrna-biogenesis-maps/analysis/2024-03-28_ASR012_HeLa-total-RNA004-rep1/ASR012_HeLa-total-RNA004-rep1_dorado_polyA_m6A_all_reads.fastq"
export genome="/g/data/lf10/as7425/genomes/human_genome/ensembl_release_110/Homo_sapiens.GRCh38.cdna.all.fa"

# map to genome
minimap2 -ax map-ont -y -k 14 -N 10 ${genome} ${reads} -t 104 | samtools view -bh > ${wd}/aligned_reads.bam

# filter for primary alignments with nonzero mapq
samtools view -b -F 2308 -@ 104 -q 1 ${wd}/aligned_reads.bam | samtools sort > ${wd}/primary.bam
samtools view -@ 104 -c ${wd}/primary.bam
samtools index ${wd}/primary.bam

# pileup mod calls
export PATH="$PATH:/g/data/lf10/as7425/apps/modkit/target/release"
modkit pileup -t 104 ${wd}/primary.bam "${wd}/pileup.bed" --log-filepath "${wd}/pileup.bed.log"

# prepare for R2Dtool 
# Prepare for R2Dtool
printf "transcript\tstart\tend\tmod\tcov\tstrand\tstart\tend\tcolor\tX1\tX2\tcov\tfrac_mod\tn_mod\tn_canonical\tn_other_mod\tn_delete\tn_fail\tn_diff\tn_no_call\n" > ${wd}/R2D_input.bed
awk 'BEGIN {FS=OFS="\t"} {gsub(",", OFS, $0); $9="."; print}' "${wd}/pileup.bed" | tr " " "\t" | awk '($12 > 9)' >> ${wd}/R2D_input.bed

# check that all mapped to 'A' nucleotide
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




module load R 

time Rscript ~/R2Dtool/scripts/R2_plotMetaTranscript.R "${wd}/methylation_calls_annotated.bed" "${wd}/metatranscript_out_Rust.png" "frac_mod" "9" "upper" -c loess -l -o "${wd}/transcript_data_Rust.tsv"

time Rscript ~/R2Dtool/scripts/R2_plotMetaJunction.R "${wd}/methylation_calls_annotated.bed" "${wd}/metajunction_out_binom_rust.png" "frac_mod" "9" "upper" -c loess -o "$wd/table"


##################################################
##################################################



# filter for high-quality mappings
samtools view -@ 104 -b -F 2308 ${wd}/all_alignments.sam > ${wd}/primary_alignemnts.bam
samtools view -c -@ 104 ${wd}/primary_alignemnts.bam
# 2767236 aln

# filter for assignments to single isoform 
samtools view -@ 104 -b -q 60 ${wd}/primary_alignemnts.bam | samtools sort > ${wd}/isoform_assignments.bam
samtools view -@ 104 -c ${wd}/isoform_assignments.bam
samtools index ${wd}/isoform_assignments.bam
# 2051031 aln

