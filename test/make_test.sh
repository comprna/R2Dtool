# navigate to cheui data directory
cd /g/data/lf10/as7425/2020-11_mouseBrain/data/2021-10-26_mouseBrain-noSplice-CHEUI/

# pick a model II dataset
export model_ii="/g/data/lf10/as7425/2020-11_mouseBrain/data/2021-10-26_mouseBrain-noSplice-CHEUI/WT_E15_read_level.txt_site_level.txt"

# get 0.1% of transcripts in the file
[ -f ./transcript_list.txt ] && echo "0" || awk 'BEGIN {srand()} !/^$/ { if (rand() <= .000017) print $0}' ${model_ii} | cut -f1 | uniq | tr "." "\t" | cut -f1  > ./transcript_list.txt

# list the transcripts 
export tx_list="/g/data/lf10/as7425/2020-11_mouseBrain/data/2021-10-26_mouseBrain-noSplice-CHEUI/transcript_list.txt"

# subset the CHEUI output
grep -f ${tx_list} ${model_ii} > ~/R2Dtool/test/CHEUI_modelII_subset.txt

# use R2Dtool utilites to convert to bed-like format 
bash ../scripts/cheui_to_bed.sh ./CHEUI_modelII_subset.txt ./out_CHEUI_modelII.bed

# Test that all positons in the input bed file correspond to 'A' nucleotides 
cd ~/R2Dtool/test/
mkdir test_input
transcriptome="/g/data/lf10/as7425/genomes/mouse_genome/GRCm39/transcriptome/GRCm39_codingPlusNoncoding.fa"
cat ./out_CHEUI_modelII.bed | tail -n +2 | cut -f1-6 > ./test_input/input_trancsriptome_sites.bed
bedtools getfasta -fi ${transcriptome} -bed ./test_input/input_trancsriptome_sites.bed | tail -n +2 | awk 'NR%2==1' | sort | uniq -c | sort -nr 
rm  -rf test_input

# subset annotation
export annotation="/g/data/lf10/as7425/genomes/mouse_genome/GRCm39/Mus_musculus.GRCm39.104.chr.gtf"
export subset_annotation="./GRCm39_subset.gtf"
grep -f ${tx_list} ${annotation} > ${subset_annotation}

# subset annotation for GFF3 
export gff3="/g/data/lf10/as7425/genomes/mouse_genome/GRCm39/Mus_musculus.GRCm39.104.chr.gff3"
export tx_list="/g/data/lf10/as7425/2020-11_mouseBrain/data/2021-10-26_mouseBrain-noSplice-CHEUI/transcript_list.txt"
grep -f ${tx_list} ${gff3} > ./GRCm39_subset.gff3
