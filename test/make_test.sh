# navigate to cheui data directory
cd /g/data/lf10/as7425/2020-11_mouseBrain/data/2021-10-26_mouseBrain-noSplice-CHEUI

# pick a model II dataset
export model_ii="/g/data/lf10/as7425/2020-11_mouseBrain/data/2021-10-26_mouseBrain-noSplice-CHEUI/WT_E15_read_level.txt_site_level.txt"

# get 0.1% of transcripts in the file
[ -f ./transcript_list.txt ] && echo "0" || awk 'BEGIN {srand()} !/^$/ { if (rand() <= .000017) print $0}' ${model_ii} | cut -f1 | uniq | tr "." "\t" | cut -f1  > ./transcript_list.txt

export tx_list="/g/data/lf10/as7425/2020-11_mouseBrain/data/2021-10-26_mouseBrain-noSplice-CHEUI/transcript_list.txt"

# subset the CHEUI output
grep -f ${tx_list} ${model_ii} > ./WT_E15_CHEUI_model_ii_subset.txt

# subset annotation
export annotation="/g/data/lf10/as7425/genomes/mouse_genome/GRCm39/Mus_musculus.GRCm39.104.chr.gtf"
export subset_annotation="./GRCm39_subset.gtf"

grep -f ${tx_list} ${annotation} > ${subset_annotation}
