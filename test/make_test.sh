# sample CHEUI model II WT_E15 data
head -n 500 ${model_ii} > ${model_ii_test_data}

# make transcript list
cat * | cut -f1 | tr "." "\t" | cut -f1 > transcript_list.txt

# subset annotation
grep -f ${tx_list} ${annotation} > ${subset_annotation}
