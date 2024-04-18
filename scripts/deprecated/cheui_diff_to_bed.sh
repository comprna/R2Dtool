#!/bin/bash

# written by AJ Sethi on 2022-05-26
# convert CHEUI model II output to 6+ column bed-like file with headers
# arg1: CHEUI model II output file
# arg2: Output file path

####################################################################################################
####################################################################################################
####################################################################################################

export cheui_model_ii_output="${1}"
export bedlike_cheui_model_ii_output="${2}"

printf "[cheui-to-bed] $(date) ..... converting ${cheui_model_ii_output##*/} to bed-like file\n"

# write header
printf "transcript\tstart\tend\tname\tscore\tstrand\tcov_1\tcov_2\tstoich_1\tstoich_2\tstoich_diff\tstat\tpval\tpadj\ttranscript\tposition\tkmer\n" > "${bedlike_cheui_model_ii_output}"

# rearrange input columns to make 6-col bed-like output

# commented code concatenates data into name columnn
# cat ${cheui_model_ii_output} | tail -n +2 | tr "_" "\t" | awk '{ print $1"\t"$2+4"\t"$2+5"\t"$2+4";"$3";"$4";"$5";"$6"\t.\t+"}' >> "${bedlike_cheui_model_ii_output}" || printf "[cheui-to-bed] $(date) ..... [error] conversion for ${cheui_model_ii_output} failed!\n"

# new code adds additional data to columns 7+
cat ${cheui_model_ii_output} | tail -n +2 | tr "_" "\t" | awk '{ print $2"\t"$3+4"\t"$3+5"\t.\t.\t+\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15}' >> "${bedlike_cheui_model_ii_output}" || printf "[cheui-to-bed] $(date) ..... [error] conversion for ${cheui_model_ii_output} failed!\n"

# check word count
count=$(cat "${bedlike_cheui_model_ii_output}" | wc -l | xargs)

# tell user output
printf "[cheui-to-bed] $(date) ..... succesfully converted ${cheui_model_ii_output} and wrote output to ${bedlike_cheui_model_ii_output##*/}\n"
printf "[cheui-to-bed] $(date) ..... [note] ${bedlike_cheui_model_ii_output} contains ${count} records\n"
