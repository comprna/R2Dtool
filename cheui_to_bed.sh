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
printf "transcript\tstart\tend\tname\tscore\tstrand\n" > "${bedlike_cheui_model_ii_output}"

# rearrange input columns to make 6-col bed-like output
cat ${cheui_model_ii_output} | tail -n +2 | tr "_" "\t" | awk '{ print $1"\t"$2+3"\t"$2+4"\t"$2+3";"$3";"$4";"$5";"$6"\t.\t+"}' >> "${bedlike_cheui_model_ii_output}" || printf "[cheui-to-bed] $(date) ..... [error] conversion for ${cheui_model_ii_output} failed!\n"

# check word count
count=$(cat "${bedlike_cheui_model_ii_output}" | wc -l | xargs)

# tell user output
printf "[cheui-to-bed] $(date) ..... succesfully converted ${cheui_model_ii_output} and wrote output to ${bedlike_cheui_model_ii_output##*/}\n"
printf "[cheui-to-bed] $(date) ..... [note] ${bedlike_cheui_model_ii_output} contains ${count} records\n"
