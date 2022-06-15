#!/bin/bash

# Written by AJ Sethi on 2022-05-14
# Pipelined script to run txannotate components to annotate and lift-over CHEUI model II output

export cheui_model_ii_output="${1}"
export annotation="${2}"
export lift_annotate_output="${3}"

##################################################

# get scriptpath
[ -n "${SCRIPTPATH+set}" ] || SCRIPTPATH="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )" || die "cannot get script path";
export SCRIPTPATH="${SCRIPTPATH}" # get the script path # taken from https://stackoverflow.com/questions/4774054/reliable-way-for-a-bash-script-to-get-the-full-path-to-itself

##################################################

# tell user output
printf "[txannotate-cheui-pipeline] $(date) ..... starting to convert ${cheui_model_ii_output}\n"

# convert CHEUI model II output to BED

bash ${SCRIPTPATH}/cheui_to_bed.sh "${cheui_model_ii_output}" "${lift_annotate_output}_temp.bed" || echo "failed to convert ${cheui_model_ii_output}"

# annotate bed-like transcriptomic sites

Rscript ${SCRIPTPATH}/annotate.R "${lift_annotate_output}_temp.bed"  "${annotation}" "${lift_annotate_output}_temp_annotated.bed"  || echo "failed to annotate ${cheui_model_ii_output}"

# liftover annotated transcriptomic sites

Rscript ${SCRIPTPATH}/lift.R "${lift_annotate_output}_temp_annotated.bed" "${annotation}" "${lift_annotate_output}" || echo "failed to liftover ${cheui_model_ii_output}"

# finish by deleting all temporary files made by the script
# rm ${lift_annotate_output}_temp*

# tell user the script is done
printf "[cheui-to-bed] $(date) ..... converting ${cheui_model_ii_output##*/} to bed-like file\n"

##################################################
