#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)!=2) {
  stop("\nUsage: Rscript lift.R /path/to/foo.bed /path/to/output.bed", call.=FALSE)
}

library(tidyverse)

input_data <- read.table(args[1], header = TRUE)

input_data[c('transcript', 'position', "kmer")] <- str_split_fixed(input_data$ID, '_', 3)
input_data$position <- as.numeric(input_data$position)

#bed format
output_data <- input_data %>%
  dplyr::mutate(name = ".", score = ".", strand = "+", start = position, end = position+1) %>%
  dplyr::select(transcript, start, end, name, score, strand, coverage_1, coverage_2, stoichiometry_1, 
                stoichiometry_2, stoichiometry_diff, statistic, pval_U, position, kmer) 

write_tsv(output_data, file = args[2], col_names = TRUE)
