#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)!=3) {
  stop("\nUsage: Rscript lift.R /path/to/foo.bed /path/to/bar.gtf /path/to/output.bed", call.=FALSE)
}

library(GenomicFeatures)
library(rtracklayer)
library(tidyverse)

# test data 
# args <- c("~/Documents/txannotate/test/liftover_output.txt.tempbed", "~/Documents/txannotate/test/liftover_output.txt.temp.gtf", "~/Documents/txannotate/test/liftover_output.txt")

################################################################################
################################################################################
################################################################################

# lift-over the bed-like tsv of transcriptomic sites to their cognate genomic sites

# read in reference transcripts
gtf <- makeTxDbFromGFF(file=args[2], format = "gtf")

# make an exon database from the reference transcripts
exons <- exonsBy(gtf, "tx", use.names=TRUE)

# prepare the exons 
exons_tib <- as_tibble(as(exons, "data.frame"))

# import bed file of transcriptome alignments
mappedLocus <- read_tsv(file = args[1], col_names = F) %>%
  dplyr::rename(transcript = 1, start = 2, end = 3, name = 4, score = 5, strand = 6) %>%
  dplyr::select(-strand) %>%
  separate(transcript, into=c("transcript_id", "transcript_version"), sep = "([.])", extra = "merge") %>%
  dplyr::select(-transcript_version)

# make lookup table for strand
print("preparing strand lookup table")
strand_lookup <- exons_tib %>%
  dplyr::rename(transcript_id = group_name) %>%
  dplyr::select(transcript_id, strand) %>% dplyr::distinct() %>% 
  mutate(transcript_id = gsub("\\..*","", transcript_id)) %>%
  dplyr::distinct()

# attach the correct strand to the bed sites
print("reparing strand")
mappedLocus_fixedStrand <- inner_join(mappedLocus, strand_lookup, by = "transcript_id")

# write out the strand-repaired file as a temporary file
print("writing strand bedfile")
write_tsv(mappedLocus_fixedStrand, args[3], col_names = F)

# read in the bed sites with corrected strand
print("importing strand bedfile")
mappedLocus <- import.bed(args[3])

# map transcript coordinates to genome
print("mapping transcript coordinates to genome")
genomeLocus <- mapFromTranscripts(x=mappedLocus, transcripts=exons)

# bind score to output
# the score column contains cheui-specific output (e.g. stoich, prob, coverage, which we aggreagate into the score column and delimit using semicolons)
print("binding output")
mcols(genomeLocus)<-cbind(mcols(genomeLocus),DataFrame(mappedLocus[genomeLocus$xHits]))

# convert output to tibble
genome_coordinates = as_tibble(as(genomeLocus, "data.frame"))

# prepare the output by selecting bed-like coordinates from
print("filtering output")
output <- genome_coordinates %>% dplyr::select(seqnames, start, end, X.name, X.seqnames, strand) %>%
  unique() %>%
  dplyr::rename(chr = seqnames, data = X.name, transcript = X.seqnames) %>%
  mutate(score = ".") %>%
  dplyr::select(chr, start, end, transcript, score, strand, data)

# separate the output 
output <- output %>%
  separate(data, sep = ";")

# based on the CHEUI runmode, separate the data column into it's actual values
#if (args[4] == "pval"){
#  output <- output %>%
#     separate(data, into=c("tx_coord", "motif", "coverage", "stoich", "prob", "p.raw"), sep = ";") %>%
#     type_convert(col_types = "fiiffficiddd")
#   output$p.adj <- p.adjust(output$p.raw, method = "fdr")
# } else if(args[4] == "II") {
#   output <- output %>%
#     separate(data, into=c("tx_coord", "motif", "coverage", "stoich", "prob"), sep = ";") %>%
#     type_convert(col_types = "fiiffficidd")
#   # output$p.adj <- p.adjust(output$p.raw, method = "fdr")
# } else if(args[4] == "diff") {
#   # to be written
#   #output <- output %>%
#   #  separate(data, into=c("tx_coord", "motif", "coverage", "stoich", "prob"), sep = ";") %>%
#   #  type_convert(col_types = "fiiffficidd")
#   #output$p.adj <- p.adjust(output$p.raw, method = "fdr")
# } else{
#   print("Invalid script run mode")
#   stop()
# }

# write the output
print("writing final output")
write_tsv(output, args[3], col_names = T, append = FALSE)

# if writing out from interactive use
# meta_dist %>% dplyr::rename("#chr" = chr) %>% write_tsv(meta_dist, args[3], col_names = T, append = FALSE)

tend <- print(paste("end time is", as.character(Sys.time())))

tstart
tend
# quit()
