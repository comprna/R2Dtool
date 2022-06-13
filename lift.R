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
# args <- c("~/Documents/txannotate/test/out_CHEUI_modelII_annotated.bed", "~/Documents/txannotate/test/GRCm39_subset.gtf", "~/Documents/txannotate/test/out_CMII_annotated_lifted.bed")

################################################################################
################################################################################
################################################################################

# import bed file of transcriptome alignments

mappedLocus <- read_tsv(file = args[1], col_names = T) %>% 
  dplyr::rename(transcript_id = 1) %>% 
  mutate(transcript_id = gsub("\\..*","",transcript_id)) %>% 
  dplyr::rename(tx_coord = 2) %>% 
  dplyr::mutate(tx_coord_start = 2, tx_coord_end = 3)

# collect the column names of columns 7+ 
targetNames <- colnames(mappedLocus)[c(4,7:length(colnames(mappedLocus)))]

# merge columns c(4,7+)
mappedLocus <- unite(mappedLocus, metaname, c(4,7:length(colnames(mappedLocus))), sep = ">_>", remove = TRUE, na.rm = FALSE)

# deselect strand 
mappedLocus <- mappedLocus %>% select(-6)

##################################################

# fetch transcript structures from transcriptome annotation 

# read in reference transcripts
gtf <- makeTxDbFromGFF(file=args[2], format = "gtf")

# make an exon database from the reference transcripts
exons <- exonsBy(gtf, "tx", use.names=TRUE)

# prepare the exons 
exons_tib <- as_tibble(as(exons, "data.frame"))

# make lookup table for strand
print("preparing strand lookup table")
strand_lookup <- exons_tib %>%
  dplyr::rename(transcript_id = group_name) %>%
  dplyr::select(transcript_id, strand) %>% dplyr::distinct() %>% 
  mutate(transcript_id = gsub("\\..*","", transcript_id)) %>%
  dplyr::distinct()

##################################################

# attach the correct strand to the bed sites
print("reparing strand")
mappedLocus_fixedStrand <- inner_join(mappedLocus, strand_lookup, by = "transcript_id")

# write out the strand-repaired file as a temporary file
print("writing strand bedfile")
write_tsv(mappedLocus_fixedStrand, args[3], col_names = F)

##################################################

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
  dplyr::select(chr, start, end, transcript, score, strand, data) %>% 
  mutate(end = end + 1) %>% 
  rename("#chr" = chr)

# separate the output 
output <- output %>% separate(data, sep = ">_>", into = targetNames) 

##################################################

# fix the genome coordinate error 

# write the output
print("writing final output")
write_tsv(output, args[3], col_names = T, append = FALSE)


