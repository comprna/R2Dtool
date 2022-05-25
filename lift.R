#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)!=4) {
  stop("\nUsage: Rscript liftoverCustomTranscriptome.R /path/to/bar.gtf /path/to/foo.bed /path/to/output.bed cheui-model", call.=FALSE)
}

library(GenomicFeatures)
library(rtracklayer)
library(tidyverse)

# interactive args for testing (B iMac)
# args <- c("/Users/asethi/localGadiData/2022-01-31_make-annotate-2/annotate_out.bed.temp.gtf", "/Users/asethi/localGadiData/2022-01-31_make-annotate-2/annotate_out.bed.tempbed", "/Users/asethi/localGadiData/2022-01-31_make-annotate-2/annotate_out.bed", "pval")

# interactive args for testing (JCSMR; Rat)
# args <- c("/Users/AJlocal/localGadiData/2022-05-05_test-annotate-rat/2022-05-05_annotated-AR3-100ksubsample.txt.temp.gtf", "/Users/AJlocal/localGadiData/2022-05-05_test-annotate-rat/2022-05-05_annotated-AR3-100ksubsample.txt.tempbed", "/Users/AJlocal/localGadiData/2022-05-05_test-annotate-rat/2022-05-05_annotated-AR3-100ksubsample.txt", "II")

# interactive args for testing (JCSMR; Mouse AR1 complete)
#
args <- c("/Users/AJlocal/localGadiData/2022-05-08_cheui-annotate-AR1-AR4/mouse_AR1_annotated.txt.temp.gtf",
          "/Users/AJlocal/localGadiData/2022-05-08_cheui-annotate-AR1-AR4/mouse_AR1_annotated.txt.tempbed",
          "/Users/AJlocal/localGadiData/2022-05-08_cheui-annotate-AR1-AR4/mouse_AR1_annotated.txt", "II")

# interactive args for testing (JCSMR; Rat AR4 complete)
# args <- c("/Users/AJlocal/localGadiData/2022-05-08_cheui-annotate-AR1-AR4/rat_AR4_annotated.txt.temp.gtf", "/Users/AJlocal/localGadiData/2022-05-08_cheui-annotate-AR1-AR4/rat_AR4_annotated.txt.tempbed", "/Users/AJlocal/localGadiData/2022-05-08_cheui-annotate-AR1-AR4/rat_AR4_annotated.txt", "II")

################################################################################
################################################################################
################################################################################

# start

tstart <- print(paste("start time is", as.character(Sys.time())))

##################################################

# process the annotation and make a reference table with transcript structure information

# read in reference transcripts
gtf <- makeTxDbFromGFF(file=args[1], format = "gtf")

# make an exon database from the reference transcripts
exons <- exonsBy(gtf, "tx", use.names=TRUE)

# convert the transcripts to a tibble
exons_tib <- as_tibble(as(exons, "data.frame"))

# fetch the length of each transcript segment from the gtf
txlen <- transcriptLengths(gtf, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE) %>%
  as_tibble() %>%
  mutate(diff = tx_len - cds_len - utr5_len - utr3_len) %>%
  dplyr::rename(transcript_id = tx_name) %>%
  dplyr::select(-tx_id, -nexon, -gene_id, -diff)

# the last command doesn't store biotype, so we read in the gtf again using another package
tx_biotype <- rtracklayer::import(args[1]) %>%
  as_tibble() %>%
  dplyr::select(transcript_id, transcript_biotype, gene_name, gene_id) %>%
  na.omit() %>%
  dplyr::distinct()

# merge the biotypes with the transcript segment lengths
merged_metadata <- inner_join(tx_biotype, txlen, by = "transcript_id")

##################################################

# lift-over the bed-like CHEUI output

# import bed file of transcriptome alignments
mappedLocus <- read_tsv(file = args[2], col_names = F) %>%
  dplyr::rename(transcript_id = 1, start = 2, end = 3, name = 4, score = 5, strand = 6) %>%
  dplyr::select(-strand) # %>%
  # separate(transcript, into=c("transcript_id", "transcript_version"), sep = "([.])", extra = "merge") %>%
  # dplyr::select(-transcript_version)

# make lookup table for strand
print("preparing strand lookup table")
strand_lookup <- exons_tib %>%
  dplyr::rename(transcript_id = group_name) %>%
  dplyr::select(transcript_id, strand) %>% dplyr::distinct()
  # mutate(transcript_id = gsub("\\..*","", transcript_id)) %>%
  # dplyr::distinct()

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

# based on the CHEUI runmode, separate the data column into it's actual values
if (args[4] == "pval"){
  output <- output %>%
    separate(data, into=c("tx_coord", "motif", "coverage", "stoich", "prob", "p.raw"), sep = ";") %>%
    type_convert(col_types = "fiiffficiddd")
  output$p.adj <- p.adjust(output$p.raw, method = "fdr")
} else if(args[4] == "II") {
  output <- output %>%
    separate(data, into=c("tx_coord", "motif", "coverage", "stoich", "prob"), sep = ";") %>%
    type_convert(col_types = "fiiffficidd")
  # output$p.adj <- p.adjust(output$p.raw, method = "fdr")
} else if(args[4] == "diff") {
  # to be written
  #output <- output %>%
  #  separate(data, into=c("tx_coord", "motif", "coverage", "stoich", "prob"), sep = ";") %>%
  #  type_convert(col_types = "fiiffficidd")
  #output$p.adj <- p.adjust(output$p.raw, method = "fdr")
} else{
  print("Invalid script run mode")
  stop()
}

##################################################

# attach the transcript lengths and biotype to output using the transcriptome coordinate
merge_out <- inner_join(output, merged_metadata %>% dplyr::rename(transcript = transcript_id), by = "transcript")

meta <- merge_out %>%
  mutate(cds_start = utr5_len,
         cds_end = utr5_len + cds_len,
         tx_end = cds_end + utr3_len) %>%
  mutate(rel_pos = ifelse(tx_coord < cds_start, # if the site is in the 5' utr
         ((tx_coord)/(cds_start)), # the relative position is simply the position of the site in the UTR
          ifelse(tx_coord < cds_end, # else, if the site is in the CDS
          (1 + (tx_coord - utr5_len)/(cds_len)), # the relative position is betwee 1 and 2 and represents the fraction of the cds the site is in
          (2 + (tx_coord - utr5_len - cds_len) / utr3_len))),  # no final condition, the site must be in the 3' utr, similar as before but the rel_pos is between 2 and 3
         abs_cds_start = tx_coord - cds_start, # absolute pos relative to CDS start
         abs_cds_end = tx_coord - cds_end) # absolute pos relative to CDS end

##################################################

# calculate distance from splice sites

# make a tibble of splice junction coordinates from the exon data frame
tx_junctions <- exons_tib %>%
  mutate(transcript_id = group_name) %>%
  group_by(transcript_id) %>%
  filter(row_number() < n()) %>% # drop the last row of each group, (including single exon transcripts) which cannot represent a bona fide splice junction
  arrange(exon_rank, .by_group = TRUE) %>%
  dplyr::select(transcript_id, width, exon_rank) %>%
  mutate(junc_coord = cumsum(width)) %>%
  dplyr::rename(transcript = transcript_id) %>%
  dplyr::rename(tx_coord = junc_coord) %>%
  dplyr::select(transcript, tx_coord)

# for each tested site, calculate the closest upstream and downstream junctions (where present) in tible tx_junctions
junc_dist <- left_join(meta %>% dplyr::select(transcript, tx_coord), tx_junctions %>% dplyr::select(transcript, tx_coord), by = "transcript", suffix = c("", ".y")) %>%
  group_by(transcript) %>%
  mutate(up_junc_dist = tx_coord - tx_coord.y, down_junc_dist = tx_coord.y - tx_coord) %>%
  pivot_longer(., cols = c(up_junc_dist, down_junc_dist), names_to = "type", values_to = "dist") %>%
  group_by(transcript, tx_coord, type) %>%
  filter(dist > 0) %>%
  dplyr::slice(which.min(dist)) %>%
  ungroup() %>%
  dplyr::select(-tx_coord.y) %>%
  pivot_wider(names_from = type, values_from = dist)

# attach this data to our previous output
meta_dist <- left_join(meta, junc_dist %>% dplyr::select(transcript, tx_coord, up_junc_dist, down_junc_dist), by = c("transcript", "tx_coord"))

##################################################

# write the output
print("writing final output")
write_tsv(meta_dist, args[3], col_names = T, append = FALSE)

# if writing out from interactive use
# meta_dist %>% dplyr::rename("#chr" = chr) %>% write_tsv(meta_dist, args[3], col_names = T, append = FALSE)

tend <- print(paste("end time is", as.character(Sys.time())))

tstart
tend
# quit()
