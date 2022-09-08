#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)!=3) {
  stop("\nUsage: Rscript annotate.R /path/to/foo.bed /path/to/bar.gtf /path/to/output.bed", call.=FALSE)
}

suppressMessages(suppressWarnings(library(GenomicFeatures, warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(library(rtracklayer, warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(library(tidyverse, warn.conflicts = F, quietly = T)))

################################################################################
################################################################################
################################################################################

# read in epitranscriptomic sites

mappedLocus <- read_tsv(file = args[1], col_names = T, guess_max = 999999999999) %>%
  dplyr::rename(transcript_id = 1) %>%
  mutate(transcript_id = gsub("\\..*","",transcript_id)) %>%
  dplyr::rename(tx_coord = 2)


##################################################

# fetch data from transcriptome annotation

# read in reference transcripts
gtf <- makeTxDbFromGFF(file=args[2], format = "gtf")

# make an exon database from the reference transcripts
exons <- exonsBy(gtf, "tx", use.names=TRUE)

# prepare the exons
exons_tib <- as_tibble(as(exons, "data.frame"))

# fetch the length of each transcript segment from the gtf
txlen <- transcriptLengths(gtf, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE) %>%
  as_tibble() %>%
  mutate(diff = tx_len - cds_len - utr5_len - utr3_len) %>%
  dplyr::rename(transcript_id = tx_name) %>%
  dplyr::select(-tx_id, -nexon, -gene_id, -diff)

# rename the gene_type field in case a gencode GTF is used (transcript_biotype) only works for Ensembl GTFs?
lookup <- c(transcript_biotype = "gene_type")

# the last command doesn't store biotype, so we read in the gtf again using another package
tx_biotype <- rtracklayer::import(args[2]) %>%
  as_tibble() %>%
  dplyr::rename(any_of(lookup)) %>%
  dplyr::select(transcript_id, transcript_biotype, gene_name, gene_id) %>%
  na.omit() %>%
  dplyr::distinct()

# merge the biotypes with the transcript segment lengths
merged_metadata <- inner_join(tx_biotype, txlen, by = "transcript_id") %>%
  mutate(transcript_id = gsub("\\..*","",transcript_id))

##################################################

# attach basal metadata to epitranscriptomic sites
merge_out <- inner_join(mappedLocus, merged_metadata, by = "transcript_id")


##################################################

# calculate metatranscipt coordinates

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
  dplyr::rename(tx_coord = junc_coord) %>%
  dplyr::select(transcript_id, tx_coord) %>%
    mutate(transcript_id = gsub("\\..*","",transcript_id))

# for each tested site, calculate the closest upstream and downstream junctions (where present) in tible tx_junctions
junc_dist <- left_join(meta %>% dplyr::select(transcript_id, tx_coord), tx_junctions %>% dplyr::select(transcript_id, tx_coord), by = "transcript_id", suffix = c("", ".y")) %>%
  group_by(transcript_id) %>%
  mutate(up_junc_dist = tx_coord - tx_coord.y, down_junc_dist = tx_coord.y - tx_coord) %>%
  pivot_longer(., cols = c(up_junc_dist, down_junc_dist), names_to = "type", values_to = "dist") %>%
  group_by(transcript_id, tx_coord, type) %>%
  filter(dist > 0) %>%
  dplyr::slice(which.min(dist)) %>%
  ungroup() %>%
  dplyr::select(-tx_coord.y) %>%
  pivot_wider(names_from = type, values_from = dist)

print("junc_dist")
print(head(junc_dist))

# attach this data to our previous output
meta_dist <- left_join(meta, junc_dist %>% dplyr::select(transcript_id, tx_coord, up_junc_dist, down_junc_dist), by = c("transcript_id", "tx_coord"))

##################################################

# rename output and print

final_output <- meta_dist %>%
  dplyr::rename(start = 2)

# write the output
print("writing final output")
write_tsv(final_output, args[3], col_names = T, append = FALSE)
