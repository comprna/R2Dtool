#!/usr/bin/env Rscript

# import positional arguments 
args = commandArgs(trailingOnly=TRUE)

# test if there are 5 arguments, as required for plotMetaTranscript 
if (length(args)!=3) {
  stop("\nUsage: Rscript R2_plotMetaTranscript.R "/path/to/annotated.bed? "/path/to/output.png", "<probability field>", "<cutoff>", "<upper/lower>", call.=FALSE)
}

# arguments for testing 
args <- c("~/localGadiData/2022-09-21_R2Dtool/methylationCalls_annotated_lifted.bed",
          "~/localGadiData/2022-09-21_R2Dtool/out.png", 
          "probability", 
          "0.9999", 
          "upper")

# load tidyverse, quietly 
suppressMessages(suppressWarnings(library(tidyverse, warn.conflicts = F, quietly = T)))

# add target column as a variable 
colName = args[3]

# read in annotated transcriptomic positions
# use a function to define significant sites based on the user's 'upper' or 'lower' call 
calls <- read_tsv(file = args[1], col_names = T, guess_max = 99999999)

# select significant set based on user input 
if (args[5] == "upper") {
  calls <- calls %>% mutate(filter = ifelse(get({{colName}}) > as.numeric(args[4]), "sig", "ns"))
} else { 
  calls <- calls %>% mutate(filter = ifelse(get({{colName}}) < as.numeric(args[4]), "sig", "ns"))
}

# get splice junction data 
sj_data <- calls  %>% 
  mutate(up_junc_dist = -up_junc_dist) %>% 
  pivot_longer(cols = c(up_junc_dist, down_junc_dist), names_to = "type_dist", values_to = "junc_dist") %>% 
  dplyr::select(-type_dist) %>% 
  dplyr::filter(junc_dist < 355 & junc_dist > -355) %>% 
  group_by(junc_dist, filter) %>% 
  summarise(n = n()) %>% 
  pivot_wider(names_from = "filter", values_from = "n") %>% 
  mutate(ratio = sig / (sig + ns))

# plot 
p <- ggplot(sj_data, aes(x = junc_dist, y = ratio)) + 
  geom_point(alpha = 0.5, color = "blue") + 
  geom_smooth(span = 0.2) + # iterate over loess span 
  geom_vline(xintercept = 0, col = "black") + # iterate vertical lines to match the breaks 
  theme(text = element_text(size=12)) + 
  theme(plot.title = element_text(hjust=0.5)) + 
  ggtitle("Methylation around splice junctions in all transcripts") + 
  xlab("Absolute distance to splice junction (NT)") + ylab("Proportion of significant sites") +
  scale_y_log10() + 
  ylim(0,0.020) + 
  theme_minimal() 

ggsave(args[2], p, scale = 1,
       width = 115,
       height = 75,
       units = c("mm"))
