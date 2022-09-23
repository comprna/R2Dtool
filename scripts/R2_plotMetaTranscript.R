#!/usr/bin/env Rscript

# import positional arguments 
args = commandArgs(trailingOnly=TRUE)

# test if there are 5 arguments, as required for plotMetaTranscript 
if (length(args)!=3) {
  stop("\nUsage: Rscript R2_plotMetaTranscript.R "/path/to/annotated.bed? "/path/to/output.png", "<probability field>", "<cutoff>", "<upper/lower>", call.=FALSE)
}

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

# define breaks for plot 
breaks <- seq(0,3,0.025) # iterate over break width 

# cut the breaks 
out_ratio <- calls %>% 
  mutate(interval = cut(rel_pos, breaks, include.lowest = TRUE, right = TRUE, labels = FALSE)) %>%
  group_by(interval, filter) %>% 
  summarise(n = n()) %>% 
  pivot_wider(names_from = "filter", values_from = "n") %>% 
  mutate(ratio = sig / (sig + ns))

# plot 
p <- ggplot(out_ratio, aes(x = interval, y = ratio)) + 
  geom_point(alpha = 0.5, color = "red") + 
  geom_smooth(span = 0.2, color = "red") + # iterate over loess span 
  geom_vline(xintercept = c(80,40), col = "black") + # iterate vertical lines to match the breaks 
  theme_minimal() + 
  theme(text = element_text(size=14)) + 
  theme(plot.title = element_text(hjust=0.5)) + 
  ggtitle("Proportion of m6A/A in metatranscript bins") + 
  xlab("Relative metatranscriptomic location") + ylab("Proportion of significant sites") + 
  coord_cartesian(xlim=c(0,120)) + 
  ylim(0,0.0125)

ggsave(args[2], p, scale = 1,
       width = 115,
       height = 75,
       units = c("mm"))