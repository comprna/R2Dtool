#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)!=3) {
  stop("\nUsage: Rscript annotate.R /path/to/annotated.bed /path/to/output.png", "<probability field>", "<cutoff>", "<upper/lower>", call.=FALSE)
}

suppressMessages(suppressWarnings(library(tidyverse, warn.conflicts = F, quietly = T)))

# testing 
args <- c("~/localGadiData/2022-09-08_R2Dtool/methylationCalls_annotated_lifted.bed",
          "/localGadiData/2022-09-08_R2Dtool/out.png", 
          "probability", 
          "0.9999", 
          "upper")

################################################################################
################################################################################
################################################################################

# read in annotated transcriptomic positions
calls <- read_tsv(file = args[1], col_names = T, guess_max = 999999999999) %>% 
  mutate(filter = ifelse(probability > 0.9999, "sig", "ns"))

# normalized abundance plot
breaks <- seq(0,3,0.025) # iterate over break width 

Sys.time()
out_ratio <- calls %>% 
  mutate(interval = cut(rel_pos, breaks, include.lowest = TRUE, right = TRUE, labels = FALSE)) %>%
  group_by(interval, filter) %>% 
  summarise(n = n()) %>% 
  pivot_wider(names_from = "filter", values_from = "n") %>% 
  mutate(ratio = sig / (sig + ns))
  Sys.time()

# pre-render plot 
ggplot(out_ratio, aes(x = interval, y = ratio)) + 
  geom_point(alpha = 0.5, color = "red") + 
  geom_smooth(span = 0.2, color = "red") + # iterate over loess span 
  geom_vline(xintercept = c(80,40), col = "black") + # iterate vertical lines to match the breaks 
  theme_minimal() + 
  theme(text = element_text(size=14)) + 
  theme(plot.title = element_text(hjust=0.5)) + 
  ggtitle("Proportion of m6A/A in metatranscript bins") + 
  xlab("Relative metatranscriptomic location") + ylab("Proportion of significant sites") + 
  ylim(0,0.010) + 
  coord_cartesian(xlim=c(0,120))

ggsave("~/plot0.svg", plot = last_plot(), scale = 1,
       width = 115,
       height = 75,
       units = c("mm"))


# write the output
print("writing final output")
write_tsv(final_output, args[3], col_names = T, append = FALSE)


################################################################################
################################################################################


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
ggplot(sj_data, aes(x = junc_dist, y = ratio)) + 
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

ggsave("~/plot1.svg", plot = last_plot(), scale = 1,
       width = 115,
       height = 75,
       units = c("mm"))
