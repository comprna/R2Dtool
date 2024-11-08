#!/usr/bin/env Rscript

# read positional arguments while ensuring that the correct number of inputs are provided 
args <- commandArgs(trailingOnly = TRUE)

# help message 
if ("-h" %in% args || length(args) == 0) {
  stop("\nUsage: Rscript R2_plotMetaJunction.R '/path/to/annotated.bed' '/path/to/output.png' '<probability field>' '<cutoff>' '<upper/lower>' [-c 'loess'/'binom'] [-o '/path/to/output_table.tsv'] [-r]", call. = FALSE)
}

# check for positional argument for the confidence interval method
# options are "loess" (default) or binomial confidence interval (-c binom)
ci_method <- "loess"
if ("-c" %in% args) {
  index <- which(args == "-c")
  ci_method <- args[index + 1]
  args <- args[-c(index, index + 1)]  # Remove the flag and its value from the args
}

# check for optional -o flag for output path
output_path <- NULL
if ("-o" %in% args) {
  index <- which(args == "-o")
  output_path <- args[index + 1]
  args <- args[-c(index, index + 1)]  # Remove the flag and its value from the args
}

# check for optional -r flag to reverse x-axis
reverse_x <- FALSE
if ("-r" %in% args) {
  reverse_x <- TRUE
  index <- which(args == "-r")
  args <- args[-index]  # Remove the flag from the args
}

# check remaining positional arguments
if (length(args) != 5) {
  stop("\nUsage: Rscript R2_plotMetaJunction.R '/path/to/annotated.bed' '/path/to/output.png' '<probability field>' '<cutoff>' '<upper/lower>' [-c 'loess'/'binom'] [-o '/path/to/output_table.tsv'] [-r]", call. = FALSE)
}

print(paste("Using", ci_method, "method for confidence intervals."))

input_file <- args[1]
output_file <- args[2]
col_name <- args[3]
cutoff <- as.numeric(args[4])
direction <- args[5]

# ensure input exists
if (!file.exists(input_file)) {
  stop("Input file does not exist")
}

# load tidyverse quietly
suppressMessages(suppressWarnings(library(tidyverse)))

# define function to label significance
filter_calls <- function(input_file, col_name, cutoff, direction) {
  
  calls <- read_tsv(file = input_file, col_names = T, guess_max = 99999999)
  
  # ensure the filter field column is present in the R2Dtool annotate output 
  if (!col_name %in% names(calls)) {
    stop(paste("Column", col_name, "does not exist in the input file. R2Dtool annotate must be run in header mode (i.e. using -H flag) for plotMetaJunction to work"))
  }
  
  # filter for significant sites 
  if (direction == "upper") {
    calls <- calls %>% mutate(filter = ifelse(.data[[col_name]] >= cutoff, "sig", "ns"))
  } else {
    calls <- calls %>% mutate(filter = ifelse(.data[[col_name]] <= cutoff, "sig", "ns"))
  }
  
  return(calls)
}

# process splice junction positional data 
get_sj_data <- function(calls, ci_method) {
  
  sj_data <- calls  %>% 
    mutate(up_junc_dist = -up_junc_dist) %>% 
    pivot_longer(cols = c(up_junc_dist, down_junc_dist), names_to = "type_dist", values_to = "junc_dist") %>% 
    dplyr::select(-type_dist) %>% 
    dplyr::filter(junc_dist < 355 & junc_dist > -355) %>%  # maximum distance shown on plot
    group_by(junc_dist, filter) %>% 
    summarise(n = n(), .groups = 'drop') %>%
    pivot_wider(names_from = "filter", values_from = "n") %>%
    mutate(ratio = sig / (sig + ns + 1e-9))
  
  if (ci_method == "binom") {
    
    sj_data <- sj_data %>% filter(sig > 0 & ns > 0)
    
    # binomial confidence intervals 
    conf_int <- binom::binom.confint(sj_data$sig, sj_data$sig + sj_data$ns, methods = "wilson")
    sj_data$lower <- conf_int$lower
    sj_data$upper <- conf_int$upper
  }
  
  return(sj_data)
}

plot_sj_data <- function(sj_data, ci_method, reverse_x = FALSE) {
  p <- ggplot(sj_data, aes(x = if(reverse_x) -junc_dist else junc_dist, y = ratio)) + 
    geom_point(alpha = 0.5, color = "blue") + 
    geom_vline(xintercept = 0, col = "black") +
    theme(text = element_text(size = 12)) + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    ggtitle("Proportion of significant sites around splice junctions in all transcripts") + 
    xlab(if(reverse_x) "Distribution of m6A sites around splice junction" else "Absolute distance to splice junction (NT)") + 
    ylab("Proportion of significant sites") +
    theme_minimal()
  
  if (ci_method == "binom") {
    # binom confidence interval
    p <- p + geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2)
  } else {
    # loess smoothing
    p <- p + geom_smooth(span = 0.2)
  }
  
  return(p)
}

# filter calls 
calls <- filter_calls(input_file, col_name, cutoff, direction)

# calculate sj data 
sj_data <- get_sj_data(calls, ci_method)

# write sj data to output, if requested 
if (!is.null(output_path)) {
  print(paste("Writing sj_data to", output_path))
  write_tsv(sj_data, output_path)
}

# generate plot 
p <- plot_sj_data(sj_data, ci_method, reverse_x)

# save plot 
ggsave(output_file, p, scale = 4, width = 600, height = 400, units = c("px"))