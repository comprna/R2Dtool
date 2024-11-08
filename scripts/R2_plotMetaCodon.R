#!/usr/bin/env Rscript

suppressMessages({
  library(tidyverse)
  library(binom)
})

cat("Arguments received:\n")
cat(paste(commandArgs(trailingOnly = TRUE), collapse = " "), "\n")

help_message <- function() {
  cat("\nUsage: Rscript script.R '/path/to/annotated.bed' '/path/to/output.png' '<probability field>' '<cutoff>' '<upper/lower>' [-c 'loess'/'binom'] [-o '/path/to/output.tsv'] [-l] (-s | -e) [-R]\n")
  cat("Options:\n")
  cat("  -c  Set confidence interval method; default is 'loess', alternative is 'binom'.\n")
  cat("  -o  Set output path for processed data table.\n")
  cat("  -l  Add metagene labels to the plot; default is no labels.\n")
  cat("  -s  Plot from the start codon using 'abs_cds_start'.\n")
  cat("  -e  Plot from the end codon using 'abs_cds_end'.\n")
  cat("  -R  Reverse x-axis sign and show distribution of m6A sites.\n")
  cat("  -h  Show this help message and exit.\n")
}

args <- commandArgs(trailingOnly = TRUE)
required_args <- args[1:5]
optional_args <- args[6:length(args)]

if (length(required_args) != 5) {
  cat("Error: Incorrect number of required arguments.\n")
  cat("Received arguments:", paste(args, collapse = " "), "\n")
  help_message()
  quit(status = 1)
}

input_file <- required_args[1]
output_file <- required_args[2]
col_name <- required_args[3]
cutoff <- as.numeric(required_args[4])
direction <- required_args[5]

options <- list(ci_method = "loess", output_path = NULL, add_labels = FALSE, plot_from = NULL, reverse_x = FALSE)

# parse optional arguments
i <- 1
while (i <= length(optional_args)) {
  flag <- optional_args[i]
  
  if (flag == "-c" && i + 1 <= length(optional_args)) {
    options$ci_method <- optional_args[i + 1]
    i <- i + 2
  } else if (flag == "-o" && i + 1 <= length(optional_args)) {
    options$output_path <- optional_args[i + 1]
    i <- i + 2
  } else if (flag == "-l") {
    options$add_labels <- TRUE
    i <- i + 1
  } else if (flag == "-s") {
    if (is.null(options$plot_from)) {
      options$plot_from <- "abs_cds_start"
    } else {
      stop("Error: Both -s and -e flags cannot be used simultaneously.")
    }
    i <- i + 1
  } else if (flag == "-e") {
    if (is.null(options$plot_from)) {
      options$plot_from <- "abs_cds_end"
    } else {
      stop("Error: Both -s and -e flags cannot be used simultaneously.")
    }
    i <- i + 1
  } else if (flag == "-R") {
    options$reverse_x <- TRUE
    i <- i + 1
  } else {
    warning(paste("Unrecognized flag:", flag))
    i <- i + 1
  }
}

# check for start or end flag 
if (is.null(options$plot_from)) {
  cat("Error: Either -s or -e flag must be specified.\n")
  help_message()
  quit(status = 1)
}

if (!file.exists(input_file)) {
  stop("Input file does not exist")
}

# filter for significant sites 
filter_calls <- function(file, col, cutoff, direction) {
  calls <- tryCatch({
    read_tsv(file = file, col_names = TRUE, guess_max = 99999999)
  }, error = function(e) {
    message("Error reading the file: ", e$message)
    stop(e)
  })
  
  # ensure the filter field column is present in the R2Dtool annotate output 
  if (!col %in% names(calls)) {
    stop(paste("Column", col, "does not exist in the input file. R2Dtool annotate must be run in header mode (i.e. with -H flag) for plotMetaCodon to work"))
  }
  
  calls <- calls %>% 
    mutate(filter = if_else(.data[[col]] >= cutoff, "sig", "ns")) %>%
    mutate(filter = if (direction == "lower") if_else(filter == "sig", "ns", "sig") else filter)
  
  return(calls)
}

# calculate ratio of significant sites 
compute_ratio <- function(calls, interval) {
  calls <- calls %>%
    filter(.data[[interval]] >= -100 & .data[[interval]] <= 100)
  
  out_ratio <- calls %>%
    group_by(.data[[interval]], filter) %>%
    summarise(n = n(), .groups = 'drop') %>%
    pivot_wider(names_from = "filter", values_from = "n", values_fill = list(n = 0)) %>%
    mutate(sig = if("sig" %in% names(.)) sig else 0,
           ns = if("ns" %in% names(.)) ns else 0,
           ratio = sig / (sig + ns + 1e-9)) %>%
    rename(interval = 1)
  
  if (options$ci_method == "binom") {
    conf_int <- binom.confint(out_ratio$sig, out_ratio$sig + out_ratio$ns, methods = "wilson")
    out_ratio <- mutate(out_ratio, lower = conf_int$lower, upper = conf_int$upper)
  }
  return(out_ratio)
}

# plot the ratio of significant sites 
plot_ratio <- function(out_ratio) {
  if (nrow(out_ratio) == 0 || is.infinite(max(out_ratio$interval))) {
    stop("out_ratio is empty or contains invalid interval data.")
  }
  
  title_prefix <- ifelse(options$plot_from == "abs_cds_start", "start", "stop")
  if (options$plot_from == "abs_cds_start") {
    labels <- c("5' UTR", "CDS")
  } else {
    labels <- c("CDS", "3' UTR")
  }
  
  x_values <- if(options$reverse_x) -out_ratio$interval else out_ratio$interval
  x_label <- if(options$reverse_x) "Distribution of m6A sites around codon" else "Distance from m6A site to codon"
  
  p <- ggplot(out_ratio, aes(x = x_values, y = ratio)) +
    geom_point(alpha = 0.5, color = "red") +
    geom_vline(xintercept = 0, color = "blue", linetype = "dashed") +
    geom_smooth(method = "loess", se = options$ci_method == "loess") +
    theme_minimal() +
    labs(title = paste("Proportion of significant sites around", title_prefix, "codon"),
         x = x_label, y = "Proportion of significant sites")
  
  if (options$add_labels) {
    max_ratio <- max(out_ratio$ratio, na.rm = TRUE)
    print(labels[1])
    print(labels[2])
    p <- p + geom_text(aes(label = labels[1]), x = -50, y = 0.9*max_ratio, vjust = -1,  size=12) +
      geom_text(aes(label = labels[2]), x = 50, y = 0.9*max_ratio, vjust = -1,  size=12)
  }
  
  if (options$ci_method == "binom") {
    p <- p + geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2)
  }
  
  return(p)
}

calls <- filter_calls(input_file, col_name, cutoff, direction)

out_ratio <- compute_ratio(calls, options$plot_from)

# optionally, write the data shown in the plot to a file 
print(options$output_path)

if (!is.null(options$output_path)) {
  write_tsv(out_ratio, options$output_path)
}

p <- plot_ratio(out_ratio)

ggsave(output_file, p, scale = 4, width = 600, height = 400, units = "px")