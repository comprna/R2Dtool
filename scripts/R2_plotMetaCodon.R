#!/usr/bin/env Rscript

suppressMessages({
  library(tidyverse)
  library(binom)
})

help_message <- function() {
  cat("\nUsage: Rscript script.R '/path/to/annotated.bed' '/path/to/output.png' '<probability field>' '<cutoff>' '<upper/lower>' [-c 'loess'/'binom'] [-o '/path/to/output.tsv'] [-l] (-s | -e)\n")
  cat("Options:\n")
  cat("  -c  Set confidence interval method; default is 'loess', alternative is 'binom'.\n")
  cat("  -o  Set output path for processed data table.\n")
  cat("  -l  Add metagene labels to the plot; default is no labels.\n")
  cat("  -s  Plot from the start codon using 'abs_cds_start'.\n")
  cat("  -e  Plot from the end codon using 'abs_cds_end'.\n")
  cat("  -h  Show this help message and exit.\n")
}

args <- commandArgs(trailingOnly = TRUE)
options <- list(ci_method = "loess", output_path = NULL, add_labels = FALSE, plot_from = NULL)

# help message 
if ("-h" %in% args || length(args) == 0) {
  help_message()
  quit(status = 0)
}

remaining_args <- list()
plot_from_specified <- FALSE

for (arg in args) {
  if (arg %in% c("-s", "-e")) {
    if (!is.null(options$plot_from)) {
      stop("Error: Both -s and -e flags cannot be used simultaneously.")
    }
    options$plot_from <- ifelse(arg == "-s", "abs_cds_start", "abs_cds_end")
    plot_from_specified <- TRUE
  } else if (arg %in% c("-c", "-o")) {
    next_index <- which(args == arg) + 1
    if (next_index > length(args)) {
      stop(paste("Error: The", arg, "flag requires a following argument."))
    }
    options[[sub("-", "", arg)]] <- args[next_index]
  } else if (arg == "-l") {
    options$add_labels <- TRUE
  } else {
    remaining_args <- c(remaining_args, arg)
  }
}

# check for codon flag 
if (!plot_from_specified) {
  help_message()
  stop("Error: Either -s or -e flag must be specified.")
}

# check remaining positional arguments
if (length(remaining_args) != 5) {
  help_message()
  stop("Error: Incorrect number of positional arguments.")
}

input_file <- as.character(remaining_args[1])
output_file <- as.character(remaining_args[2])
col_name <- as.character(remaining_args[3])
cutoff <- as.numeric(remaining_args[4])
direction <- as.character(remaining_args[5])

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
    mutate(filter = if_else(.data[[col]] > cutoff, "sig", "ns")) %>%
    mutate(filter = if (direction == "lower") if_else(filter == "sig", "ns", "sig") else filter)
  
  return(calls)
}

# calculate ratio of significant sites 
compute_ratio <- function(calls, interval) {
  # Filter calls to include only those within the specified range
  calls <- calls %>%
    filter(.data[[interval]] >= -100 & .data[[interval]] <= 100)

  # Calculate the count of 'sig' and 'ns' at each position
  out_ratio <- calls %>%
    group_by(.data[[interval]], filter) %>%
    summarise(n = n(), .groups = 'drop') %>%
    pivot_wider(names_from = "filter", values_from = "n", values_fill = list(n = 0)) %>%
    mutate(sig = if("sig" %in% names(.)) sig else 0,
           ns = if("ns" %in% names(.)) ns else 0,
           ratio = sig / (sig + ns + 1e-9)) %>%
    rename(interval = 1)

  # Calculate confidence intervals if the binom method is chosen
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

  p <- ggplot(out_ratio, aes(x = interval, y = ratio)) +
    geom_point(alpha = 0.5, color = "red") +
    geom_vline(xintercept = 0, color = "blue", linetype = "dashed") +
    geom_smooth(method = "loess", se = options$ci_method == "loess") +
    theme_minimal() +
    labs(title = paste("Proportion of significant sites around", title_prefix, "codon"),
         x = "Absolute metatranscriptomic location", y = "Proportion of significant sites")

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
if (!is.null(options$output_path)) {
  write_tsv(out_ratio, options$output_path)
}

p <- plot_ratio(out_ratio)

ggsave(output_file, p, scale = 4, width = 600, height = 400, units = "px")
