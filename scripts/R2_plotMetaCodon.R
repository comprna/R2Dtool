#!/usr/bin/env Rscript

# Load required libraries quietly
suppressMessages({
  library(tidyverse)
  library(binom)
})

# Define help message
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

# Initialize parameters
args <- commandArgs(trailingOnly = TRUE)
options <- list(ci_method = "loess", output_path = NULL, add_labels = FALSE, plot_from = NULL)

if ("-h" %in% args || length(args) == 0) {
  help_message()
  quit(status = 0)
}

remaining_args <- list()
plot_from_specified <- FALSE

# Parse the command line arguments
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

# Validate required flags
if (!plot_from_specified) {
  help_message()
  stop("Error: Either -s or -e flag must be specified.")
}

# Validate number of remaining positional arguments
if (length(remaining_args) != 5) {
  help_message()
  stop("Error: Incorrect number of positional arguments.")
}

# Assign remaining arguments
input_file <- as.character(remaining_args[1])
output_file <- as.character(remaining_args[2])
col_name <- as.character(remaining_args[3])
cutoff <- as.numeric(remaining_args[4])
direction <- as.character(remaining_args[5])

# Check for file existence
if (!file.exists(input_file)) {
  stop("Input file does not exist")
}


# Define function to classify calls as significant or non-significant
filter_calls <- function(file, col, cutoff, direction) {
  calls <- tryCatch({
    read_tsv(file = file, col_names = TRUE, guess_max = 99999999)
  }, error = function(e) {
    message("Error reading the file: ", e$message)
    stop(e)
  })

  if (!col %in% names(calls)) {
    stop(paste("Column", col, "does not exist in the input file."))
  }

  calls <- calls %>% mutate(filter = ifelse(.data[[col]] > cutoff, "sig", "ns"), if (direction == "lower") filter = !filter)
  return(calls)
}

# Define function to calculate ratios and confidence intervals if needed
compute_ratio <- function(calls, interval) {
  # Filter calls to include only those within the specified range
  calls <- calls %>%
    filter(.data[[interval]] >= -100 & .data[[interval]] <= 100)

  # Calculate the count of 'sig' and 'ns' at each position
  out_ratio <- calls %>%
    group_by(.data[[interval]], filter) %>%
    summarise(n = n(), .groups = 'drop') %>%
    pivot_wider(names_from = "filter", values_from = "n", values_fill = list(n = 0)) %>%
    mutate(ratio = sig / (sig + ns + 1e-9)) %>%
    rename(interval = 1)

  # Calculate confidence intervals if the binom method is chosen
  if (options$ci_method == "binom") {
    conf_int <- binom.confint(out_ratio$sig, out_ratio$sig + out_ratio$ns, methods = "wilson")
    out_ratio <- mutate(out_ratio, lower = conf_int$lower, upper = conf_int$upper)
  }
  return(out_ratio)
}



# Define function to plot the ratio of significant sites
# Define function to plot the ratio of significant sites
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

# Processing the data
calls <- filter_calls(input_file, col_name, cutoff, direction)
out_ratio <- compute_ratio(calls, options$plot_from)

# Optional writing of results
if (!is.null(options$output_path)) {
  write_tsv(out_ratio, options$output_path)
}

# Plot and save
p <- plot_ratio(out_ratio)
ggsave(output_file, p, scale = 4, width = 850, height = 750, units = "px")
