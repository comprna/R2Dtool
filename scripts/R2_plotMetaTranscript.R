#!/usr/bin/env Rscript

## read positional arguments while ensuring that the correct number of inputs are provided 
args <- commandArgs(trailingOnly = TRUE)

## check for positional arguments 

# -c flag for selecting confidence interval method, options are loess (default), or binom confint (-c binom)
ci_method <- "loess"
if ("-c" %in% args) {
  index <- which(args == "-c")
  ci_method <- args[index + 1]
  args <- args[-c(index, index + 1)]  # Remove the flag and its value from the args
}

# -o, optional output of processed data table 
output_path <- NULL
if ("-o" %in% args) {
  index <- which(args == "-o")
  output_path <- args[index + 1]
  args <- args[-c(index, index + 1)]  # Remove the flag and its value from the args
}

# -l flag for adding metagene labels to the plot (default: NO)
add_labels <- FALSE
if ("-l" %in% args) {
  add_labels <- TRUE
  index <- which(args == "-l")
  args <- args[-index]  # Remove the flag from the args
}

# Check number of positional arguments
if (length(args) != 5) {
  stop("\nUsage: Rscript R2_plotMetaTranscript.R '/path/to/annotated.bed' '/path/to/output.png' '<probability field>' '<cutoff>' '<upper/lower>' [-c 'loess'/'binom'] [-o '/path/to/output.tsv']", call. = FALSE)
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
suppressMessages(suppressWarnings(library(binom)))

# define function to classify calls as significant or non-significant 
filter_calls <- function(input_file, col_name, cutoff, direction) {
  
  # read in annotated modification positions while supressing most output except errors 
  calls <- tryCatch({
    read_tsv(file = input_file, col_names = T, guess_max = 99999999)
  }, error = function(e) {
    message("Error reading the file: ", e$message)
    stop(e)
  })
  
  # Ensure the column exists in the dataset
  if (!col_name %in% names(calls)) {
    stop(paste("Column", col_name, "does not exist in the input file."))
  }
  
  # Filter based on direction
  if (direction == "upper") {
    calls <- calls %>% mutate(filter = ifelse(.data[[col_name]] > cutoff, "sig", "ns"))
  } else {
    calls <- calls %>% mutate(filter = ifelse(.data[[col_name]] < cutoff, "sig", "ns"))
  }
  
  return(calls)
}

# define function to calculate the proportion of significant sites in metatranscript bins based on the "transcript_metacoordinate" column from R2D_annotate 
# Modify the compute_ratio function to ensure intervals are generated correctly
compute_ratio <- function(calls, ci_method) {
  breaks <- seq(0, 3, 0.025)
  out_ratio <- calls %>%
    mutate(interval = cut(transcript_metacoordinate, breaks, include.lowest = TRUE, right = TRUE, labels = FALSE)) %>%
    group_by(interval, filter) %>%
    summarise(n = n(), .groups = 'drop') %>%
    pivot_wider(names_from = "filter", values_from = "n", values_fill = list(n = 0)) %>%
    mutate(ratio = sig / (sig + ns + 1e-9))

  if (ci_method == "binom") {
    out_ratio <- out_ratio %>% filter(sig > 0 & ns > 0)
    conf_int <- binom::binom.confint(out_ratio$sig, out_ratio$sig + out_ratio$ns, methods = "wilson")
    out_ratio$lower <- conf_int$lower
    out_ratio$upper <- conf_int$upper
  }

  return(out_ratio)
}

# Modify the plot_ratio function to check for valid data
plot_ratio <- function(out_ratio, ci_method, add_labels) {
  if (nrow(out_ratio) == 0 || is.infinite(max(out_ratio$interval))) {
    stop("out_ratio is empty or contains invalid interval data.")
  }

  x_positions <- seq(0, max(out_ratio$interval), length.out = 7)
  y_position <- if ("upper" %in% names(out_ratio)) {
    0.87 * max(out_ratio$upper, na.rm = TRUE)
  } else {
    1.03 * max(out_ratio$ratio, na.rm = TRUE)
  }
  p <- ggplot(out_ratio, aes(x = interval, y = ratio)) +
    geom_point(alpha = 0.5, color = "red") +
    theme_minimal() +
    theme(text = element_text(size = 14)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle("Proportion of significant sites across metatranscript bins") +
    xlab("Relative metatranscriptomic location") + ylab("Proportion of significant sites") +
    geom_vline(xintercept = c(80,40), col = "black") + 
    coord_cartesian(xlim = c(0, max(out_ratio$interval)), ylim = c(0, y_position))

  if (add_labels) {
    p <- p +
      geom_text(label = "5' UTR", x = x_positions[2], y = y_position, vjust = -0.5, size = 10) +
      geom_text(label = "CDS", x = x_positions[4], y = y_position, vjust = -0.5, size = 10) +
      geom_text(label = "3'UTR", x = x_positions[6], y = y_position, vjust = -0.5, size = 10)
  }

  if (ci_method == "binom") {
    p <- p + geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
      geom_smooth(span = 0.2, color = "red", se = FALSE)
  } else {
    p <- p + geom_smooth(span = 0.2, color = "red", se = TRUE)
  }

  return(p)
}


# define the significant sites 
calls <- filter_calls(input_file, col_name, cutoff, direction)

# compute the ratio of significant sites across each bin
out_ratio <- compute_ratio(calls, ci_method)

out_ratio <- out_ratio %>% filter(!is.na(interval))

# for diagnostic; print the output ratio 
# print(out_ratio)

# write out_ratio to path set by -o flag 
if (!is.null(output_path)) {
  print(paste("Writing out_ratio to", output_path))
  write_tsv(out_ratio, output_path)
}

# plot the ratio of significant sites at each metatranscript bin 
p <- plot_ratio(out_ratio, ci_method, add_labels)

# save the plot
ggsave(output_file, p, scale = 4, width = 850, height = 750, units = c("px"))
