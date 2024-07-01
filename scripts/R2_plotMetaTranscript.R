#!/usr/bin/env Rscript

# read positional arguments 
args <- commandArgs(trailingOnly = TRUE)

# help message 
if ("-h" %in% args || length(args) == 0) {
  stop("\nUsage: Rscript R2_plotMetaTranscript.R '/path/to/annotated.bed' '/path/to/output.png' '<probability field>' '<cutoff>' '<upper/lower>' [-c 'loess'/'binom'] [-o '/path/to/output.tsv']", call. = FALSE)
}

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
  args <- args[-index]  
}

# check remaining positional arguments
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
  
  # ensure the filter field column is present in the R2Dtool annotate output 
  if (!col_name %in% names(calls)) {
    stop(paste("Column", col_name, "does not exist in the input file. R2Dtool annotate must be run in header mode (i.e. using -H flag) for plotMetaTranscript to work"))
  }
 
  # ensure transcript metacoordinate and filter field are defined 
  calls <- calls %>% filter(!is.na(transcript_metacoordinate) & 
    is.finite(transcript_metacoordinate) & 
    !is.na(!!sym(col_name)) & 
    is.finite(!!sym(col_name)))
        
  # filter for significant sites
  if (direction == "upper") {
    calls <- calls %>% mutate(filter = ifelse(.data[[col_name]] >= cutoff, "sig", "ns"))
  } else {
    calls <- calls %>% mutate(filter = ifelse(.data[[col_name]] <= cutoff, "sig", "ns"))
  }
  

  return(calls)
}

# calculate the proportion of significant sites in metatranscript bins based on the "transcript_metacoordinate" column from R2D_annotate 
compute_ratio <- function(calls, ci_method) {
  breaks <- seq(0, 3, 0.025) # bin interval can be changed here, default, 120 bins of width 0.025
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

  out_ratio <- out_ratio %>% mutate(interval = interval / 40)
  return(out_ratio)
}

# make metatranscript plot 
plot_ratio <- function(out_ratio, ci_method, add_labels) {
  if (nrow(out_ratio) == 0 || is.infinite(max(out_ratio$interval))) {
    stop("out_ratio is empty or contains invalid interval data.")
  }

  x_positions <- seq(0, max(out_ratio$interval), length.out = 7)
  max_y <- if ("upper" %in% names(out_ratio)) {
    max(out_ratio$upper, na.rm = TRUE)
  } else {
    max(out_ratio$ratio, na.rm = TRUE)
  }
  
  y_position <- 1.1 * max_y
  y_limit <- 1.2 * max_y
  print(paste("Y limit is", y_limit))

  p <- ggplot(out_ratio, aes(x = interval, y = ratio)) +
    geom_point(alpha = 0.5, color = "red") +
    theme_minimal() +
    theme(text = element_text(size = 14)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle("Proportion of significant sites across metatranscript bins") +
    xlab("Relative metatranscriptomic location") + ylab("Proportion of significant sites") +
    geom_vline(xintercept = c(1,2), col = "black") + 
    expand_limits(y = y_limit)  

  if (add_labels) {
    p <- p +
      geom_text(label = "5' UTR", x = x_positions[2], y = y_position, vjust = -0.5, size = 5) +
      geom_text(label = "CDS", x = x_positions[4], y = y_position, vjust = -0.5, size = 5) +
      geom_text(label = "3'UTR", x = x_positions[6], y = y_position, vjust = -0.5, size = 5)
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

# write out_ratio to path set by -o flag 
if (!is.null(output_path)) {
  print(paste("Writing out_ratio to", output_path))
  write_tsv(out_ratio, output_path)
}

# plot the ratio of significant sites at each metatranscript bin 
p <- plot_ratio(out_ratio, ci_method, add_labels)

# save the plot
ggsave(output_file, p, scale = 4, width = 600, height = 400, units = c("px"))
