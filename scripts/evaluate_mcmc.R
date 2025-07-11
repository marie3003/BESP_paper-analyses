library(coda)
library(beastio)

mcmc_successful <- function(log_path, burnin_frac = 0.1, cutoff = 200) {
  mcmc <- readLog(log_path, burnin = burnin_frac)
  low_ess <- checkESS(mcmc, cutoff = cutoff, value = TRUE)
  if (length(low_ess) == 0) {
    trees_path <- sub("\\.log$", ".trees", log_path)
    return(data.frame(trees_path = trees_path, burnin = floor(nrow(mcmc) * burnin_frac)))
  } else {
    return(NULL)
  }
}

args <- commandArgs(trailingOnly = FALSE)
script.path <- normalizePath(sub("--file=", "", args[grep("--file=", args)]))
script.dir <- dirname(script.path)
parent_folder <- file.path(script_dir, "../results/pop_size_simulations/simulation_results")
output_file <- file.path(script_dir, "successful_mcmc_runs.csv")

# Find all .log files recursively
log_files <- list.files(parent_folder, pattern = "\\.log$", full.names = TRUE, recursive = TRUE)

# Evaluate each log file
results <- do.call(rbind, lapply(log_files, mcmc_successful))

# Write result to CSV if any successful runs
if (!is.null(results)) {
  write.table(results, file = output_file, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = ",")
}

# Filter for successful logs
#successful_logs <- log_files[sapply(log_files, mcmc_successful)]

# Write the result into the output file
#writeLines(successful_logs, output_file)

