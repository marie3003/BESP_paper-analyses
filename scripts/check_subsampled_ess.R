#!/usr/bin/env Rscript
# Checks ESS (cutoff = 200) for every non-empty *.subsampled.log file under
# a BEAST_DIR, to confirm that subsampling to ~1000 samples didn't drop ESS
# below the threshold used in the original make_summary_tree ESS check.
#
# Usage:
#   Rscript scripts/check_subsampled_ess.R <beast_dir> <output_tsv>
# Example:
#   Rscript scripts/check_subsampled_ess.R results/run1/beast_inference results/run1/beast_inference/subsampled_ess_check.tsv

suppressMessages({
  library(coda)
  library(beastio)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript check_subsampled_ess.R <beast_dir> <output_tsv>")
}
beast_dir  <- args[1]
output_tsv <- args[2]

log_files <- list.files(beast_dir, pattern = "\\.subsampled\\.log$",
                         recursive = TRUE, full.names = TRUE)

cat(sprintf("Found %d subsampled log files.\n", length(log_files)))

results <- data.frame(
  file       = character(0),
  size_bytes = integer(0),
  status     = character(0),
  n_low_ess  = integer(0),
  min_ess    = numeric(0),
  low_params = character(0),
  stringsAsFactors = FALSE
)

for (f in log_files) {
  sz <- file.info(f)$size
  if (is.na(sz) || sz == 0) {
    results <- rbind(results, data.frame(
      file = f, size_bytes = 0, status = "empty (ESS failed upstream)",
      n_low_ess = NA, min_ess = NA, low_params = NA
    ))
    next
  }

  status <- tryCatch({
    mcmc <- readLog(f, burnin = 0)
    ess_vals <- effectiveSize(mcmc)
    low <- checkESS(mcmc, cutoff = 200, value = TRUE)
    list(
      status = if (length(low) > 0) "FAIL" else "PASS",
      n_low_ess = length(low),
      min_ess = min(ess_vals),
      low_params = if (length(low) > 0) paste(names(low), collapse = ";") else ""
    )
  }, error = function(e) {
    list(status = paste("ERROR:", conditionMessage(e)),
         n_low_ess = NA, min_ess = NA, low_params = NA)
  })

  results <- rbind(results, data.frame(
    file = f, size_bytes = sz, status = status$status,
    n_low_ess = status$n_low_ess, min_ess = status$min_ess,
    low_params = status$low_params
  ))
}

write.table(results, output_tsv, sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n--- Summary ---\n")
print(table(results$status))
cat(sprintf("\nFull results written to: %s\n", output_tsv))

failed <- results[results$status == "FAIL", , drop = FALSE]
if (nrow(failed) > 0) {
  cat("\nFiles with ESS < 200 after subsampling:\n")
  print(failed$file)
}
