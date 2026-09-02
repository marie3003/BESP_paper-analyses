#!/usr/bin/env Rscript
# Debug helper: inspect what readLog() actually returns for one file.
# Usage: Rscript scripts/debug_readlog.R <path_to_subsampled_log>

suppressMessages({
  library(coda)
  library(beastio)
})

f <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(f)) stop("Usage: Rscript debug_readlog.R <path_to_subsampled_log>")

cat("File:", f, "\n")
cat("Lines in file (excluding # / Sample / State headers):\n")
system(paste0("grep -Evc '^(#|Sample|State)' '", f, "'"))

mcmc <- readLog(f, burnin = 0)
cat("\nclass(mcmc):", paste(class(mcmc), collapse=", "), "\n")
cat("is.mcmc.list:", coda::is.mcmc.list(mcmc), "\n")
cat("dim(mcmc):", paste(dim(mcmc), collapse=" x "), "\n")
cat("nvar:", tryCatch(coda::nvar(mcmc), error=function(e) "ERR"), "\n")
cat("niter:", tryCatch(coda::niter(mcmc), error=function(e) "ERR"), "\n")
cat("varnames:\n")
print(coda::varnames(mcmc))

cat("\nSummary of first few columns:\n")
mat <- as.matrix(mcmc)
print(dim(mat))
print(summary(mat[, 1:min(5, ncol(mat))]))

cat("\nAny all-NA columns?\n")
na_cols <- colnames(mat)[apply(mat, 2, function(x) all(is.na(x)))]
print(na_cols)

cat("\neffectiveSize result:\n")
print(tryCatch(coda::effectiveSize(mcmc), error = function(e) conditionMessage(e)))
