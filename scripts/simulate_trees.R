library(phylodyn)

# file needs to be run from BESP_paper-analysis project folder
#script_dir <- dirname(normalizePath(sys.frame(1)$ofile))
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  
  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  }
  
  return(getwd())
}

script_dir <- get_script_dir()
source(file.path(script_dir, "SimUtils.R"))

# Global settings
nreplicates <- 5
samp_start  <- 0
samp_end    <- 500
nlimit      <- 10
nrsamples   <- 250
outputbase <- file.path(script_dir, "../results/pop_size_simulations_new_param/")

trajectories <- data.frame(row.names  =c("expgrowth_fast", "expgrowth_slow", "uniform", "bottleneck"),
                           names      =c("Fast exponential growth", "Slow exponential growth","Uniform", "Bottleneck"),
                           # only used for plotting
                           maxdensity = c(0.4, 0.4, 0.05, 0.05),
                           maxlineages = c(300, 300, 100, 100))


# simulate trees with independent homogeneous sampling
set.seed(9)
outputpath <- paste0(outputbase,"independent_homochronous/")
for (i in 1:nrow(trajectories)) {
  
  trajname <- rownames(trajectories)[i]
  cat(paste0("## ",trajectories$names[i],"\n"))
  
  traj <- lapply(rep(trajname, nreplicates), get_trajectory)
  sims <- lapply(traj, simulate_genealogy, samp_type="independent", nrsamples=nrsamples, samp_start=samp_start, samp_end=samp_start, nlimit=nlimit)
  sims <- save_simulation(sims, basename=trajname, path=paste0(outputpath,trajname,"/"), RData=FALSE, csv=FALSE, newick=TRUE, json=TRUE)
}

# simulate trees with independent heterochroneous sampling
set.seed(9)
outputpath <- paste0(outputbase,"independent_heterochroneous/")
for (i in 1:nrow(trajectories)) {
  
  trajname <- rownames(trajectories)[i]
  cat(paste0("## ",trajectories$names[i],"\n"))
  
  traj <- lapply(rep(trajname, nreplicates), get_trajectory)
  sims <- lapply(traj, simulate_genealogy, samp_type="independent", nrsamples=nrsamples, samp_start=samp_start, samp_end=samp_end, nlimit=nlimit)
  sims <- save_simulation(sims, basename=trajname, path=paste0(outputpath,trajname,"/"), RData=FALSE, csv=FALSE, newick=TRUE, json=TRUE)
}

# simulate trees with preferential heterogeneous sampling, 1 epoch
set.seed(9)
outputpath <- paste0(outputbase,"linear_constant/")
for (i in 1:nrow(trajectories)) {
  
  trajname <- rownames(trajectories)[i]
  cat(paste0("## ",trajectories$names[i],"\n"))
  
  traj <- lapply(rep(trajname, nreplicates), get_trajectory)
  sims <- lapply(traj, simulate_genealogy, samp_type="preferential", nrsamples=nrsamples, samp_start=samp_start, samp_end=samp_end, nlimit=nlimit)
  sims <- save_simulation(sims, basename=trajname, path=paste0(outputpath,trajname,"/"), RData=FALSE, csv=FALSE, newick=TRUE, json=TRUE)
}

# simulate trees with preferential heterogeneous sampling, 24 epochs
set.seed(9)
outputpath <- paste0(outputbase,"linear_epoch_24/")
for (i in 1:nrow(trajectories)) {
  
  trajname <- rownames(trajectories)[i]
  cat(paste0("## ",trajectories$names[i],"\n"))
  
  traj <- lapply(rep(trajname, nreplicates), get_trajectory)
  sims <- lapply(traj, simulate_genealogy, samp_type="preferential", nrsamples=nrsamples, samp_start=samp_start, samp_end=samp_end, nlimit=nlimit, nrepochs=24)
  sims <- save_simulation(sims, basename=trajname, path=paste0(outputpath,trajname,"/"), RData=FALSE, csv=FALSE, newick=TRUE, json=TRUE)
}