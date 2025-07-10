library(phylodyn)

# file needs to be run from BESP_paper-analysis project folder
#script_dir <- dirname(normalizePath(sys.frame(1)$ofile))
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
script_dir <- if (length(script_path) > 0) dirname(normalizePath(script_path)) else getwd()
setwd(script_dir)
source("SimUtils.R")

# Global settings
nreplicates <- 5
samp_start  <- 0
samp_end    <- 50
nlimit      <- 10
nrsamples   <- 250
outputbase  <- "../results/pop_size_simulations/"

trajectories <- data.frame(row.names  =c("expgrowth_fast", "expgrowth_slow", "uniform"),
                           names      =c("Fast exponential growth", "Slow exponential growth","Uniform"),
                           maxdensity = c(0.4, 0.4, 0.05),
                           maxlineages = c(300, 300, 100))

# simulate trees with independent homogeneous sampling
set.seed(9)
outputpath <- paste0(outputbase,"independent_homochronous/")
for (i in 1:nrow(trajectories)) {
  
  trajname <- rownames(trajectories)[i]
  cat(paste0("## ",trajectories$names[i],"\n"))
  
  traj <- lapply(rep(trajname, nreplicates), get_trajectory)
  sims <- lapply(traj, simulate_genealogy, samp_type="independent", nrsamples=nrsamples, samp_start=samp_start, samp_end=samp_start, nlimit=nlimit)
  sims <- save_simulation(sims, basename=trajname, path=paste0(outputpath,trajname,"/"), RData=TRUE, csv=TRUE, newick=TRUE, json=TRUE)
}

# simulate trees with preferential heterogeneous sampling
#set.seed(9)
#outputpath <- paste0(outputbase,"linear_constant/")
#for (i in 1:nrow(trajectories)) {
#  
#  trajname <- rownames(trajectories)[i]
#  cat(paste0("## ",trajectories$names[i],"\n"))
  
#  traj <- lapply(rep(trajname, nreplicates), get_trajectory)
#  sims <- lapply(traj, simulate_genealogy, samp_type="preferential", nrsamples=nrsamples, samp_start=samp_start, samp_end=samp_end, nlimit=nlimit)
#  sims <- save_simulation(sims, basename=trajname, path=paste0(outputpath,trajname,"/"), RData=TRUE, csv=TRUE, newick=TRUE, json=TRUE)
#}