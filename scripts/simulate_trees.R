library(phylodyn)
library(bdskytools)

source("SimUtils.R")

# Global settings
nreplicates <- 3
samp_start  <- 0
samp_end    <- 100
nlimit      <- 10
nrsamples   <- 200
outputbase  <- "../results/pop_size_simulations/"

trajectories <- data.frame(row.names  =c("expgrowth","uniform","bottleneck"),
                           names      =c("Exponential growth","Uniform","Bottleneck"),
                           maxdensity = c(0.4, 0.05, 0.05),
                           maxlineages = c(300, 100, 100))

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
set.seed(9)
outputpath <- paste0(outputbase,"linear_constant/")
for (i in 1:nrow(trajectories)) {
  
  trajname <- rownames(trajectories)[i]
  cat(paste0("## ",trajectories$names[i],"\n"))
  
  traj <- lapply(rep(trajname, nreplicates), get_trajectory)
  sims <- lapply(traj, simulate_genealogy, samp_type="preferential", nrsamples=nrsamples, samp_start=samp_start, samp_end=samp_end, nlimit=nlimit)
  sims <- save_simulation(sims, basename=trajname, path=paste0(outputpath,trajname,"/"), RData=TRUE, csv=TRUE, newick=TRUE, json=TRUE)
}