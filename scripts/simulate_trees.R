library(phylodyn)
library(bdskytools)

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

# Optional CLI arguments: popmodel sampling_type
# e.g. Rscript simulate_trees.R bottleneckmid100 linearconstant
cli_args      <- commandArgs(trailingOnly = TRUE)
filter_pop    <- if (length(cli_args) >= 1) cli_args[1] else NULL
filter_samp   <- if (length(cli_args) >= 2) cli_args[2] else NULL

# Global settings
nreplicates <- 10
samp_start  <- 0
#samp_end    <- 100
nlimit      <- 10
nrsamples   <- 250
outputbase <- file.path(script_dir, "../results/run_bottleneck/simulated_data/")

# usual version
# trajectories <- data.frame(row.names         =c("expgrowthfast", "expgrowthslow", "uniform", "bottleneck"),
#                            names             =c("Fast exponential growth", "Slow exponential growth","Uniform", "Bottleneck"),
#                            # only used for plotting
#                            maxdensity        = c(0.1, 0.1, 0.05, 0.1),
#                            maxlineages       = c(300, 300, 100, 100),
#                            samp_end          = c(100, 100, 100, 30),
#                            maxtime_homo      = c(NA, NA, 500, 50),   # homochronous: cut to 100 years
#                            maxtime_hetero    = c(NA, NA, 500, 300))     # heterochronous: NULL=auto for expgrowth, fixed for others

# trying out new bottleneck scenarios
trajectories <- data.frame(row.names         =c("bottleneck20", "bottleneck50", "bottleneck100", "bottleneck200", "bottlenecklate", "bottlenecklatesampling", "bottleneckmid50", "bottleneckmid100"),
                          names             =c("Bottleneck (Depth 20)", "Bottleneck (Depth 50)", "Bottleneck (Depth 100)", "Bottleneck (Depth 200)", "Later Bottleneck", "Later Bottleneck (sampling through)", "Mid Bottleneck (Depth 50)", "Mid Bottleneck (Depth 100)"),
                          # only used for plotting
                          maxdensity        = c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1),
                          maxlineages       = c(100, 100, 100, 100, 100, 100, 100, 100),
                          samp_end          = c(30, 30, 30, 30, 30, 100, 30, 30),
                          maxtime_homo      = c(50, 50, 50, 50, 100, 100, 50, 50),
                          maxtime_hetero    = c(300, 300, 300, 300, 300, 300, 300, 300))


# --- Plotting helper functions ------------------------------------------------

getLTT <- function(tree, grid, plot=TRUE, col=pal.dark(cblue, 0.25)) {
  ltt <- ltt.plot.coords(tree)
  lttinterp <- approx(ltt[,1]*-1, ltt[,2], xout=grid, method='constant', rule=2)
  if (plot) {
    lines(ltt[,1]*-1, ltt[,2], col=col, type='S')
  }
  return(lttinterp$y)
}

getHPD <- function(sims, par, plot=TRUE, col=pal.dark(cblue), fill=pal.dark(cblue, 0.5), ...) {
  densities   <- lapply(sims, function(x) density(x[[par]], ...))
  densities_y <- sapply(densities, function(d) d$y)
  par_hpd     <- getMatrixHPD(t(densities_y))
  if (plot) {
    polygon(c(densities[[1]]$x, rev(densities[[1]]$x)), c(par_hpd[1,], rev(par_hpd[3,])), border=NA, col=fill)
    lines(densities[[1]]$x, par_hpd[2,], col=col)
  }
  return(par_hpd)
}

plotPointProcess <- function(sim, par, col=pal.dark(cblue)) {
  axis(1, lwd=0, lwd.ticks=1)
  abline(h=0)
  abline(v=sim[[par]], col=col)
  abline(v=max(sim[[par]]), lty=3)
  points(sim[[par]], rep(1, length(sim[[par]])), pch=20, col=col)
}

summariseSims <- function(sims, maxdensity=0.05, maxlineages=100, maxtime=NULL) {

  if (is.null(maxtime) || is.na(maxtime)) maxtime <- max(sims[[1]]$coal_times)

  layout(matrix(c(1,2,3,4,5,4,6,4), byrow=TRUE, ncol=2), heights=c(0.5,0.5,0.2,0.2))
  par(mar=c(5,4,2,3)+0.1)
  plot(sims[[1]]$trajgrid[,1], sims[[1]]$trajgrid[,2], type='l', lwd=2,
       xlim=c(0,maxtime), xlab="Time", ylab=expression("True N"[e]), las=1)
  abline(v=(sims[[1]]$samp_end), lty=3, lwd=1)

  if (!is.na(sims[[1]]$samp_intensity[1])) {
    par(new=TRUE)
    epoch_times <- seq(sims[[1]]$samp_start, sims[[1]]$samp_end, length.out=(sims[[1]]$nrepochs+1))
    plot(epoch_times, c(sims[[1]]$samp_intensity[1], sims[[1]]$samp_intensity),
         type='S', col=pal.dark(cred), lty=1, lwd=2,
         xlim=c(0,maxtime), ylim=c(0,max(sims[[1]]$samp_intensity)*1.25),
         bty='n', axes=FALSE, xlab="", ylab="")
    axis(4, las=1)
    legend('top', inset=c(0,-0.1), legend="Sampling intensity",
           col=pal.dark(cred), lty=1, lwd=2, bty='n', xpd=TRUE)
  }

  plot(1, type='n', xlim=c(0,maxtime), ylim=c(0,maxlineages),
       xlab="Time", ylab="LTT (95% HPD)", las=1)
  timegrid <- seq(0, maxtime, length.out=100)
  ltt_y    <- sapply(sims, function(x) getLTT(x$phylo, timegrid, col=pal.dark(cblue,0.1)))
  ltt_hpd  <- getMatrixHPD(t(ltt_y))
  lines(timegrid, ltt_hpd[1,], lty=2, col=pal.dark(cred))
  lines(timegrid, ltt_hpd[2,], lty=1, col=pal.dark(cred))
  lines(timegrid, ltt_hpd[3,], lty=2, col=pal.dark(cred))

  plot(1, type='n', xlim=c(0,maxtime), ylim=c(0,maxdensity),
       xlab="Time", ylab="Density (95% HPD)", las=1)
  samp_time_hpd <- getHPD(sims, par="samp_times", to=maxtime, bw=1,
                           col=pal.dark(cred), fill=pal.dark(cred,0.5))
  coal_time_hpd <- getHPD(sims, "coal_times", to=maxtime, bw=1,
                           col=pal.dark(cblue), fill=pal.dark(cblue,0.5))
  legend('topright', legend=c("Coalescent events","Sampling events"),
         border=pal.dark(c(cblue, cred)), fill=pal.dark(c(cblue,cred),0.5), bty='n')

  plot(ladderize(sims[[1]]$phylo), show.tip.label=FALSE,
       edge.width=0.5, direction="leftwards")
  mtext(side=1, "Example tree")

  plot(1, type='n', xlim=c(0,maxtime), ylim=c(0,1.1),
       xlab="Sampling events (example)", ylab="", las=1, bty='n', axes=FALSE)
  plotPointProcess(sims[[1]], "samp_times", col=pal.dark(cred,0.25))

  plot(1, type='n', xlim=c(0,maxtime), ylim=c(0,1.1),
       xlab="Coalescent events (example)", ylab="", las=1, bty='n', axes=FALSE)
  plotPointProcess(sims[[1]], "coal_times", col=pal.dark(cblue,0.25))
}


# --- Simulations --------------------------------------------------------------

# simulate trees with independent homogeneous sampling

if (is.null(filter_samp) || filter_samp == "independenthomochronous") {
  set.seed(9)
  outputpath <- file.path(outputbase, "independenthomochronous")
  traj_subset <- if (!is.null(filter_pop)) trajectories[rownames(trajectories) == filter_pop, , drop=FALSE] else trajectories
  for (i in 1:nrow(traj_subset)) {

    trajname <- rownames(traj_subset)[i]
    cat(paste0("## ",traj_subset$names[i],"\n"))

    traj <- lapply(rep(trajname, nreplicates), get_trajectory)
    sims <- lapply(traj, simulate_genealogy, samp_type="independent", nrsamples=nrsamples, samp_start=samp_start, samp_end=samp_start, nlimit=nlimit)
    sims <- save_simulation(sims, basename=trajname, path=paste0(file.path(outputpath, trajname), "/"), RData=TRUE, csv=FALSE, newick=TRUE, json=FALSE)

    pdf(file.path(outputpath, trajname, paste0(trajname, ".pdf")), width=12, height=10)
    summariseSims(sims, maxdensity=traj_subset$maxdensity[i], maxlineages=traj_subset$maxlineages[i], maxtime=traj_subset$maxtime_homo[i])
    title(main=traj_subset$names[i],
          sub=paste0("Nr. replicates: ",nreplicates, ", Nr. samples: ",nrsamples, ", Sampling period: [",samp_start,",",samp_start,"]"),
          outer=TRUE, line=-1.0)
    dev.off()
  }
}

# simulate trees with independent heterochroneous sampling
# set.seed(9)
# outputpath <- file.path(outputbase, "independent_heterochroneous")
# for (i in 1:nrow(trajectories)) {
# 
#   trajname <- rownames(trajectories)[i]
#   cat(paste0("## ",trajectories$names[i],"\n"))
# 
#   traj <- lapply(rep(trajname, nreplicates), get_trajectory)
#   sims <- lapply(traj, simulate_genealogy, samp_type="independent", nrsamples=nrsamples, samp_start=samp_start, samp_end=samp_end, nlimit=nlimit)
#   sims <- save_simulation(sims, basename=trajname, path=paste0(file.path(outputpath, trajname), "/"), RData=FALSE, csv=FALSE, newick=TRUE, json=FALSE)
# 
#   pdf(file.path(outputpath, trajname, paste0(trajname, ".pdf")), width=12, height=10)
#   summariseSims(sims, maxdensity=trajectories$maxdensity[i], maxlineages=trajectories$maxlineages[i])
#   title(main=trajectories$names[i],
#         sub=paste0("Nr. replicates: ",nreplicates, ", Nr. samples: ",nrsamples, ", Sampling period: [",samp_start,",",samp_end,"]"),
#         outer=TRUE, line=-1.0)
#   dev.off()
# }

#simulate trees with preferential heterogeneous sampling, 1 epoch
if (is.null(filter_samp) || filter_samp == "linearconstant") {
  set.seed(9)
  outputpath <- file.path(outputbase, "linearconstant")
  traj_subset <- if (!is.null(filter_pop)) trajectories[rownames(trajectories) == filter_pop, , drop=FALSE] else trajectories
  for (i in 1:nrow(traj_subset)) {

    trajname <- rownames(traj_subset)[i]
    cat(paste0("## ",traj_subset$names[i],"\n"))

    traj <- lapply(rep(trajname, nreplicates), get_trajectory)
    sims <- lapply(traj, simulate_genealogy, samp_type="preferential", nrsamples=nrsamples, samp_start=samp_start, samp_end=traj_subset$samp_end[i], nlimit=nlimit)
    sims <- save_simulation(sims, basename=trajname, path=paste0(file.path(outputpath, trajname), "/"), RData=TRUE, csv=FALSE, newick=TRUE, json=FALSE)

    pdf(file.path(outputpath, trajname, paste0(trajname, ".pdf")), width=12, height=10)
    summariseSims(sims, maxdensity=traj_subset$maxdensity[i], maxlineages=traj_subset$maxlineages[i], maxtime=traj_subset$maxtime_hetero[i])
    title(main=traj_subset$names[i],
          sub=paste0("Nr. replicates: ",nreplicates, ", Nr. samples: ~",nrsamples, ", Sampling period: [",samp_start,",",traj_subset$samp_end[i],"]"),
          outer=TRUE, line=-1.0)
    dev.off()
  }
}

# simulate trees with bottleneck trajectories (varying depth and width)
bottleneck_depths <- c(10, 20, 30, 40, 50)
bottleneck_widths <- c(3, 5, 10, 20)
# bottleneck starts at t=10, max Ne=1000, samp_end=max(stop)+10 to capture recovery
bn_samp_end    <- 50
bn_maxtime_homo <- 50
bn_maxtime_hetero <- 300

# set.seed(9)
# outputpath <- file.path(outputbase, "bottleneck_homochronous")
# for (depth in bottleneck_depths) {
#   for (width in bottleneck_widths) {
#     trajname <- paste0("bottleneck_d", depth, "_w", width)
#     cat(paste0("## Bottleneck depth=", depth, " width=", width, "\n"))
#     traj <- lapply(rep(trajname, nreplicates), get_trajectory)
#     sims <- lapply(traj, simulate_genealogy, samp_type="independent", nrsamples=nrsamples, samp_start=samp_start, samp_end=samp_start, nlimit=nlimit)
#     sims <- save_simulation(sims, basename=trajname, path=paste0(file.path(outputpath, trajname), "/"), RData=TRUE, csv=FALSE, newick=TRUE, json=FALSE)

#     pdf(file.path(outputpath, trajname, paste0(trajname, ".pdf")), width=12, height=10)
#     summariseSims(sims, maxdensity=0.1, maxlineages=100, maxtime=bn_maxtime_homo)
#     title(main=paste0("Bottleneck depth=", depth, " width=", width),
#           sub=paste0("Nr. replicates: ", nreplicates, ", Nr. samples: ", nrsamples, ", Sampling period: [", samp_start, ",", samp_start, "]"),
#           outer=TRUE, line=-1.0)
#     dev.off()
#   }
# }

# set.seed(9)
# outputpath <- file.path(outputbase, "bottleneck_heterochronous")
# for (depth in bottleneck_depths) {
#   for (width in bottleneck_widths) {
#     trajname <- paste0("bottleneck_d", depth, "_w", width)
#     cat(paste0("## Bottleneck depth=", depth, " width=", width, "\n"))
#     traj <- lapply(rep(trajname, nreplicates), get_trajectory)
#     sims <- lapply(traj, simulate_genealogy, samp_type="preferential", nrsamples=nrsamples, samp_start=samp_start, samp_end=bn_samp_end, nlimit=nlimit)
#     sims <- save_simulation(sims, basename=trajname, path=paste0(file.path(outputpath, trajname), "/"), RData=TRUE, csv=FALSE, newick=TRUE, json=FALSE)

#     pdf(file.path(outputpath, trajname, paste0(trajname, ".pdf")), width=12, height=10)
#     summariseSims(sims, maxdensity=0.1, maxlineages=100, maxtime=bn_maxtime_hetero)
#     title(main=paste0("Bottleneck depth=", depth, " width=", width),
#           sub=paste0("Nr. replicates: ", nreplicates, ", Nr. samples: ~", nrsamples, ", Sampling period: [", samp_start, ",", bn_samp_end, "]"),
#           outer=TRUE, line=-1.0)
#     dev.off()
#   }
# }

# simulate trees with preferential heterogeneous sampling, 24 epochs
# set.seed(9)
# outputpath <- file.path(outputbase, "linear_epoch_24")
# for (i in 1:nrow(trajectories)) {
# 
#   trajname <- rownames(trajectories)[i]
#   
#   
#   cat(paste0("## ",trajectories$names[i],"\n"))
# 
#   traj <- lapply(rep(trajname, nreplicates), get_trajectory)
#   sims <- lapply(traj, simulate_genealogy, samp_type="preferential", nrsamples=nrsamples, samp_start=samp_start, samp_end=samp_end, nlimit=nlimit, nrepochs=24)
#   sims <- save_simulation(sims, basename=trajname, path=paste0(file.path(outputpath, trajname), "/"), RData=FALSE, csv=FALSE, newick=TRUE, json=FALSE)
# 
#   pdf(file.path(outputpath, trajname, paste0(trajname, ".pdf")), width=12, height=10)
#   summariseSims(sims, maxdensity=trajectories$maxdensity[i], maxlineages=trajectories$maxlineages[i])
#   title(main=trajectories$names[i],
#         sub=paste0("Nr. replicates: ",nreplicates, ", Nr. samples: ~",nrsamples, ", Sampling period: [",samp_start,",",samp_end,"]"),
#         outer=TRUE, line=-1.0)
#   dev.off()
# }
