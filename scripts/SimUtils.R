require(phylodyn)
require(ape)

bottleneck_traj_param <- function(t, start=0.5, stop=1, min=0.1, max=1) {
    result = rep(0,length(t))
    result[t <= start] <- max
    result[t > start & t < stop] <- min
    result[t >= stop] <- max
    return(result)
}

exp_traj_lower_limit = function(t, scale=1000, rate=1, lower_limit = 10)
{
  return(pmax(scale * exp(-t * rate), lower_limit))
}


#' Return trajectory function
#' 
#' @param type One of {uniform, expgrowth, boombust, logistic, cyclic, bottleneck}
get_trajectory <- function(type) {
  
    unif_upper       <- 1000
    unif_lower       <- 10
    exp_scale        <- 2000
    exp_lower_limit   <- 10
    exp_rate_fast    <- 0.02
    exp_rate_slow    <- 0.01
    bottleneck_start <- 10
    bottleneck_end   <- 13
    
  
    # if (grepl("^bottleneck_d[0-9]+_w[0-9]+$", type)) {
    #   params <- as.numeric(regmatches(type, gregexpr("[0-9]+", type))[[1]])
    #   bn_depth <- params[1]
    #   bn_width <- params[2]
    #   return(function(t) bottleneck_traj_param(t, min=bn_depth, max=unif_upper, start=bottleneck_start, stop=bottleneck_start+bn_width))
    # }

    switch(type,
          "uniform" = function(t) unif_traj(t, level = unif_upper),
          "expgrowthfast"  = function(t) exp_traj_lower_limit(t, scale=exp_scale, rate=exp_rate_fast, lower_limit=exp_lower_limit),
          "expgrowthslow"  = function(t) exp_traj_lower_limit(t, scale=exp_scale, rate=exp_rate_slow, lower_limit=exp_lower_limit),
          "bottleneck" = function(t) bottleneck_traj_param(t, min=unif_lower, max=unif_upper, start=bottleneck_start, stop=bottleneck_end),
          "bottleneck20" = function(t) bottleneck_traj_param(t, min=20, max=unif_upper, start=bottleneck_start, stop=bottleneck_end),
          "bottleneck50" = function(t) bottleneck_traj_param(t, min=50, max=unif_upper, start=bottleneck_start, stop=bottleneck_end),
          "bottleneck100" = function(t) bottleneck_traj_param(t, min=100, max=unif_upper, start=bottleneck_start, stop=bottleneck_end),
          "bottleneck200" = function(t) bottleneck_traj_param(t, min=200, max=unif_upper, start=bottleneck_start, stop=bottleneck_end),
          "bottlenecklate" = function(t) bottleneck_traj_param(t, min=unif_lower, max=unif_upper, start=50, stop=53),
          "bottlenecklatesampling" = function(t) bottleneck_traj_param(t, min=unif_lower, max=unif_upper, start=50, stop=53),
    )
}

traj_beta <- function(t, traj, beta) {
  pmax(traj(t), .Machine$double.eps)^beta
}


#'
#' constant power law
#' epoch times, evenly distributed at power law
#' independent
#' 
#' 
#' Constant linear sampling is identical to epoch sampling with just 1 epoch, so no need to implement both
simulate_genealogy <- function(traj, nrsamples=500, samp_start=0, samp_end=50, samp_type="preferential", nlimit=1, beta=1, nrepochs=1) {

    samp_end <- min(samp_end, tryCatch(uniroot(function(t) traj(t)-nlimit, interval = c(samp_start,samp_end))$root, error = function(x) NA), na.rm=TRUE)
    samp_end <- floor(samp_end)
    
    
    if (samp_type == "independent") {
        # Uniform probability of sampling times between samp_start and samp_end
        samp_intensity <- NA
        beta           <- NA
        samp_times     <- sort(c(samp_start, runif((nrsamples-1), min=samp_start, max=samp_end)))
        epoch_times    <- c(samp_start, samp_end)
    } else
    if (samp_type == "preferential") {
        # Distribute samples evenly in a number of epochs evenly spaced between samp_start and samp_end
        # Assumes a constant sampling intensity during each epoch
        # (Constant sampling intensity through time is equivalent to having only 1 epoch)
        epoch_times    <- seq(samp_start, samp_end, length.out=(nrepochs+1))
        epoch_samples  <- nrsamples/nrepochs
        samp_intensity <- c()
        samp_times     <- c()
        for (i in 1:nrepochs) {
            epoch_intensity <- epoch_samples/integrate(traj_beta, epoch_times[i], epoch_times[i+1], traj=traj, beta=beta, subdivisions = 1000, rel.tol = 1e-6, abs.tol = 1e-12)$value    
            samp_times      <- c(samp_times, pref_sample(traj, c=epoch_intensity, lim=c(epoch_times[i], epoch_times[i+1]), beta=beta))
            samp_intensity  <- c(samp_intensity, epoch_intensity)
        }
    } else {
        stop(paste0("Error! Unknown sampling type: ", samp_type))
    }

    # Simulate the genealogy under the population size trajectory and sampling times
    genealogy      <- coalsim(samp_times = samp_times, 
                              n_sampled  = rep(1,length(samp_times)), 
                              traj       = traj, 
                              lower_bound= 0.1, 
                              method     = "thin")   
    
    
    # Sanity check
    if (!all.equal(samp_times, genealogy$samp_times)) {
        stop("Error! Sampling times and simulated sampling times don't match...")
    }
    
    # Get a phylogenetic tree that corresponds to the genealogy
    phylogeny      <- generate_newick(genealogy)$newick
    
    return(list(nrsamples      = length(samp_times),
                samp_start     = samp_start,
                samp_end       = samp_end,
                nrepochs       = nrepochs,
                epoch_times    = epoch_times,
                samp_intensity = samp_intensity,
                beta           = beta,
                coal_times     = genealogy$coal_times,
                samp_times     = genealogy$samp_times,
                phylo          = phylogeny,
                traj           = traj))
}


#' Deprecated
simulate_genealogy_save <- function(traj, nrsamples=500, samp_end=50) {
  
  # Get the sampling intensity and sampling times (sampling intensity is a constant of proportionality)
  # Assumes a constant sampling intensity through time
  samp_intensity <- nrsamples/integrate(traj_beta, 0, samp_end, traj=traj, beta=1, subdivisions = 1000, rel.tol = 1e-6, abs.tol = 1e-12)$value    
  samp_times     <- pref_sample(traj, c=samp_intensity, lim=c(0,samp_end), beta=1)
  
  # Simulate the genealogy under the population size trajectory and sampling times
  genealogy      <- coalsim(samp_times = samp_times, 
                            n_sampled  = rep(1,length(samp_times)), 
                            traj       = traj, 
                            lower_bound= 10, 
                            method     = "thin")    
  
  # Sanity check
  if (!all.equal(samp_times, genealogy$samp_times)) {
    stop("Error! Sampling times and simulated sampling times don't match...")
  }
  
  # Get a phylogenetic tree that corresponds to the genealogy
  phylogeny      <- generate_newick(genealogy)$newick
  
  return(list(nrsamples      = length(samp_times), 
              samp_intensity = samp_intensity,
              coal_times     = genealogy$coal_times,
              samp_times     = genealogy$samp_times,
              phylo          = phylogeny,
              traj           = traj))
}


#' Sanitise a single simulation and grid the trajectory
prepare_simulation <- function(simresult, gridsize=1000, sanitise_times=TRUE, sanity_check=TRUE) {
    
    simresult$timeshift <- 0
    if (sanitise_times) {
      # Arbitrarily shift times so that most recent sample is taken at t=0.
      # Also removes the most recent time from the tip labels of the tree (so tip labels are generic)
      simresult$timeshift <- min(simresult$samp_times)
      
      # Shift coalescent and sampling times
      simresult$coal_times <- simresult$coal_times - simresult$timeshift
      simresult$samp_times <- simresult$samp_times - simresult$timeshift
      simresult$samp_end   <- simresult$samp_end   - simresult$timeshift
      simresult$samp_start <- simresult$samp_start - simresult$timeshift
      
      
      # Change tree tip labels
      for (i in 1:length(simresult$phylo$tip.label)) {
          simresult$phylo$tip.label[i] <- strsplit(simresult$phylo$tip.label[i],"_")[[1]][1]
      }
      
      # Sanity check (comparing sampling and coalescent times to tree times)
      # requires TreeSim (should get rid of dependency)
      if (sanity_check) {
        treetimes  <- TreeSim::getx(simresult$phylo, sersampling=TRUE)
        samp_times <- sort(treetimes[which(treetimes[,2]==0),1])
        coal_times <- sort(treetimes[which(treetimes[,2]==1),1])
        if (!all.equal(samp_times, simresult$samp_times) || !all.equal(coal_times, simresult$coal_times)) {
           stop("Error! Something wrong after shifting times to 0...")
        }
      }
    } 
    
    # Get trajectory on a timegrid
    timegrid <- seq(0,max(simresult$coal_times,simresult$samp_times), length.out=gridsize)
    trajgrid <- simresult$traj((simresult$timeshift+timegrid))
    simresult$trajgrid <- cbind(timegrid, trajgrid)
    
    return(simresult)
}


save_simulation_csv <- function(sim, path, basename) {
    #print(paste(path,basename))
    write.table(sim$coal_times,     file=paste0(path,basename,"_coaltimes.csv"),     quote=FALSE, row.names=FALSE, col.names="coalescent_times",   sep=",")
    write.table(sim$samp_times,     file=paste0(path,basename,"_samptimes.csv"),     quote=FALSE, row.names=FALSE, col.names="sampling_times",     sep=",")
    write.table(sim$epoch_times,    file=paste0(path,basename,"_epochtimes.csv"),    quote=FALSE, row.names=FALSE, col.names="epoch_times",        sep=",")
    write.table(sim$trajgrid,       file=paste0(path,basename,"_popsize.csv"),       quote=FALSE, row.names=FALSE, col.names=c("t","popsize"),     sep=",")
    write.table(sim$samp_intensity, file=paste0(path,basename,"_sampintensity.csv"), quote=FALSE, row.names=FALSE, col.names="sampling_intensity", sep=",")
    write.table(sim$timeshift,      file=paste0(path,basename,"_timeshift.csv"),     quote=FALSE, row.names=FALSE, col.names="first_sample",       sep=",")
    #write.tree(sim$phylo, file=paste0(path,basename,".tree"))
}

#' Save simulation result in a variety of formats
#' 
#' @param sanitise_times If true
save_simulation <- function(simresult, basename, path, gridsize=1000,
                            sanitise_times=TRUE, RData=TRUE, csv=TRUE,
                            newick=TRUE, json=TRUE) {
  
  dir.create(path, showWarnings=FALSE, recursive=TRUE)
  
  if (is.null(names(simresult))) {
    simresult <- lapply(simresult, prepare_simulation,
                        gridsize=gridsize,
                        sanitise_times=sanitise_times,
                        sanity_check=TRUE)
  } else {
    simresult <- prepare_simulation(simresult,
                                    gridsize=gridsize,
                                    sanitise_times=sanitise_times,
                                    sanity_check=TRUE)
  }
  
  if (RData) {
    save(simresult, file=paste0(path, basename, ".RData"))
  }
  
  if (csv) {
    if (is.null(names(simresult))) {
      i <- 0
      lapply(simresult, function(x) {
        save_simulation_csv(x, path, paste(basename, i, sep="_"))
        i <<- i + 1
      })
    } else {
      save_simulation_csv(simresult, path, basename)
    }
  }
  
  # define JSON-safe version in both cases
  if (is.null(names(simresult))) {
    simresultJ <- lapply(simresult, function(x) {
      x$phylo <- NULL
      x
    })
  } else {
    simresultJ <- simresult
    simresultJ$phylo <- NULL
  }
  
  if (newick) {
    trees_file <- paste0(path, basename, ".trees")
    if (file.exists(trees_file)) file.remove(trees_file)
    
    if (is.null(names(simresult))) {
      trees <- lapply(simresult, function(x) x$phylo)
      invisible(lapply(trees, write.tree,
                       file=trees_file,
                       append=TRUE))
    } else {
      write.tree(simresult$phylo, file=paste0(path, basename, ".tree"))
    }
  }
  
  if (json) {
    simresultJSON <- jsonlite::toJSON(simresultJ)
    jsonfile <- file(paste0(path, basename, ".json"), "w")
    cat(simresultJSON, file=jsonfile)
    close(jsonfile)
  }
  
  return(simresult)
}


simulation_example <- function() {
    
    # Example 1 - single simulation trajectory
    trajname <- "expgrowth_fast"
    traj <- get_trajectory(trajname)
    sim  <- simulate_genealogy(traj)
    save_simulation(sim, basename="exp_growth_test", path="../results/test/", RData=TRUE, csv=TRUE)
     
    
    # Example 2 - 100 simulation trajectories
    #traj <- lapply(rep(trajname, 100), get_trajectory)
    #sims <- lapply(traj, simulate_genealogy)
    #save_simulation(sims, basename=trajname, path="../results/test/", RData=TRUE, csv=TRUE)
    
}

