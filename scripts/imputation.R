# parallel tools
require(Rmpi, lib.loc = c('singularity/libs', .libPaths()))
require(snow, lib.loc = c('singularity/libs', .libPaths()))
# modeling tools
library(dsdive, lib.loc = c('singularity/libs', .libPaths()))
# configuration tools
library(composr, lib.loc = c('singularity/libs', .libPaths()))
library(yaml, lib.loc = c('singularity/libs', .libPaths()))


#
# Set up environment
#

# read in configuration groups
args = commandArgs(TRUE)
if(length(args)>0) {
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
} else {
  groups = NULL
}
rm(args,i)

# build configuration
cfg = compose_cfg(file = file.path('conf', 'config.yaml'), groups = groups)
rm(groups)

# output paths
out.dir = file.path(cfg$base_paths$fit, cfg$data$name, cfg$subset$name, 
                    cfg$validation$name, cfg$observation_model$name, 
                    cfg$priors$name)
dir.create(out.dir, recursive = TRUE)

# load data and utility functions
source(file.path('scripts', 'utils', 'datafns.R'))


#
# Set up parallel environment
#

# get or create cluster
cl = getMPIcluster()
if(is.null(cl)) {
  cl = makeCluster(spec = parallel::detectCores() - 1, type = 'SOCK')
}

# initialize RNG streams across cluster
parallel::clusterSetRNGStream(cl = cl, NULL)

# load modeling tools on nodes 
clusterEvalQ(cl = cl, expr = library(dsdive, lib.loc = c('singularity/libs', 
                                                         .libPaths())))


#
# load data
#

dives.obs = dives.load(path = cfg$data$path, 
                       dive_pattern = cfg$data$file_patterns$dive,
                       depth_pattern = cfg$data$file_patterns$depths)

load(file.path(out.dir, cfg$base_names$fit_inds))


#
# extract lists
#

depth.bins = dives.obs[[1]]$depth.bins
dives.obs.list = lapply(dives.obs, function(d) d$dive)

tstep = diff(dives.obs[[1]]$dive$times[1:2])


#
# load posterior samples locally, and on nodes
#

load(file.path(out.dir, cfg$base_names$fit))
clusterExport(cl, c('out.dir', 'cfg'))
clusterEvalQ(cl, load(file.path(out.dir, cfg$base_names$fit)))


#
# impute dives
#

# export additional data to nodes
clusterExport(cl = cl, c('dives.obs', 'depth.bins', 'tstep', 'fit.inds'))

# partition posterior samples to distribute across nodes; use more partitions 
# than clusters as a way to checkpoint output and control memory demands
dive.nodes = balance_inds(n.inds = nrow(state$theta), 
                          n.partitions = 4*length(cl))

if(is.null(cfg$sub_paths$imputations)) {
  defaults = compose_cfg(file = 'conf/config.yaml')
  cfg$sub_paths$imputations = defaults$sub_paths$imputations
}


# create output directories
for(dive.id in fit.inds$fit) {
  dir.create(path = file.path(out.dir, cfg$sub_paths$imputations, 
                              paste('dive', dive.id, sep='')), recursive = TRUE)
}

# impute dives in batches (node:dive:sample)
dives.imputed = clusterApply(cl = cl, x = dive.nodes, fun = function(inds) {

  # loop over posterior samples for model parameters
  samples = lapply(inds, function(ind) {
    
    # extract posterior sample for model parameters
    beta = state$theta[ind, 1:2]
    lambda = state$theta[ind,3:5]
    
    # compute uniformized transition rate
    rate.unif = max(outer(lambda, 2 * depth.bins[,2], '/'))
    
    # probability transition matrix for observations
    P.raw = lapply(1:3, function(s) {
      dsdive.obstx.matrix(depth.bins = depth.bins, beta = beta, 
                          lambda = lambda, s0 = s, tstep = tstep, 
                          include.raw = TRUE)
    })
    
    # probability transition matrix for uniformized DTMC
    P.tx = lapply(1:3, function(s) {
      dsdive.tx.matrix.uniformized(depth.bins = depth.bins, beta = beta, 
                                   lambda = lambda, s0 = s, 
                                   rate.uniformized = rate.unif)
    })
    
    # loop over dives
    x = lapply(1:length(fit.inds$fit), function(dive.ind) {
      # convert index to dive id
      dive.id = fit.inds$fit[dive.ind]
      # extract dive observation
      d = dives.obs[[dive.id]]$dive
      # extract random effects
      t.stages = state$trace.t.stages[[ind]][[dive.ind]]
      offset = state$trace.offsets[ind, dive.ind]
      # align dive observations with offset
      d.aligned = dsdive.align.obs(depths = d$depths, times = d$times, 
                                   t.stages = t.stages, offset = offset)
      # impute dive
      d.imputed.aligned = dsdive.impute(
        depths = d.aligned$depths, times = d.aligned$times, t.stages = t.stages, 
        rate.unif = rate.unif, P.raw = P.raw, P.tx = P.tx, 
        n.bins = nrow(depth.bins), max.tx = 100)
      # align imputed dive to observation time-frame
      d.imputed.aligned$times = d.imputed.aligned$times + offset
      if(offset > 0) {
        d.imputed.aligned$times[1] = 0
        d.imputed.aligned$durations[1] = d.imputed.aligned$durations[1] + offset
      }
      # return result
      d.imputed.aligned
    })
  })
  
  # save outputs by dive id
  for(dive.ind in 1:length(fit.inds$fit)) {
    dive.id = fit.inds$fit[dive.ind]
    imputed = lapply(samples, function(s) { s[[dive.ind]] })
    o = file.path(out.dir, cfg$sub_paths$imputations, 
                  paste('dive', dive.id, sep=''))
    save(imputed, 
         file = file.path(o, paste('dive', dive.id, 'imputed', 
                                   inds[1], inds[length(inds)], 
                                   '.RData', sep = '_')))
  }
  
  # return 1 for success
  1
})