# configuration tools
library(composr)
library(yaml)
# parallel tools
require(Rmpi, lib.loc = c('.', .libPaths()))
require(snow, lib.loc = c('.', .libPaths()))
# modeling tools
library(dsdive, lib.loc = c('.', .libPaths()))


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

# groups = list(
#   data = 'zc84_800',
#   observation_model = 'uniform_systematic',
#   priors = 'tyack_priors_fixed_stage',
#   sampler = 'prod',
#   subset = 'all_dives',
#   validation= 'holdout_half'
# )

# build configuration
cfg = compose_cfg(file = file.path('conf', 'config.yaml'), groups = groups)
rm(groups)

# output paths
out.dir = file.path(cfg$base_paths$fit, cfg$data$name, cfg$subset$name, 
                    cfg$validation$name, cfg$observation_model$name, 
                    cfg$priors$name)

# get or create cluster
cl = getMPIcluster()
if(is.null(cl)) {
  cl = makeCluster(spec = parallel::detectCores() - 1, type = 'SOCK')
}

# initialize RNG streams across cluster
parallel::clusterSetRNGStream(cl = cl, NULL)

# load modeling tools on nodes 
clusterEvalQ(cl = cl, expr = library(dsdive, lib.loc = c(.libPaths(), '.')))

# load data and utility functions
source(file.path('scripts', 'utils', 'datafns.R'))


#
# load data
#

dives.obs = dives.load(path = cfg$data$path, 
                       dive_pattern = cfg$data$file_patterns$dive,
                       depth_pattern = cfg$data$file_patterns$depths)


#
# load priors
#

load(file.path(out.dir, cfg$base_names$stage_priors))
load(file.path(out.dir, cfg$base_names$ctmc_priors))


#
# extract lists
#

depth.bins = dives.obs[[1]]$depth.bins
dives.obs.list = lapply(dives.obs, function(d) d$dive)

tstep = diff(dives.obs[[1]]$dive$times[1:2])


#
# sample dives from prior predictive distribution
#

# number of dives to sample
mcit = 5e4

# create and empty output directories
pred.dir = file.path(out.dir, cfg$sub_paths$prior_predictions)
dir.create(path = pred.dir, recursive = TRUE)
file.remove(dir(pred.dir, full.names = TRUE))

# export additional data to nodes
clusterExport(cl = cl, c('dives.obs', 'depth.bins', 'tstep', 'pred.dir',
                         'lambda1.prior', 'lambda2.prior', 'lambda3.prior',
                         'pi1.prior', 'pi2.prior', 'T1.prior.params', 
                         'T2.prior.params', 'cfg', 'out.dir'))

# partition posterior samples to distribute across nodes
dive.nodes = balance_inds(n.inds = mcit, n.partitions = length(cl))

# impute dives in batches (node:dive:sample)
dives.imputed = clusterApply(cl = cl, x = dive.nodes, fun = function(inds) {
  
  # loop over posterior samples for model parameters
  samples = lapply(inds, function(ind) {
    
    # generate prior sample for model parameters
    beta = c(
      rbeta(n = 1, shape1 = pi1.prior[1], shape2 = pi1.prior[2]),
      rbeta(n = 1, shape1 = pi2.prior[1], shape2 = pi2.prior[2])
    )
    lambda = c(
      rgamma(n = 1, shape = lambda1.prior[1], rate = lambda1.prior[2]),
      rgamma(n = 1, shape = lambda2.prior[1], rate = lambda2.prior[2]),
      rgamma(n = 1, shape = lambda3.prior[1], rate = lambda3.prior[2])
    )
    
    # sample dive stage durations
    stages.dur = c(
      rgamma(n = 1, shape = T1.prior.params[1], rate = T1.prior.params[2]),
      rgamma(n = 1, shape = T2.prior.params[1], rate = T2.prior.params[2])
    )
    
    # sample and return dive
    d = dsdive.fwdsample.dive(depth.bins = depth.bins, beta = beta, 
                              lambda = lambda, t0 = 0, steps.max = 1e3, 
                              T1 = stages.dur[1], T2 = stages.dur[2])
    d
  })
  
  # construct output file name
  f = gsub(pattern = cfg$base_names$prior_predictions$index_pattern, 
           replacement = paste(c(inds[1], inds[length(inds)]), collapse = '_'), 
           x = cfg$base_names$prior_predictions$file)
  
  # save outputs by posterior sample number
  save(samples, file = file.path(pred.dir, f))
  
  # return 1 for success
  1
})

stopCluster(cl)
