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

groups = list(
  data = 'zc84_800',
  observation_model = 'uniform_systematic',
  priors = 'tyack_priors',
  sampler = 'prod',
  subset = 'all_dives',
  validation= 'holdout_half'
)

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
# extract lists
#

depth.bins = dives.obs[[1]]$depth.bins
dives.obs.list = lapply(dives.obs, function(d) d$dive)

tstep = diff(dives.obs[[1]]$dive$times[1:2])


#
# load posterior samples locally, and on nodes
#

load(file.path(out.dir, cfg$base_names$fit))
clusterExport(cl, c('cfg', 'out.dir'))
clusterEvalQ(cl, load(file.path(out.dir, cfg$base_names$fit)))

if(is.null(state$trace)) {
  state$trace = state$theta
  state$theta = NULL
}

if(!('trace.offset' %in% names(state))) {
  state$trace.offset = state$trace.offsets
  state$trace.offsets = NULL
}

clusterExport(cl, 'state')

#
# sample dives from posterior predictive distribution
#

# create and empty output directories
pred.dir = file.path(out.dir, cfg$sub_paths$posterior_predictions)
dir.create(path = pred.dir, recursive = TRUE)
file.remove(dir(pred.dir, full.names = TRUE))

# export additional data to nodes
clusterExport(cl = cl, c('dives.obs', 'depth.bins', 'tstep', 'pred.dir'))

# partition posterior samples to distribute across nodes
dive.nodes = balance_inds(n.inds = nrow(state$trace), n.partitions = length(cl))

# impute dives in batches (node:dive:sample)
dives.imputed = clusterApply(cl = cl, x = dive.nodes, fun = function(inds) {

  attach(state)
  
  # loop over posterior samples for model parameters
  samples = lapply(inds, function(ind) {
    
    # extract posterior sample for model parameters
    beta = trace[ind, 1:2]
    lambda = trace[ind,3:5]
    
    # sample dive stage durations
    stages.dur = diff(c(0, unlist(sample(x = trace.t.stages[[ind]], size = 1))))
    
    # sample dive offset
    offset = sample(x = trace.offset[ind,], size = 1)
    
    # sample and return dive
    d = dsdive.fwdsample.dive(depth.bins = depth.bins, beta = beta, 
                              lambda = lambda, t0 = 0, steps.max = 1e3, 
                              T1 = stages.dur[1], T2 = stages.dur[2])
    d$t0.offset = offset
    d
  })
  
  # construct output file name
  f = gsub(pattern = cfg$base_names$posterior_predictions$index_pattern, 
           replacement = paste(c(inds[1], inds[length(inds)]), collapse = '_'), 
           x = cfg$base_names$posterior_predictions$file)
  
  # save outputs by posterior sample number
  save(samples, file = file.path(pred.dir, f))
  
  # return 1 for success
  1
})

stopCluster(cl)
