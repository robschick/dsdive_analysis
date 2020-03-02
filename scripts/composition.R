# parallel tools
require(Rmpi, lib.loc = c('.', .libPaths()))
require(snow, lib.loc = c('.', .libPaths()))
# modeling tools
library(dsdive, lib.loc = c('.', .libPaths()))


#
# Set up environment
#

# get or create cluster
cl = getMPIcluster()
if(is.null(cl)) {
  cl = makeCluster(spec = parallel::detectCores() - 1, type = 'SOCK')
}

# initialize RNG streams across cluster
parallel::clusterSetRNGStream(cl = cl, NULL)

# load modeling tools on nodes 
clusterEvalQ(cl = cl, expr = library(dsdive, lib.loc = c(.libPaths(), '.')))


#
# load data and prior parameters for stage
#

# identify dive files
obs.path = file.path('data')
dives = dir(path = obs.path, pattern = 'dive', full.names = TRUE)
depths = dir(path = obs.path, pattern = 'depths', full.names = TRUE)

# storage for observed dives
dives.obs = vector('list', length(dives))

# load dives 
for(i in 1:length(dives)) {
  
  d = read.csv(file = dives[i], header = TRUE)
  db = read.csv(file = depths[i], header = TRUE)
  
  # center dive times
  d$times = d$times - d$times[1]
  
  dive = as.list(d)
  class(dive) = 'dsobs'
  
  # save
  dives.obs[[i]] = list(
    dive = dive,
    depth.bins = db
  )
  
}


# clean up 
rm(d, db, dive, i, depths, dives, obs.path)


#
# extract lists
#

depth.bins = dives.obs[[1]]$depth.bins
dives.obs.list = lapply(dives.obs, function(d) d$dive)

tstep = diff(dives.obs[[1]]$dive$times[1:2])


#
# load posterior samples locally, and on nodes
#

load(file.path('fit', 'trace_theta.RData'))
clusterEvalQ(cl, load(file.path('fit', 'trace_theta.RData')))



#
# impute dives
#

# export additional data to nodes
clusterExport(cl = cl, c('dives.obs', 'depth.bins', 'tstep'))

# partition posterior samples to distribute across nodes
dive.nodes = balance_inds(n.inds = nrow(trace), n.partitions = length(cl))

# create output directories
for(dive.id in 1:length(dives.obs)) {
  dir.create(path = file.path('fit', 'imputations', 
                              paste('dive', dive.id, sep='')), recursive = TRUE)
}

# impute dives in batches (node:dive:sample)
dives.imputed = clusterApply(cl = cl, x = dive.nodes, fun = function(inds) {

  # loop over posterior samples for model parameters
  samples = lapply(inds, function(ind) {
    
    # extract posterior sample for model parameters
    beta = trace[ind, 1:2]
    lambda = trace[ind,3:5]
    
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
    x = lapply(1:length(dives.obs), function(dive.id) {
      d = dives.obs[[dive.id]]$dive
      dsdive.impute(
        depths = d$depths, times = d$times, 
        t.stages = trace.t.stages[[ind]][[dive.id]], rate.unif = rate.unif, 
        P.raw = P.raw, P.tx = P.tx, n.bins = nrow(depth.bins), max.tx = 100)
    })
  })
  
  # save outputs by dive id
  for(dive.id in 1:length(dives.obs)) {
    imputed = lapply(samples, function(s) { s[[dive.id]] })
    o = file.path('fit', 'imputations', 
                  paste('dive', dive.id, sep=''))
    save(imputed, 
         file = file.path(o, paste('dive', dive.id, 'imputed', 
                                   inds[1], inds[length(inds)], 
                                   '.RData', sep = '_')))
  }
  
  # # list structure: dive.id : dive.sample
  # flattened = vector('list', length(dives.obs))
  # for(dive.id in 1:length(dives.obs)) {
  #   flattened[[dive.id]] = lapply(samples, function(s) {
  #     s[[dive.id]]
  #   })
  # }
  
  # return 1 for success
  1
})