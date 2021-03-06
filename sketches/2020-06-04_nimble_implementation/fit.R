# configuration tools
library(composr, lib.loc = c('singularity/libs', .libPaths()))
library(yaml, lib.loc = c('singularity/libs', .libPaths()))
# modeling tools
library(dsdive, lib.loc = c('singularity/libs', .libPaths()))
library(MASS)

# clear workspace
rm(list = ls())


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

# groups = list(
#   data = 'sim_tyack_more_known_end_30',
#   observation_model = 'exact_systematic',
#   priors = 'tyack_simulation_priors',
#   sampler = 'prod_restart',
#   subset = 'all_dives',
#   validation= 'no_validation'
# )

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
source(file.path('scripts', 'utils', '85pct_rule.R'))


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

t.stages.list = lapply(dives.obs.list, function(d) {
  seq(from = d$times[1], to = d$times[length(d$times)], length.out = 4)[2:3]
})


#
# select dives for fitting
#

fit.inds = fit.ind.fn(dives.obs = dives.obs,
                      duration_min = cfg$subset$duration_min,
                      duration_max = cfg$subset$duration_max,
                      holdout = cfg$validation$holdout,
                      seed = cfg$validation$seed,
                      holdout_prop = cfg$validation$proportion,
                      bin_start_max = cfg$subset$bin_start_max,
                      bin_end_max = cfg$subset$bin_end_max)


#
# initial parameters
#

if(cfg$sampler$restart) {
  # load existing output
  load(file.path(out.dir, cfg$base_names$fit))
  load(file.path(out.dir, cfg$base_names$fit_inds))
  # reconfigure remaining samples to draw
  it = nrow(state$theta)
  cfg$sampler$iterations = cfg$sampler$iterations - it
  # extract last params
  params = list(
    beta = state$theta[1:2],
    lambda = state$theta[3:5]
  )
  t.stages = state$trace.t.stages[[it]]
  offsets = state$trace.offsets[it,]
  offsets.tf = state$trace.offsets.tf[it,]
  # label state for merging
  state.bak = state
} else {
  params = list(
    beta = c(.95, .05),
    lambda = c(1.25, .3, .5)
  )
  t.stages = t.stages.list[fit.inds$fit]
  n = length(t.stages)
  offsets = rep(0, n)
  offsets.tf = rep(0, n)
}


#
# Specify priors
#

beta.param = function(mu, sd) {
  # Compute shape1 and shape2 parameters for beta distribution when the
  # distribution is specified via its mean and standard deviation
  v = sd^2
  v.partial = (v + mu * (mu - 1)) / v
  params = c(shape1 = - mu * v.partial, shape2 = (mu - 1) * v.partial)
  if(any(params<0)) {
    stop('No beta distribution may have this combination of mu and sd.')
  } else {
    params
  }
}

gamma.param = function(mu, sd) {
  # Compute shape and rate parameters for gamma distribution when the
  # distribution is specified via its mean and standard deviation
  v = sd^2
  params = c(shape = mu^2/v, rate = mu/v)
  if(any(params<0)) {
    stop('No gamma distribution may have this combination of mu and sd.')
  } else {
    params
  }
}

mu = cfg$priors$descent_preference$mu
b = cfg$priors$descent_preference$shape2
pi1.prior = c(mu/(1-mu) * b, b)

mu = 1-cfg$priors$ascent_preference$mu
b = cfg$priors$ascent_preference$shape1
pi2.prior = rev(c(mu/(1-mu) * b, b))

lambda1.prior = gamma.param(mu = cfg$priors$descent_speed$mu,
                            sd = cfg$priors$descent_speed$sd)
lambda2.prior = gamma.param(mu = cfg$priors$forage_speed$mu,
                            sd = cfg$priors$forage_speed$sd)
lambda3.prior = gamma.param(mu = cfg$priors$ascent_speed$mu,
                            sd = cfg$priors$ascent_speed$sd)

times.stages.est = times.stages(dives.obs = dives.obs[fit.inds$fit])

if(grepl(pattern = 'simulation', x = cfg$priors$name)) {
  # load stage transition time priors for simulations from disk
  load(file.path(cfg$data$path, '..', 'params', 'params.RData'))
  T1.prior = list(estimate = params$T1.params)
  T2.prior = list(estimate = params$T2.params)
  # convert scale of prior parameters from minutes to seconds
  T1.prior.params = T1.prior$estimate / c(1, 60)
  T2.prior.params = T2.prior$estimate / c(1, 60)
} else {

  if(identical(cfg$priors$stage1_tx, 'empirical')) {
    # use 85% rule to determine stage transition time priors
    T1.prior = fitdistr(x = times.stages.est$sub.time.min, densfun = 'gamma')
    # convert scale of prior parameters from minutes to seconds
    T1.prior.params = T1.prior$estimate / c(1, 60)
  } else {
    T1.prior.params = gamma.param(mu = cfg$priors$stage1_tx$mu,
                                  sd = cfg$priors$stage1_tx$sd)
  }

  if(identical(cfg$priors$stage2_tx, 'empirical')) {
    # use 85% rule to determine stage transition time priors
    T2.prior = fitdistr(x = times.stages.est$bottom.time.min, densfun = 'gamma')
    # convert scale of prior parameters from minutes to seconds
    T2.prior.params = T2.prior$estimate / c(1, 60)
  } else {
    T2.prior.params = gamma.param(mu = cfg$priors$stage2_tx$mu,
                                  sd = cfg$priors$stage2_tx$sd)
  }

}


#
# flatten format for nimble
#

depths = do.call(c, lapply(dives.obs.list[fit.inds$fit], function(d) d$depths))

times = do.call(c, lapply(dives.obs.list[fit.inds$fit], function(d) d$times))

start.inds = c(1, 1 + cumsum(
  do.call(c, lapply(dives.obs.list[fit.inds$fit], function(d) length(d$depths)))
))

n.dives = length(start.inds) - 1


#
# construct nimble model, etc.
#

source('sketches/2020-06-04_nimble_implementation/methods/nimble_tools.R')