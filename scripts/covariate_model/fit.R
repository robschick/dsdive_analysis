# configuration tools
library(composr, lib.loc = c('singularity/libs', .libPaths()))
library(yaml, lib.loc = c('singularity/libs', .libPaths()))
# modeling tools
library(dsdive, lib.loc = c('singularity/libs', .libPaths()))
library(MASS)
# # parallel tools
require(parallel, lib.loc = c('singularity/libs', .libPaths()))
require(Rdsm, lib.loc = c('singularity/libs', .libPaths()))
require(snow, lib.loc = c('singularity/libs', .libPaths()))

# clear workspace
rm(list = ls())


#
# Set up nodes
#

# get or create cluster
cl = getMPIcluster()
if(is.null(cl)) {
  cl = makeCluster(spec = detectCores(), type = 'SOCK')
}


# get cluster size
nodes = length(cl)

# initialize RNG streams across cluster
parallel::clusterSetRNGStream(cl, NULL)

# load analysis package and support on nodes
clusterEvalQ(cl, library(dsdive, lib.loc = c('singularity/libs', .libPaths())))
clusterEvalQ(cl, library(Rdsm, lib.loc = c('singularity/libs', .libPaths())))

# Rdsm-initialize cluster
mgrinit(cl)


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
#   data = 'zc84_800_covariates',
#   observation_model = 'uniform_systematic',
#   priors = 'tyack_cov_priors',
#   sampler = 'prod',
#   subset = 'all_dives',
#   validation= 'holdout_half'
# )

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
                       depth_pattern = cfg$data$file_patterns$depths, 
                       covariates = cfg$data$file_patterns$covariates)


#
# extract lists
#

depth.bins = dives.obs$dives[[1]]$depth.bins
dives.obs.list = lapply(dives.obs$dives, function(d) d$dive)

t.stages.list = lapply(dives.obs.list, function(d) {
  seq(from = d$times[1], to = d$times[length(d$times)], length.out = 4)[2:3]
})


#
# select dives for fitting
#

fit.inds = fit.ind.fn(dives.obs = dives.obs$dives, 
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
  # load(file.path(out.dir, cfg$base_names$fit))
  # load(file.path(out.dir, cfg$base_names$fit_inds))
  # # reconfigure remaining samples to draw
  # it = nrow(state$theta)
  # cfg$sampler$iterations = cfg$sampler$iterations - it
  # # extract last params
  # params = list(
  #   beta = state$theta[1:2],
  #   lambda = state$theta[3:5]
  # )
  # t.stages = state$trace.t.stages[[it]]
  # offsets = state$trace.offsets[it,]
  # offsets.tf = state$trace.offsets.tf[it,]
  # # label state for merging
  # state.bak = state
} else {
  params = list(
    beta = list(
      as.numeric(cfg$priors$logit_descent_preference$mu),
      as.numeric(cfg$priors$logit_ascent_preference$mu)
    ),
    alpha = list(
      as.numeric(cfg$priors$log_descent_speed$mu),
      as.numeric(cfg$priors$log_forage_speed$mu),
      as.numeric(cfg$priors$log_ascent_speed$mu)
    )
  )
  t.stages = t.stages.list[fit.inds$fit]
  n = length(t.stages)
  offsets = rep(0, n)
  offsets.tf = rep(0, n)
}


#
# Specify priors and extract formulas
#

pi.formula = list(
  formula(cfg$priors$logit_descent_preference$formula),
  formula(cfg$priors$logit_ascent_preference$formula)
)

lambda.formula = list(
  formula(cfg$priors$log_descent_speed$formula),
  formula(cfg$priors$log_forage_speed$formula),
  formula(cfg$priors$log_ascent_speed$formula)
)

beta.priors = list(
  list(mu = as.numeric(cfg$priors$logit_descent_preference$mu),
       sd = as.numeric(cfg$priors$logit_descent_preference$sd)),
  list(mu = as.numeric(cfg$priors$logit_ascent_preference$mu),
       sd = as.numeric(cfg$priors$logit_ascent_preference$sd))
)

alpha.priors = list(
  list(mu = as.numeric(cfg$priors$log_descent_speed$mu),
       sd = as.numeric(cfg$priors$log_descent_speed$sd)),
  list(mu = as.numeric(cfg$priors$log_forage_speed$mu),
       sd = as.numeric(cfg$priors$log_forage_speed$sd)),
  list(mu = as.numeric(cfg$priors$log_ascent_speed$mu),
       sd = as.numeric(cfg$priors$log_ascent_speed$sd))
)

times.stages.est = times.stages(dives.obs = dives.obs$dives[fit.inds$fit])

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
# gibbs sample
#
  
save(beta.priors, alpha.priors,
     file = file.path(out.dir, cfg$base_names$ctmc_priors))
save(T1.prior.params, T2.prior.params, 
     file = file.path(out.dir, cfg$base_names$stage_priors))
save(fit.inds, file = file.path(out.dir, cfg$base_names$fit_inds))
write_yaml(cfg, file = file.path(out.dir, 'cfg.yaml'))

dump.state = function(state) {
  
  if(cfg$sampler$restart) {
    # state$theta = rbind(state.bak$theta, state$theta)
    # state$trace.t.stages = c(state.bak$trace.t.stages, state$trace.t.stages)
    # state$trace.offsets = rbind(state.bak$trace.offsets, state$trace.offsets)
    # state$trace.offsets.tf = rbind(state.bak$trace.offsets.tf, 
    #                                state$trace.offsets.tf)
  }
  
  save.time = date()
  save(state, save.time, params, file = file.path(out.dir, cfg$base_names$fit))
}

# Save crash info to file last.dump.rda
dump_on_error <- function() {
  dump.frames(dumpto = file.path(out.dir, 'last.dump'), to.file = TRUE, 
              include.GlobalEnv = TRUE)
}
options(error = dump_on_error)

fit = dsdive.gibbs.obs.cov(
  dsobs.list = dives.obs.list[fit.inds$fit], 
  t.stages.list = t.stages, 
  beta.init = params$beta, alpha.init = params$alpha, 
  verbose = cfg$sampler$verbose, 
  maxit = cfg$sampler$iterations, checkpoint.fn = dump.state, 
  checkpoint.interval = cfg$sampler$checkpoint_interval, 
  beta1.prior = beta.priors[[1]], beta2.prior = beta.priors[[2]], 
  alpha1.prior = alpha.priors[[1]], alpha2.prior = alpha.priors[[2]],
  alpha3.prior = alpha.priors[[3]], tstep = cfg$data$tstep, 
  depth.bins = depth.bins, T1.prior.params = T1.prior.params, 
  T2.prior.params = T2.prior.params, max.width = 100, max.width.offset = 60, 
  t0.prior.params = unlist(cfg$observation_model$parameters),
  tf.prior.params = unlist(cfg$observation_model$parameters_tf), 
  offsets = offsets, offsets.tf = offsets.tf, covs = dives.obs$covariates, 
  pi.formula = pi.formula, lambda.formula = lambda.formula,
  warmup = as.numeric(cfg$sampler$warmup), cl = cl, gapprox = TRUE)

options(error = NULL)

if(exists('fit')) {
  dump.state(state = fit)
}

stopCluster(cl)

# clear shared memory artifacts
file.remove(dir(pattern = '.desc'))
