# configuration tools
library(composr, lib.loc = c('singularity/libs', .libPaths()))
library(yaml, lib.loc = c('singularity/libs', .libPaths()))
# modeling tools
library(dsdive, lib.loc = c('singularity/libs', .libPaths()))
library(MASS)


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

params = list(
  beta = c(.95, .05),
  lambda = c(1.25, .3, .5)
)


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

# use 85% max depth rule to determine time in stages
times.stages = do.call(rbind, lapply(dives.obs, function(d) {
  # extract depths
  depths = d$depth.bins$center[d$dive$depths]
  # find max depth, and stage threshold
  max.depth = max(depths)
  stage.thresh = .85 * max.depth
  # compute observed stage vector
  bottom.range = range(which(depths >= stage.thresh))
  if(length(unique(bottom.range))==1) {
    bottom.range[2] = bottom.range[1] + 1
  }
  stages = stagevec(length.out = length(depths), breaks = bottom.range)
  # linearly interpolate to find stage transition times
  t.inds = which(diff(stages)==1)
  t.stages = sapply(t.inds, function(ind) {
    # get start and end times/depths
    d0 = depths[ind]
    t0 = d$dive$times[ind]
    df = depths[ind+1]
    tf = d$dive$times[ind+1]
    # compute time at which stage.thresh is crossed
    if(df==d0) {
      mean(c(t0,tf))
    } else {
      (stage.thresh - d0)/(df-d0) * (tf-t0) + t0
    }
  })
  # return results
  data.frame(sub.time.min = t.stages[1]/60, 
             bottom.time.min = diff(t.stages)/60)
}))


if(grepl(pattern = 'simulation', x = cfg$priors$name)) {
  # load stage transition time priors for simulations from disk
  load(file.path(cfg$data$path, '..', 'params', 'params.RData'))
  T1.prior = list(estimate = params$T1.params)
  T2.prior = list(estimate = params$T2.params)
} else {
  # use 85% rule to determine stage transition time priors
  T1.prior = fitdistr(x = times.stages$sub.time.min, densfun = 'gamma')
  T2.prior = fitdistr(x = times.stages$bottom.time.min, densfun = 'gamma')
}

# # plot priors
# curve(dbeta(x = x, shape1 = pi1.prior[1], shape2 = pi1.prior[2]), 
#       xlab = expression(pi[1]), ylab = expression(f(pi[1])))
# curve(dbeta(x = x, shape1 = pi2.prior[1], shape2 = pi2.prior[2]),
#       xlab = expression(pi[2]), ylab = expression(f(pi[2])))
# curve(dgamma(x = x, shape = lambda1.prior[1], rate = lambda1.prior[2]), 
#       xlab = expression(lambda[1]), ylab = expression(f(lambda[1])), 
#       from = 0, to = 4)
# curve(dgamma(x = x, shape = lambda2.prior[1], rate = lambda2.prior[2]), 
#       xlab = expression(lambda[2]), ylab = expression(f(lambda[2])), 
#       from = 0, to = 1)
# curve(dgamma(x = x, shape = lambda3.prior[1], rate = lambda3.prior[2]), 
#       xlab = expression(lambda[3]), ylab = expression(f(lambda[3])), 
#       from = 0, to = 2)
# curve(dgamma(x = x, shape = T1.prior$estimate[1], rate = T1.prior$estimate[2]),
#       xlab = expression(T^(1)~~(min.)), ylab = expression(f(T^(1))),
#       from = 0, to = 25)
# curve(dgamma(x = x, shape = T2.prior$estimate[1], rate = T2.prior$estimate[2]),
#       xlab = expression(T^(2)~~(min.)), ylab = expression(f(T^(2))),
#       from = 0, to = 60)

# convert scale of prior parameters from minutes to seconds
T1.prior.params = T1.prior$estimate / c(1, 60)
T2.prior.params = T2.prior$estimate / c(1, 60)


#
# gibbs sample
#

save(T1.prior.params, T2.prior.params, 
     file = file.path(out.dir, cfg$base_names$stage_priors))
save(fit.inds, file = file.path(out.dir, cfg$base_names$fit_inds))
write_yaml(cfg, file = file.path(out.dir, 'cfg.yaml'))

dump.state = function(state) {
  save.time = date()
  save(state, save.time, params, file = file.path(out.dir, cfg$base_names$fit))
}

# Save crash info to file last.dump.rda
dump_on_error <- function() {
  dump.frames(dumpto = file.path(out.dir, 'last.dump'), to.file = TRUE)
}
options(error = dump_on_error)

fit = dsdive.gibbs.obs(
  dsobs.list = dives.obs.list[fit.inds$fit], 
  t.stages.list = t.stages.list[fit.inds$fit], 
  beta.init = params$beta, lambda.init = params$lambda, 
  verbose = cfg$sampler$verbose, 
  maxit = cfg$sampler$iterations, checkpoint.fn = dump.state, 
  checkpoint.interval = cfg$sampler$checkpoint_interval, 
  pi1.prior = pi1.prior, pi2.prior = pi2.prior, lambda1.prior = lambda1.prior, 
  lambda2.prior = lambda2.prior, lambda3.prior = lambda3.prior, 
  tstep = cfg$data$tstep, 
  depth.bins = depth.bins, T1.prior.params = T1.prior.params, 
  T2.prior.params = T2.prior.params, max.width = 100, max.width.offset = 60, 
  t0.prior.params = unlist(cfg$observation_model$parameters))

options(error = NULL)

if(exists('fit')) {
  dump.state(state = fit)
}
