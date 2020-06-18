# configuration tools
library(composr, lib.loc = c('singularity/libs', .libPaths()))
library(yaml, lib.loc = c('singularity/libs', .libPaths()))
# modeling tools
library(dsdive, lib.loc = c('singularity/libs', .libPaths()))
library(MASS)


#
# set up environment
#

# clear workspace
rm(list = ls())

# load data and utility functions
source(file.path('scripts', 'utils', '85pct_rule.R'))


#
# load data
#

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

# dive observations, priors, other data structures, and base-level 
# inclusion/exclusion from model fitting
dives.processed = readRDS(file.path('sketches', '2020-06-17_nimble_fulltag',
                                    'dives_processed.rds'))

# depth bins
depth.bins = read.csv(file.path('data', 'depth_template.csv'))


#
# extract basic information
#

n.dives = length(dives.processed$dive.flags)


#
# initial parameters
#


params.deep = list(
  pi = c(.95, .5, .05),
  lambda = c(1.25, .3, .5)
)

params.shallow = list(
  pi = c(.95, .05),
  lambda = c(.6, .6)
)


# stage transition times
Tmat = do.call(rbind, lapply(1:n.dives, function(dive.id) {
  
  T.init = rep(0,4)
  
  # find the dive record in the processed information
  ranges_row = which(dives.processed$dive.ranges$dive.id == dive.id)
  
  if(length(ranges_row) > 0) {
    
    # set start time random effect
    start_row = which(dives.processed$endpoint.inds$dive_start == dive.id)
    T.init[1] = mean(unlist(
      dives.processed$endpoint.inds[start_row, c('t_lwr', 't_upr')]
    ))
    
    # set end time random effect
    end_row = which(dives.processed$endpoint.inds$dive_end == dive.id)
    T.init[4] = mean(unlist(
      dives.processed$endpoint.inds[end_row, c('t_lwr', 't_upr')]
    ))
    
    # get ranges of depth bins associated with dive
    start_ind = dives.processed$dive.ranges$start.ind[ranges_row]
    end_ind = dives.processed$dive.ranges$end.ind[ranges_row]
    
    # set initial stage transition times...
    if(dives.processed$dive.ranges$type[ranges_row] == 1) {
      # ...for deep dives
      
      # use 85% rule to estimate stage durations
      stage_durations = times.stages(list(list(
        dive = list(
          depths = dives.processed$depths[start_ind:end_ind],
          times = dives.processed$times[start_ind:end_ind] - 
            dives.processed$times[start_ind]
        ),
        depth.bins = depth.bins
      ))) * 60
      
      # convert to stage transition times
      T.init[2] = T.init[1] + stage_durations$sub.time.min
      T.init[3] = T.init[2] + stage_durations$bottom.time.min
      
    } else {
      # ...for shallow dives
  
      # stop descent halfway through dive
      T.init[2] = mean(T.init[c(1,4)])
      # copy dive-end time
      T.init[3] = T.init[4]
    }
    
  }
  
  T.init
}))

colnames(Tmat) = paste('T', 0:3, sep='')

# stage durations
xi = t(apply(Tmat, 1, function(r) {
  diff(r)
}))


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

lambda1.prior.shallow = gamma.param(mu = cfg$priors$descent_speed$mu/2,
                            sd = cfg$priors$descent_speed$sd)
lambda2.prior.shallow = gamma.param(mu = cfg$priors$forage_speed$mu/2,
                            sd = cfg$priors$forage_speed$sd)


#
# empirically determine stage duration priors for deep dives
#

deep.inds = which(sapply(1:n.dives, function(dive.id) {
  range_row = which(dives.processed$dive.ranges$dive.id == dive.id)
  type = dives.processed$dive.ranges$type[range_row]
  ifelse(length(type)==0, 0, type)
}) == 1)

shallow.inds = which(sapply(1:n.dives, function(dive.id) {
  range_row = which(dives.processed$dive.ranges$dive.id == dive.id)
  type = dives.processed$dive.ranges$type[range_row]
  ifelse(length(type)==0, 0, type)
}) == 2)


T1.prior = fitdistr(x = apply(Tmat[deep.inds,c('T0', 'T1')], 1, diff), 
                    densfun = 'gamma')$estimate

T2.prior = fitdistr(x = apply(Tmat[deep.inds,c('T1', 'T2')], 1, diff), 
                    densfun = 'gamma')$estimate

T1.prior.shallow = fitdistr(x = apply(Tmat[shallow.inds,c('T0', 'T1')], 1, 
                                      diff), 
                            densfun = 'gamma')$estimate


#
# construct nimble model, etc.
#

source('sketches/2020-06-17_nimble_fulltag/methods/nimble_tools.R')
