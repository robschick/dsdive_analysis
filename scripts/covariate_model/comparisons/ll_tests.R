# configuration tools
library(composr)
library(yaml)
# modeling tools
library(dsdive, lib.loc = c('.', .libPaths()))
# validation tools
library(dplyr)
library(ggplot2)
library(ggthemes)
library(tidyr)
library(forcats)
library(coda)

# clear workspace
rm(list = ls())

#
# Set up nodes
#

# get or create cluster
cl = getMPIcluster()
if(is.null(cl)) {
  # cl = makeCluster(spec = detectCores(), type = 'SOCK')
  cl = makeCluster(spec = 4, type = 'SOCK')
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


groups.compare = list(
  intercept = list(
    data = 'zc84_800_covariates',
    observation_model = 'uniform_systematic',
    priors = 'tyack_cov_intercept_priors',
    sampler = 'prod',
    subset = 'all_dives',
    validation= 'holdout_half'
  ),
  daynight = list(
    data = 'zc84_800_covariates',
    observation_model = 'uniform_systematic',
    priors = 'tyack_cov_daynight_priors',
    sampler = 'prod',
    subset = 'all_dives',
    validation= 'holdout_half'
  ),
  duration = list(
    data = 'zc84_800_covariates',
    observation_model = 'uniform_systematic',
    priors = 'tyack_cov_duration_priors',
    sampler = 'prod',
    subset = 'all_dives',
    validation= 'holdout_half'
  ),
  surfactivity = list(
    data = 'zc84_800_covariates',
    observation_model = 'uniform_systematic',
    priors = 'tyack_cov_surfactivity_priors',
    sampler = 'prod',
    subset = 'all_dives',
    validation= 'holdout_half'
  ),
  depthvar = list(
    data = 'zc84_800_covariates',
    observation_model = 'uniform_systematic',
    priors = 'tyack_cov_depthvar_priors',
    sampler = 'prod',
    subset = 'all_dives',
    validation= 'holdout_half'
  ),
  full = list(
    data = 'zc84_800_covariates',
    observation_model = 'uniform_systematic',
    priors = 'tyack_cov_priors2',
    sampler = 'prod',
    subset = 'all_dives',
    validation= 'holdout_half'
  )
)

# build configuration
cfg.list = lapply(groups.compare, function(g) {
  compose_cfg(file = file.path('conf', 'config.yaml'), groups = g)
})
rm(groups.compare)

names.group = names(cfg.list)

# output paths
out.dir.list = sapply(cfg.list, function(cfg) {
  file.path(cfg$base_paths$fit, cfg$data$name, cfg$subset$name, 
            cfg$validation$name, cfg$observation_model$name, 
            cfg$priors$name)
})

# load utility functions
source(file.path('scripts', 'utils', 'LoadToEnvironment.R'))


# Capitalize the first letter in a word string
CapStr <- function(y) {
  sapply(y, function(y) {
    c <- strsplit(y, " ")[[1]]
    paste(toupper(substring(c, 1,1)), substring(c, 2),
          sep="", collapse=" ")
  })
}


#
# load validation files to a list of environments 
#

LoadValidationData = function(filename) {
  lapply(1:length(cfg.list), function(i) {
    cfg = cfg.list[[i]]
    o = out.dir.list[i]
    LoadToEnvironment(file.path(o, filename))
  })
}


#
# likelihood evaluated at posterior mean 
#




#
# compare depth-by-time validations
#

valcmp = LoadValidationData('sampler_trace.RData')


ll.post = sapply(valcmp, function(env) {
  burn = 1:1e3
  # 
  inds = which(!is.na(env$state$theta$beta1[,1]))
  
  # 
  # theta.mat = cbind(env$state$theta$beta1, env$state$theta$beta2, 
  #                   env$state$theta$alpha1, env$state$theta$alpha2, 
  #                   env$state$theta$alpha3)[inds,][-burn,]
  # 
  
  # dsdive.ob
  
  finite = is.finite(env$state$trace.ll)
  
  c(ll.post = mean(env$state$trace.ll[finite][-burn]))
})

# in-sample model comparison
data.frame(Model = names.group, ll.post = ll.post) %>% arrange(-ll.post)

# in-sample likelihood comparisons (but won't have p-values b/c support 
# violated, and also unsure how to select d.o.f.)

# out-of-sample log scores on validation data (TBD, may be slow)