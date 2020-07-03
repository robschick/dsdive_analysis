# configuration tools
library(composr, lib.loc = c('singularity/libs', .libPaths()))
library(yaml, lib.loc = c('singularity/libs', .libPaths()))
# modeling tools
library(dsdive, lib.loc = c('singularity/libs', .libPaths()))
library(MASS)
# plotting tools
library(ggplot2)
library(ggthemes)

# clear workspace
rm(list = ls())


#
# Set up nodes
#

# # get or create cluster
# cl = getMPIcluster()
# if(is.null(cl)) {
#   cl = makeCluster(spec = parallel::detectCores() - 1, type = 'SOCK')
# }
# 
# # get cluster size
# nodes = length(cl)
# 
# # initialize RNG streams across cluster
# parallel::clusterSetRNGStream(cl, NULL)
# 
# # load analysis package on nodes
# clusterEvalQ(cl, library(dsdive, lib.loc = c('singularity/libs', .libPaths())))


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
# EDA for stages
#

times.stages.est = times.stages(dives.obs = dives.obs)

ggplot(times.stages.est, aes(x=sub.time.min, y=bottom.time.min)) + 
  geom_point() + 
  xlab('Descent duration (min)') + 
  ylab('Foraging duration (min)') + 
  ggtitle('EDA via 85% rule') + 
  geom_smooth(method = 'lm') + 
  theme_few()


fit = lm(bottom.time.min ~ sub.time.min, times.stages.est)

summary(fit)

# pretty decent residuals
plot(fit)

cor(times.stages.est)