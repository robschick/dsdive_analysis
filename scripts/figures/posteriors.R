# configuration tools
library(composr)
library(yaml)
# parallel tools
require(Rmpi, lib.loc = c('.', .libPaths()))
require(snow, lib.loc = c('.', .libPaths()))
# modeling tools
library(dsdive, lib.loc = c('.', .libPaths()))
# validation tools
library(dplyr)
library(ggplot2)
library(ggthemes)
library(coda)


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

cfg = compose_cfg(file = 'output/sim_deeper_300/all_dives/no_validation/exact_systematic/simulation_priors/cfg.yaml')

# output paths
out.dir = file.path(cfg$base_paths$fit, cfg$data$name, cfg$subset$name, 
                    cfg$validation$name, cfg$observation_model$name, 
                    cfg$priors$name)

# load data and utility functions
source(file.path('scripts', 'utils', 'datafns.R'))


#
# load data
#

dives.obs = dives.load(path = cfg$data$path, 
                       dive_pattern = cfg$data$file_patterns$dive,
                       depth_pattern = cfg$data$file_patterns$depths)


#
# load posterior output
#

load(file.path(out.dir, cfg$base_names$fit))

m = mcmc(state$theta[-(1:cfg$sampler$burn),])


#
# posterior summaries
#

if(is.null(cfg$sub_paths$posteriors)) {
  defaults = compose_cfg(file = 'conf/config.yaml')
  cfg$sub_paths$posteriors = defaults$sub_paths$posteriors
}

o = file.path(out.dir, cfg$sub_paths$figures, cfg$sub_paths$posteriors)

dir.create(o, recursive = TRUE)


# traceplots and density estimates
for(i in 1:ncol(m)) {
  nom = colnames(m)[i]
  png(file.path(o, paste(nom, '.png', sep='')), width = 480*2, height = 480)
  plot(m[,i], main = colnames(m)[i])
  dev.off()
}

# text-based summaries
sink(file.path(o, 'diving_parameters_summary.txt'))
summary(m)
cat('\nHPD Intervals\n\n')
HPDinterval(m)
cat('\nEffective sample size\n\n')
effectiveSize(m)
sink()
