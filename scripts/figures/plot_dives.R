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
library(tidyr)
library(forcats)
library(KSgeneral)


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

cfg = read_yaml(file = 'output/sim_deeper_300/all_dives/no_validation/exact_systematic/simulation_priors/cfg.yaml')

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
# plot dives
#

if(is.null(cfg$sub_paths$dives)) {
  defaults = compose_cfg(file = 'conf/config.yaml')
  cfg$sub_paths$dives = defaults$sub_paths$dives
}

o = file.path(out.dir, cfg$sub_paths$figures, cfg$sub_paths$dives)

dir.create(o, recursive = TRUE)

for(i in 1:length(dives.obs)) {
  d = dives.obs[[i]]
  pl = plot(x = d$dive, depth.bins = d$depth.bins, errorbars = TRUE, 
            time.as.POSIXct = TRUE)
  ggsave(pl, filename = file.path(o, paste('dive', i, '.png', sep ='')), 
         dpi = 'print')
}