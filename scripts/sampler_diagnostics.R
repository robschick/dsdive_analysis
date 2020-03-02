# configuration tools
library(composr)
library(yaml)
# modeling tools
library(dsdive, lib.loc = c('.', .libPaths()))
# mcmc tools
library(coda)
# plotting tools
library(ggplot2)
library(ggthemes)


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

groups = list(validation = 'holdout_half', sampler = 'prod')

# build configuration
cfg = compose_cfg(file = file.path('conf', 'config.yaml'), groups = groups)
rm(groups)

# output paths
out.dir = file.path(cfg$base_paths$fit, cfg$data$name, cfg$subset$name, 
                    cfg$validation$name, cfg$priors$name)

# load data and utility functions
source(file.path('scripts', 'utils', 'datafns.R'))


#
# load posterior samples and indices
#

load(file.path(out.dir, cfg$base_names$fit))
load(file.path(out.dir, cfg$base_names$fit_inds))


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
# plot dives
#

o = file.path(out.dir, cfg$sub_paths$figures, 'dives', 'training')
dir.create(o, recursive = TRUE)

for(i in 1:length(fit.inds$fit)) {
  pl = plot(x = dives.obs[[i]]$dive, depth.bins = dives.obs[[i]]$depth.bins, 
            errorbars = TRUE)
  ggsave(pl, filename = file.path(o, paste('dive', i, '.png', sep='')), 
         dpi = 'print')
}


#
# plot posteriors for trace variables
#


o = file.path(out.dir, cfg$sub_paths$figures, 'sampler_diagnostics')
dir.create(o, recursive = TRUE)

for(i in 1:ncol(state$trace.offset)) {
  png(file.path(o, paste('dive', i, 'offset.png', sep = '_')),
      width = 480*2, height = 480)
  plot(mcmc(state$trace.offset[-(1:cfg$sampler$burn),i]))
  dev.off()
}


#
# plot posteriors for model parameters
#

state$trace = state$trace[-nrow(state$trace),]

for(i in 1:ncol(state$trace)) {
  png(file.path(o, paste(colnames(state$trace)[i], '.png', sep = '')),
      width = 480*2, height = 480)
  plot(mcmc(state$trace[-(1:cfg$sampler$burn),i]))
  dev.off()
}
