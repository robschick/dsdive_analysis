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

cfg = compose_cfg(file = 'output/zc84/all_dives/no_validation/uniform_systematic/standard_priors/cfg.yaml')

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


#
# compare priors and posteriors for stage transition times
#

# load parameters for stage transition priors
load(file.path(out.dir, cfg$base_names$stage_priors))

# MC approximation of stage 2->3 transition time (vs. stage 2 duration)
mc.it = 1e5
T2.prior.samples = rgamma(n = mc.it, shape = T1.prior.params[1], 
                          rate = T1.prior.params[2]) + 
  rgamma(n = mc.it, shape = T2.prior.params[1], rate = T2.prior.params[2])

# extract posterior samples of stage 1 transition times
T1.mcmc = sapply(state$trace.t.stages, function(t.stages.list) {
  sapply(t.stages.list, function(t.stages) { t.stages[1] })
})

# extract posterior samples of stage 2 transition times
T2.mcmc = sapply(state$trace.t.stages, function(t.stages.list) {
  sapply(t.stages.list, function(t.stages) { t.stages[2] })
})

# extract posterior samples of stage 1 durations
stage1.dur.mcmc = sapply(state$trace.t.stages[-(1:cfg$sampler$burn)], 
                         function(t.stages.list) {
  sapply(t.stages.list, function(t.stages) { t.stages[1] })
})

# extract posterior samples of stage 2 durations
stage2.dur.mcmc = sapply(state$trace.t.stages[-(1:cfg$sampler$burn)], 
                         function(t.stages.list) {
  sapply(t.stages.list, function(t.stages) { diff(t.stages) })
})

# compare prior and posterior for stage 1 transition times
png(file.path(o, 'T1_learning.png'))
plot(density(as.numeric(T1.mcmc[,-(1:cfg$sampler$burn)])/60),
     xlab = expression(T^(1)~'(min)'), 
     main = 'Posterior (black) vs. Prior (red)')
curve(60*dgamma(x*60, shape = T1.prior.params[1], rate = T1.prior.params[2]), 
      from = min(T1.mcmc)/60-1, to = max(T1.mcmc)/60+1, col = 2, add = TRUE)
dev.off()

# compare prior and posterior for stage 2 transition times
png(file.path(o, 'T2_learning.png'))
plot(density(as.numeric(T2.mcmc[,-(1:cfg$sampler$burn)])/60),
     xlab = expression(T^(2)~'(min)'), 
     main = 'Posterior (black) vs. Prior (red)')
lines(density(T2.prior.samples/60), col = 2)
dev.off()

# plot posterior for time in stage 1
png(file.path(o, 'stage1_dur.png'), width = 480*2)
plot(density(stage1.dur.mcmc/60), xlab = 'Descent duration (min)', 
     main = 'Posterior density')
dev.off()

# plot posterior for time in stage 2
png(file.path(o, 'stage2_dur.png'), width = 480*2)
plot(density(stage2.dur.mcmc/60), xlab = 'Foraging duration (min)', 
     main = 'Posterior density')
dev.off()