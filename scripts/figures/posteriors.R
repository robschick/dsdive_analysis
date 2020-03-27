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
library(xtable)
library(tidyr)


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
  validation= 'no_validation'
)

# build configuration
cfg = compose_cfg(file = file.path('conf', 'config.yaml'), groups = groups)
rm(groups)

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
# load posterior predictive samples, if they exist
#

pred.dir = file.path(out.dir, cfg$sub_paths$posterior_predictions)
if(dir.exists(pred.dir)) {
  postpred.files = dir(file.path(out.dir, cfg$sub_paths$posterior_predictions), 
                       pattern = '*.RData', full.names = TRUE)
  postpred.samples = do.call(c, lapply(postpred.files, function(f) {
    load(f)
    samples
  }))
  
}


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

# publication-formatted summaries
sink(file.path(o, 'diving_parameters_latex.txt'))
df= data.frame(Mean = round(colMeans(m), 2),
              `Sd.` = round(apply(m,2,sd), 2),
              `HPD Interval` = apply(round(HPDinterval(m),2), 1, function(r) {
                paste('(', r[1], ', ', r[2], ')', sep ='')
              }))
rownames(df) = c('$\\pi^{(1)}$', '$\\pi^{(3)}$', '$\\lambda^{(1)}$', 
                 '$\\lambda^{(2)}$', '$\\lambda^{(3)}$')
print(xtable(df, caption = paste(
  'Posterior means, standard deviations, and 95\\% Highest posterior density', 
  '(HPD) intervals for model parameters $\\v\\Theta$ when fit to satellite tag',
  'data.', sep = ' '), 
  label = 'table:post_parameters'), 
  booktabs = TRUE, sanitize.text.function = function(x) x, timestamp = NULL, 
  comment = FALSE
)
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

# posterior correlations for model parameters
sink(file.path(o, 'parameter_correlations.txt'))
round(cor(state$theta),2)
sink()

# posterior t0 offests
png(file.path(o, 't0_offset.png'), width = 480*2)
plot(density(state$trace.offsets[-(1:cfg$sampler$burn),]), 
     xlab = expression(epsilon), main = 'Posterior density')
dev.off()

# posterior tf offests
png(file.path(o, 'tf_offset.png'), width = 480*2)
plot(density(state$trace.offsets.tf[-(1:cfg$sampler$burn),]), 
     xlab = expression(epsilon), main = 'Posterior density')
dev.off()


#
# processs posterior predictive distribution, if it exists
#

if(exists('postpred.samples')) {
  
  # set plotting directory
  o = file.path(out.dir, cfg$sub_paths$figures, 
                cfg$sub_paths$posterior_predictions)
  
  # create directory, if necessary
  dir.create(o, recursive = TRUE)
  
  # build non-overlapping depth bins, for display
  depth.bins = dives.obs[[1]]$depth.bins %>% 
    mutate(bin.min = center - halfwidth, 
           bin.max = center + halfwidth,
           bin.ind = 1:nrow(depth.bins))
  depth.bins$bin.min[-1] = depth.bins$bin.max[1:(nrow(depth.bins)-1)]
  
  # format depth bin display labels
  depth.bins = depth.bins %>% mutate(
    bin.range = paste('[', round(bin.min), ', ', round(bin.max), ')', 
                      sep =''),
    bin.range = ordered(bin.range, bin.range[bin.ind])
  )
  
  # only output information about posterior predictive dives that exceed this 
  # threshold
  max_depth = 800
  
  burn = 1:cfg$sampler$burn
  
  # extract summary information for posterior predictive dives
  postpred.summaries = do.call(rbind, lapply(postpred.samples[-burn], 
                                             function(d) {

    # compute stage durations
    t.stages = c(0,
                 d$times[c(FALSE, diff(d$stages)==1)],
                 d$times[length(d$times)])
    stage.durations = diff(t.stages)
    
    # max depth bin
    max.bin.ind = max(d$depths)
    
    # summary information
    data.frame(
      max.depth = depth.bins$center[max.bin.ind],
      bin.range = depth.bins$bin.range[max.bin.ind],
      duration = t.stages[3],
      'Stage 1' = stage.durations[1],
      'Stage 2' = stage.durations[2],
      'Stage 3' = stage.durations[3]
    )
  }))
  
  # dive durations
  pl = postpred.summaries %>% 
    filter(max.depth >= max_depth) %>% 
    ggplot(aes(x = duration/60)) + 
    stat_density(geom = 'line') + 
    xlab('Dive duration (min.)') + 
    ylab('Post. Pred. Density') + 
    theme_few() + 
    theme(panel.border = element_blank())
  
  ggsave(pl, filename = file.path(o, 'durations.png'), dpi = 'print', 
         width = 5, height = 4)
  
  # max depth
  pl = postpred.summaries %>% 
    filter(max.depth >= max_depth) %>%
    mutate(total = length(max.depth)) %>% 
    group_by(bin.range) %>%
    summarise(prob = n() / total[1]) %>%
    ggplot(aes(x = bin.range, y = prob)) + 
    geom_bar(stat = 'identity') +
    xlab('Maximum depth range (m)') + 
    ylab('Post. Pred. Probability') + 
    theme_few() + 
    theme(panel.border = element_blank())
  
  ggsave(pl, filename = file.path(o, 'max_depths.png'), dpi = 'print', 
         width = 5, height = 4)
  
  label_sub = function (labels, multi_line = TRUE) {
    labels <- lapply(labels, function(lab) {
      gsub('\\.', ' ', as.character(lab)) 
      })
    if (multi_line) {
      labels
    }
    else {
      collapse_labels_lines(labels)
    }
  }
  
  # stage durations
  pl = postpred.summaries %>% 
    filter(max.depth >= max_depth) %>% 
    pivot_longer(cols = c('Stage.1', 'Stage.2', 'Stage.3'), 
                 names_to = 'stage', 
                 values_to = 'stage_duration') %>% 
    ggplot(aes(x = stage_duration/60)) + 
    stat_density(geom = 'line') + 
    xlab('Stage duration (min)') + 
    ylab('Post. Pred. Density') + 
    facet_wrap(~stage, labeller = label_sub) + 
    theme_few()
  
  ggsave(pl, filename = file.path(o, 'stage_durations.png'), dpi = 'print', 
         width = 8, height = 4)
  
  #
  # textual summaries
  #
  
  sink(file.path(o, 'post_pred_summaries.txt'))
  
  summary(mcmc(postpred.summaries %>% 
                 filter(max.depth >= max_depth) %>% 
                 select(-bin.range))/60)
  
  round(HPDinterval(mcmc(postpred.summaries %>% 
                filter(max.depth >= max_depth) %>% 
                  select(-bin.range))/60))
  
  sink()
  
}