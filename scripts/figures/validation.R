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

groups=list(validation="holdout_half",observation_model="uniform_systematic",sampler="prod")


# build configuration
cfg = compose_cfg(file = file.path('conf', 'config.yaml'), groups = groups)
rm(groups)

cfg = read_yaml(file = 'output/zc84_bak/all_dives/holdout_half/exact_systematic/standard_priors/cfg.yaml')

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

# identify validation dives
load(file.path(out.dir, cfg$base_names$fit_inds))

# extract validation dives
depth.bins = dives.obs[[1]]$depth.bins
dives.obs.list = lapply(dives.obs, function(d) d$dive)[fit.inds$validate]

tstep = diff(dives.obs[[1]]$dive$times[1:2])
tstep.seq = tstep

# extract information about validation dives
n.dives = length(dives.obs.list)
validation.obs = do.call(rbind, lapply(1:n.dives, function(ind) {
  # extract dive 
  d = dives.obs.list[[ind]]
  # extract estimates of stage durations
  stages.dur = diff(c(0, d$times[c(FALSE, diff(d$stages.est)==1)], 
                      d$times[length(d$times)]))
  # loop over all observation times
  do.call(rbind, lapply(tstep.seq, function(tstep) {
    # observed duration of dive
    duration.obs = d$times[length(d$times)]
    # extract summary features of dives
    data.frame(
      dive.id = ind,
      duration.obs = d$times[length(d$times)],
      max.depth.obs = depth.bins$center[max(d$depths)],
      n.obs = length(d$times),
      n.tx.obs = sum(diff(d$depths) != 0),
      tstep = tstep,
      dur.s1 = stages.dur[1],
      dur.s2 = stages.dur[2], 
      dur.s3 = stages.dur[3]
    )
  }))
}))

# compute depth bin distributions over time
validation.depthseq = do.call(rbind, lapply(1:n.dives, function(ind) {
  # extract dive
  d = dives.obs.list[[ind]]
  # return tidy data frame with depth bins over time
  data.frame(dive.id = ind, depth = depth.bins$center[d$depths], time = d$times)
}))


#
# load posterior predictive samples
#

postpred.files = dir(file.path(out.dir, cfg$sub_paths$posterior_predictions), 
                     pattern = '*.RData', full.names = TRUE)
postpred.samples = do.call(c, lapply(postpred.files, function(f) {
  load(f)
  samples
}))


#
# observe posterior predictive samples
#

# extract information about posterior predictive samples
n.samples = length(postpred.samples)
postpred.samples.obs = do.call(rbind, lapply(1:n.samples, function(ind) {
  # extract dive 
  d = postpred.samples[[ind]]
  # loop over all observation times
  do.call(rbind, lapply(tstep.seq, function(tstep) {
    # exact duration of dive
    duration = d$times[length(d$times)]
    # build sequence of observation times
    t.obs.reported = seq(from = 0, to = duration + 2*tstep, by = tstep)
    # shift observation times according to offset
    t.obs = t.obs.reported - d$t0.offset
    # observe dive
    obs = dsdive.observe(depths = d$depths, times = d$times, 
                         stages = d$stages, t.obs = t.obs)
    # correct observations for offset
    obs$depths = c(rep(1, sum(t.obs<0)), obs$depths)
    obs$stages = c(rep(1, sum(t.obs<0)), obs$stages)
    obs$times = t.obs.reported
    # remove trailing surface observations
    trailing.obs = which(obs$depths==1)[-(1:2)]
    if(length(trailing.obs) > 0) {
      obs$depths = obs$depths[-trailing.obs]
      obs$stages = obs$stages[-trailing.obs]
      obs$times = obs$times[-trailing.obs]
    }
    # compute observed stage durations
    stages.dur = diff(c(0, obs$times[c(FALSE, diff(obs$stages)==1)], 
                        obs$times[length(obs$times)]))
    # extract summary features of dives
    data.frame(
      sample = ind,
      duration = duration,
      duration.obs = t.obs.reported[length(t.obs.reported)],
      depth.start.obs = obs$depths[1],
      max.depth = depth.bins$center[max(d$depths)],
      max.depth.obs = depth.bins$center[max(obs$depths)],
      n.obs = length(t.obs),
      n.tx = sum(diff(d$depths) != 0),
      n.tx.obs = sum(diff(obs$depths) != 0),
      tstep = tstep,
      dur.s1 = stages.dur[1],
      dur.s2 = stages.dur[2], 
      dur.s3 = stages.dur[3]
    )
  }))
}))

# compute depth bin distributions over time
postpred.depthseq = do.call(rbind, lapply(1:n.samples, function(ind) {
  # extract dive 
  d = postpred.samples[[ind]]
  # exact duration of dive
  duration = d$times[length(d$times)]
  # build sequence of observation times
  t.obs.reported = seq(from = 0, to = duration + 2*tstep, by = tstep)
  # shift observation times according to offset
  t.obs = t.obs.reported - d$t0.offset
  # observe dive
  obs = dsdive.observe(depths = d$depths, times = d$times, 
                       stages = d$stages, t.obs = t.obs)
  # correct observations for offset
  obs$depths = c(rep(1, sum(t.obs<0)), obs$depths)
  obs$stages = c(rep(1, sum(t.obs<0)), obs$stages)
  obs$times = t.obs.reported
  # remove trailing surface observations
  trailing.obs = which(obs$depths==1)[-(1:2)]
  if(length(trailing.obs) > 0) {
    obs$depths = obs$depths[-trailing.obs]
    obs$stages = obs$stages[-trailing.obs]
    obs$times = obs$times[-trailing.obs]
  }
  # return tidy data frame with depth bins over time
  data.frame(sample = ind, 
             depth = depth.bins$center[obs$depths], 
             duration.obs = t.obs.reported[length(t.obs.reported)],
             depth.start.obs = obs$depths[1],
             time = obs$times,
             max.depth.obs = depth.bins$center[max(obs$depths)])
}))



#
# discrete KS gof test formatting
#

ks.gof = function(df, var, verbose = TRUE) {
  
  # extract estimate of posterior predictive distribution
  prob = df %>% filter(series=='Post. Predictive') %>% ungroup()
  prob = prob[order(prob[var] %>% unlist()),]
  
  # extract validation samples
  obs = df %>% filter(series=='Empirical Validation') %>% ungroup() 
  
  # compute support, and merge data; replace missing items with 0
  df.ks = data.frame(support = sort(unique(c(prob[var] %>% unlist(), 
                                                obs[var] %>% unlist())))) %>% 
    left_join(prob %>% select(var, prob), by = c('support' = var)) %>% 
    left_join(obs %>% select(var, prob), by = c('support' = var)) %>% 
    mutate_all(~replace(., which(is.na(.)), 0))
  
  colnames(df.ks)[2:3] = c('prob.posterior', 'prob.observed')
  
  # expand validation samples to raw counts
  samples = rep(obs[var] %>% unlist(), obs$n)
  
  # run ks test
  res = disc_ks_test(x = samples, 
                     y = stepfun(x = df.ks$support, 
                                 y = c(0,cumsum(df.ks$prob.posterior))),
                     exact = TRUE)
  
  if(verbose) {

    cat("==================================\n")
    cat("K-S goodness of fit test\n")
    cat("==================================\n")
    
    cat("\n")
    
    print(c(res$statistic, p=res$p.value))
    
    cat("\n")
    
    print(df.ks)
    
    cat("\n")
    
    print(c(n=length(samples)))
  }
  
  res
}


#
# chisq gof test formatting
#

chisq.gof = function(df, var, verbose = TRUE, collapse = TRUE) {
  # Parameters:
  #  collapse - TRUE to merge cells with expected counts less than 5.  The 
  #    merging will be done working from the tails toward the center of the 
  #    distribution's support.
  
  # extract estimate of posterior predictive distribution
  prob = df %>% filter(series=='Post. Predictive') %>% ungroup() 
  
  # extract validation samples
  obs = df %>% filter(series=='Empirical Validation') %>% ungroup() 
  
  # compute support, and merge data; replace missing items with 0
  df.chisq = data.frame(support = sort(unique(c(prob[var] %>% unlist(), 
                                                obs[var] %>% unlist())))) %>% 
    left_join(prob %>% select(var, prob), by = c('support' = var)) %>% 
    left_join(obs %>% select(var, n), by = c('support' = var)) %>% 
    mutate_all(~replace(., which(is.na(.)), 0))
  
  colnames(df.chisq)[3] = 'observed'
  
  # run chisq test
  res = chisq.test(x = df.chisq$observed, p = df.chisq$prob)
  
  if(verbose) {
    df.chisq$expected = round(res$expected)
    
    cat("==================================\n")
    cat("Chi-squared goodness of fit test\n")
    cat("==================================\n")
    
    cat("\n")
    
    print(c(res$statistic, res$parameter, p=res$p.value))
    
    cat("\n")
    
    print(df.chisq)
  }
  
  res
}


#
# plots and information
#

o = file.path(out.dir, cfg$sub_paths$figures, 
              cfg$sub_paths$posterior_predictions)

dir.create(o, recursive = TRUE)

burn = 1:cfg$sampler$burn

# DKW CDF bounds, Source:
# https://en.wikipedia.org/wiki/CDF-based_nonparametric_confidence_interval

min_dur = cfg$subset$duration_min
max_dur = cfg$subset$duration_max
bin_start_max = cfg$subset$bin_start_max


#
# max observed depth distributions
#

df = rbind(
  postpred.samples.obs %>% 
    filter(!(sample %in% burn),
           max.depth.obs >= 1e3,
           duration.obs >= min_dur,
           duration.obs <= max_dur,
           depth.start.obs <= bin_start_max) %>% 
    mutate(series = 'Post. Predictive',
           total = length(unique(sample)),
           eps = sqrt(log(2/.05)/(2*total))) %>% 
    group_by(max.depth.obs, series) %>% 
    summarise(prob = n() / total[1],
              n = n(),
              eps = eps[1]) %>% 
    group_by(series) %>% 
    mutate(cdf = cumsum(prob),
           cdf.lwr = pmax(cdf - eps, 0),
           cdf.upr = pmin(cdf + eps, 1)),
  validation.obs %>% 
    mutate(series = 'Empirical Validation', 
           total = length(unique(dive.id)),
           eps = sqrt(log(2/.05)/(2*total))) %>% 
    group_by(max.depth.obs, series) %>% 
    summarise(prob = n() / total[1],
              n = n(),
              eps = eps[1]) %>% 
    group_by(series) %>% 
    mutate(cdf = cumsum(prob),
           cdf.lwr = pmax(cdf - eps, 0),
           cdf.upr = pmin(cdf + eps, 1))
)

pl = ggplot(df, aes(x = max.depth.obs, y = cdf, ymin = cdf.lwr, ymax = cdf.upr,
                    fill = series, col = series)) + 
  geom_ribbon(alpha = .05, col = NA) + 
  geom_point() + 
  geom_line(lty = 3, alpha = .6) + 
  scale_color_brewer('Distribution', type = 'qual', palette = 'Dark2') + 
  scale_fill_brewer('Distribution', type = 'qual', palette = 'Dark2') + 
  xlab('Max. observed depth (m)') + 
  ylab('CDF') + 
  theme_few() + 
  theme(panel.border = element_blank())

ggsave(pl, filename = file.path(o, 'max_observed_depth_cdf.png'), 
       dpi = 'print')

sink(file.path(o, paste('max_observed_depth_chisq.txt')))
r = chisq.gof(df, 'max.depth.obs')
sink()

sink(file.path(o, paste('max_observed_depth_ks.txt')))
r = ks.gof(df, 'max.depth.obs')
sink()

save(df, r, file = file.path(o, paste('max_observed_depth.RData')))

#
# observed duration distributions
#

df = rbind(
  postpred.samples.obs %>% 
    filter(!(sample %in% burn),
           max.depth.obs >= 1e3,
           duration.obs >= min_dur,
           duration.obs <= max_dur,
           depth.start.obs <= bin_start_max) %>%
    mutate(series = 'Post. Predictive',
           total = length(unique(sample)),
           eps = sqrt(log(2/.05)/(2*total))) %>% 
    group_by(duration.obs, series) %>% 
    summarise(prob = n() / total[1],
              eps = eps[1]) %>% 
    group_by(series) %>% 
    mutate(cdf = cumsum(prob),
           cdf.lwr = pmax(cdf - eps, 0),
           cdf.upr = pmin(cdf + eps, 1)),
  validation.obs %>% 
    mutate(series = 'Empirical Validation', 
           total = length(unique(dive.id)),
           eps = sqrt(log(2/.05)/(2*total))) %>% 
    group_by(duration.obs, series) %>% 
    summarise(prob = n() / total[1],
              n = n(),
              eps = eps[1]) %>% 
    group_by(series) %>% 
    mutate(cdf = cumsum(prob),
           cdf.lwr = pmax(cdf - eps, 0),
           cdf.upr = pmin(cdf + eps, 1))
)

pl = ggplot(df, aes(x = duration.obs/60, y = cdf, col = series, fill = series,
                    ymin = cdf.lwr, ymax = cdf.upr)) + 
  geom_ribbon(alpha = .05, col = NA) + 
  geom_point() + 
  geom_line(lty = 3, alpha = .6) + 
  scale_color_brewer('Distribution', type = 'qual', palette = 'Dark2') + 
  scale_fill_brewer('Distribution', type = 'qual', palette = 'Dark2') + 
  xlab('Observed dive duration (min)') + 
  ylab('CDF') + 
  theme_few() + 
  theme(panel.border = element_blank())

ggsave(pl, filename = file.path(o, 'observed_duration_cdf.png'), 
       dpi = 'print')

sink(file.path(o, paste('observed_duration_chisq.txt')))
r = chisq.gof(df, 'duration.obs')
sink()

sink(file.path(o, paste('max_observed_duration_ks.txt')))
r = ks.gof(df, 'duration.obs')
sink()

save(df, r, file = file.path(o, paste('observed_duration.RData')))


#
#  distribution of depths at observation times
#

df = rbind(
  postpred.depthseq %>% 
    filter(!(sample %in% burn),
           max.depth.obs >= 1e3,
           duration.obs >= min_dur,
           duration.obs <= max_dur,
           depth.start.obs <= bin_start_max) %>%
    group_by(time) %>% 
    mutate(total = length(unique(sample)),
           eps = sqrt(log(2/.05)/(2*total)),
           Distribution = 'Post. Predictive') %>% 
    group_by(time, depth, Distribution) %>% 
    summarise(prob = n() / total[1],
              n = n(),
              eps = eps[1]) %>%
    group_by(Distribution, time) %>% 
    mutate(cdf = cumsum(prob),
           cdf.lwr = pmax(cdf - eps, 0),
           cdf.upr = pmin(cdf + eps, 1)),
  validation.depthseq %>% 
    group_by(time) %>% 
    mutate(total = length(unique(dive.id)),
           eps = sqrt(log(2/.05)/(2*total)),
           Distribution = 'Empirical Validation') %>%
    group_by(time, depth, Distribution) %>% 
    summarise(prob = n() / total[1],
              n = n(),
              eps = eps[1]) %>%
    group_by(Distribution, time) %>% 
    mutate(cdf = cumsum(prob),
           cdf.lwr = pmax(cdf - eps, 0),
           cdf.upr = pmin(cdf + eps, 1))
)

pl = ggplot(df, aes(x = (depth), y = cdf, ymin = cdf.lwr, ymax = cdf.upr,
                    fill = Distribution, col = Distribution, 
                    group = Distribution)) + 
  geom_ribbon(alpha = .05, col = NA) + 
  geom_point() + 
  geom_line(lty = 3, alpha = .6) + 
  scale_color_brewer(type = 'qual', palette = 'Dark2') + 
  scale_fill_brewer(type = 'qual', palette = 'Dark2') +
  xlab('Depth (m)') + 
  ylab('CDF') + 
  facet_wrap(~factor(time/60)) +  
                     # labels = paste(unique(time)/60, 'min', sep = ' '))) + 
  theme_few()

ggsave(pl, filename = file.path(o, 'depths_by_time_cdf.png'), 
       dpi = 'print', width = 14, height = 14)

sink(file.path(o, paste('depths_by_time_chisq.txt')))
r = lapply(sort(unique(df$time))[-1], function(s) {
  res = chisq.gof(df %>% mutate(series = Distribution) %>% 
                    filter(time==s), 'depth')
  cat(paste("\n(t=", s/60, " min. results)\n\n\n", sep=''))
  res
})
sink()

sink(file.path(o, paste('depths_by_time_ks.txt')))
r = lapply(sort(unique(df$time))[-1], function(s) {
  res = ks.gof(df %>% mutate(series = Distribution) %>% 
                 filter(time==s), 'depth')
  cat(paste("\n(t=", s/60, " min. results)\n\n\n", sep=''))
  res
})
sink()

save(df, r, file = file.path(o, paste('depths_by_time.RData')))


#
# observed stage duration distributions
#

df = rbind(
  postpred.samples.obs %>% 
    filter(!(sample %in% burn),
           max.depth.obs >= 1e3,
           duration.obs >= min_dur,
           duration.obs <= max_dur,
           depth.start.obs <= bin_start_max) %>%
    mutate(series = 'Post. Predictive',
           total = length(unique(sample)),
           eps = sqrt(log(2/.05)/(2*total))) %>% 
    pivot_longer(cols = starts_with("dur."), 
                 names_to = 'stage', values_to = 'stage.duration') %>%
    filter(!is.na(stage.duration)) %>%
    mutate(stage = fct_recode(stage, 
                              'Stage 1' = 'dur.s1',
                              'Stage 2' = 'dur.s2',
                              'Stage 3' = 'dur.s3')) %>% 
    group_by(stage, stage.duration, series) %>% 
    summarise(prob = n() / total[1],
              n = n(),
              eps = eps[1]) %>% 
    group_by(stage, series) %>%
    mutate(cdf = cumsum(prob),
           cdf.lwr = pmax(cdf - eps, 0),
           cdf.upr = pmin(cdf + eps, 1)),
  validation.obs %>% 
    mutate(series = 'Empirical Validation', 
           total = length(unique(dive.id)),
           eps = sqrt(log(2/.05)/(2*total))) %>% 
    pivot_longer(cols = starts_with("dur."), 
                 names_to = 'stage', values_to = 'stage.duration') %>%
    mutate(stage = fct_recode(stage, 
                              'Stage 1' = 'dur.s1',
                              'Stage 2' = 'dur.s2',
                              'Stage 3' = 'dur.s3')) %>% 
    group_by(stage, stage.duration, series) %>% 
    summarise(prob = n() / total[1],
              n = n(),
              eps = eps[1]) %>% 
    group_by(stage, series) %>%
    mutate(cdf = cumsum(prob),
           cdf.lwr = pmax(cdf - eps, 0),
           cdf.upr = pmin(cdf + eps, 1))
) 


pl = ggplot(df, aes(x = stage.duration/60, y = cdf, ymin = cdf.lwr, 
                    ymax = cdf.upr, col = series, fill = series)) + 
  geom_ribbon(alpha = .05, col = NA) + 
  geom_point() + 
  geom_line(lty = 3, alpha = .6) + 
  scale_color_brewer('Distribution', type = 'qual', palette = 'Dark2') + 
  scale_fill_brewer('Distribution', type = 'qual', palette = 'Dark2') + 
  xlab('Observed stage duration (min)') + 
  facet_wrap(~stage) + 
  ylab('CDF') + 
  theme_few() 

ggsave(pl, filename = file.path(o, 'stage_duration_cdfs.png'), 
       dpi = 'print')

sink(file.path(o, paste('stage_duration_chisq.txt')))
r = lapply(levels(df$stage), function(s) {
  df2 = df %>% filter(stage==s) %>% group_by(series) %>% 
    mutate(prob = prob/sum(prob)) %>% ungroup()
  res = chisq.gof(df2, 'stage.duration')
  cat(paste("\n(", s," results)\n\n\n", sep=''))
  res
})
sink()

sink(file.path(o, paste('stage_duration_ks.txt')))
r = lapply(levels(df$stage), function(s) {
  df2 = df %>% filter(stage==s) %>% group_by(series) %>% 
    mutate(prob = prob/sum(prob)) %>% ungroup()
  res = ks.gof(df2, 'stage.duration')
  cat(paste("\n(", s," results)\n\n\n", sep=''))
  res
})
sink()

save(df, r, file = file.path(o, paste('stage_duration.RData')))

