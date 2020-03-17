# assess parameter recovery in simulation

library(coda)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(dsdive)
library(dplyr)
library(stringr)
library(scoringRules)
library(composr)


#
# set up environment
#

# load data and utility functions
source(file.path('scripts', 'utils', 'datafns.R'))
source(file.path('scripts', 'utils', 'LoadToEnvironment.R'))


#
# define simulation configurations to post-process
#

sim.series = 'tyack_more_known_end'
data.groups = gsub(pattern = '\\.yaml', replacement = '', 
                   x = dir(path = file.path('conf', 'data'), 
                           pattern = sim.series))

cfg.list = sweep_cfg(file = file.path('conf', 'config.yaml'), 
                     groups = list(
                       data = data.groups,
                       observation_model = 'exact_systematic',
                       priors = 'tyack_simulation_priors',
                       sampler = 'prod'))

rm(sim.series, data.groups)


#
# load data and model parameters
#

cfg = cfg.list[[1]]

dives.obs = dives.load(path = cfg$data$path, 
                       dive_pattern = cfg$data$file_patterns$dive,
                       depth_pattern = cfg$data$file_patterns$depths)

dives.truth = dives.load(path = file.path(cfg$data$path, '..', 'truth'), 
                         dive_pattern = cfg$data$file_patterns$dive,
                         depth_pattern = cfg$data$file_patterns$depths)

load(file.path(cfg$data$path, '..', 'params', 'params.RData'))

rm(cfg)

# add stage transition times to true dives
dives.truth = lapply(dives.truth, function(d) {
  d$t.stages = d$dive$times[c(FALSE, diff(d$dive$stages)==1)]
  d
})

# determine output directories, and build plot directories
out.dir.list = lapply(cfg.list, function(cfg) {
  out.dir = file.path(cfg$base_paths$fit, cfg$data$name, cfg$subset$name, 
                      cfg$validation$name, cfg$observation_model$name, 
                      cfg$priors$name)
  dir.create(file.path(out.dir, cfg$sub_paths$figure, 
                       cfg$sub_paths$comparisons), recursive = TRUE)
  out.dir
})


#
# load and preprocess mcmc outputs
#

# load mcmc outputs
trace.theta = lapply(1:length(cfg.list), function(i) {
  cfg = cfg.list[[i]]
  out.dir = out.dir.list[[i]]
  # load mcmc samples
  r = LoadToEnvironment(file.path(out.dir, cfg$base_names$fit))
  # add timestep and burnin period
  r$tstep = cfg$data$tstep
  cfg$sampler$burn = 250
  r$burn = cfg$sampler$burn
  # return
  r
})

# compute posterior summaries
trace.summaries = lapply(trace.theta, function(trace) {
  if(nrow(trace$state$theta) > trace$burn) {
    burn = 1:trace$burn
    o = mcmc(trace$state$theta[-burn,])
    data.frame(mean = colMeans(o), HPDinterval(o), 
               truth = unlist(params)[1:5],
               param = colnames(o), 
               tstep = trace$tstep)
  } else {
    NULL
  }
})

# arrange posterior samples for t.stages by dive.ind
trace.t.stages = lapply(trace.theta, function(trace) {
  n = length(trace$state$trace.t.stages[[1]])
  lapply(1:n, function(dive.ind) {
    t(sapply(trace$state$trace.t.stages, function(t.stages.family) {
      t.stages.family[[dive.ind]]
    }))
  })
})

# posterior summaries of t.stages recovery
trace.t.stages.summaries = lapply(trace.t.stages, function(trace) {
  n = length(trace)
  if(length(trace[[1]]) > 1e1) {
    burn = 1:500
    do.call(rbind, lapply(1:n, function(dive.ind) {
      # extract true stage transition times
      t.stages.true = dives.truth[[dive.ind]]$t.stages
      # convert to stage durations
      dur.stages.true = c(t.stages.true[1], diff(t.stages.true))
      # extract posterior samples for stage transition ties
      m = mcmc(trace[[dive.ind]][-burn,])
      # extract posterior samples for stage durations
      m.dur = mcmc(t(apply(m, 1, function(r) c(r[1], diff(r)))))
      # posterior summaries of stage transition times
      df = data.frame(post.mean = colMeans(m), 
                      HPDinterval(m),
                      truth = t.stages.true,
                      variable = c('T1', 'T2'),
                      dive.ind = dive.ind,
                      crps = crps_sample(t.stages.true/60, t(m)/60)) %>% 
        mutate(covered = lower <= truth & truth <= upper)
      # posterior summaries of stage durations
      df2 = data.frame(post.mean = colMeans(m.dur), 
                       HPDinterval(m.dur),
                       truth = dur.stages.true,
                       variable = c('T1.dur', 'T2.dur'),
                       dive.ind = dive.ind,
                       crps = crps_sample(dur.stages.true/60, t(m.dur)/60)) %>% 
        mutate(covered = lower <= truth & truth <= upper)
      # return both sets of posterior summaries
      rbind(df, df2)
    }))
  } else {
    NULL
  }
})

# plots to demonstrate recovery of t.stages for each dive.
#   we see that recovery of the t.stage values becomes more difficult as the 
#   parameters themselves become more extreme.
t.stages.pl = lapply(1:length(cfg.list), function(i) {
  df = trace.t.stages.summaries[[i]]
  cfg = cfg.list[[i]]
  # build plot
  pl = ggplot(df, aes(x = truth/60, y = post.mean/60, ymin = lower/60, 
                      ymax = upper/60)) + 
    geom_abline(slope = 1, intercept = 0) +
    geom_pointrange() + 
    facet_wrap(~variable, scales = 'free') + 
    xlab('True stage transition time (min.)') + 
    ylab('Posterior estimate') + 
    theme_few() + 
    theme(panel.border = element_blank())
  # save plot
  ggsave(pl + ggtitle(paste(cfg$data$tstep, 'sec. between observations', 
                            sep = ' ')),
         filename = file.path(out.dir.list[[i]], cfg$sub_paths$figures,
                              'tstages.png'), 
         dpi = 'print', width = 15, height = 10)
  # return plot
  pl
})


# empirical coverage rate of posterior credible intervals for T1, T2 is good;
# crps scores show relatively small errors in stage duration recovery
lapply(1:length(cfg.list), function(i) {
  cfg = cfg.list[[i]]
  df = trace.t.stages.summaries[[i]]
  d = df %>% 
    group_by(variable) %>% 
    summarise(coverage = mean(covered), 
              crps = mean(crps))
  write.csv(d, file = file.path(out.dir.list[[i]], cfg$sub_paths$figures,
                                'tstages_coverage.csv'), row.names = FALSE)
})

# plot to demonstrate recovery of model parameter for varying choice of 
# model parameter, holding all others fixed
series.plot = function(param.name) {
  
  df = do.call(rbind, lapply(trace.summaries, function(trace.sum){ 
    if(!is.null(trace.sum)) {
      p.name = gsub('beta', 'pi', param.name)
      trace.sum %>% filter(param == p.name)
    } else {
      NULL
    }
  })) 
  
  df %>% ggplot(aes(x = factor(tstep), y = mean, ymin = lower, ymax = upper)) + 
    geom_pointrange() + 
    xlab('Time between observations (sec.)') + 
    ylab('Posterior estimate') + 
    theme_few() + 
    theme(panel.border = element_blank())
}


sapply(1:length(cfg.list), function(i) {
  cfg = cfg.list[[i]]
  ggsave(series.plot('lambda1') + 
           geom_hline(yintercept = params$lambda[1], lty = 3),
         filename = file.path(out.dir.list[[i]], cfg$sub_paths$figures, 
                              cfg$sub_paths$comparisons, 'lambda1_series.png'), 
         dpi = 'print')
})

sapply(1:length(cfg.list), function(i) {
  cfg = cfg.list[[i]]
  ggsave(series.plot('lambda2') + 
           geom_hline(yintercept = params$lambda[2], lty = 3),
         filename = file.path(out.dir.list[[i]], cfg$sub_paths$figures, 
                              cfg$sub_paths$comparisons, 'lambda2_series.png'), 
         dpi = 'print')
})

sapply(1:length(cfg.list), function(i) {
  cfg = cfg.list[[i]]
  ggsave(series.plot('lambda3') + 
           geom_hline(yintercept = params$lambda[3], lty = 3),
         filename = file.path(out.dir.list[[i]], cfg$sub_paths$figures, 
                              cfg$sub_paths$comparisons, 'lambda3_series.png'), 
         dpi = 'print')
})

sapply(1:length(cfg.list), function(i) {
  cfg = cfg.list[[i]]
  ggsave(series.plot('pi1') + 
           geom_hline(yintercept = params$beta[1], lty = 3),
         filename = file.path(out.dir.list[[i]], cfg$sub_paths$figures, 
                              cfg$sub_paths$comparisons, 'pi1_series.png'), 
         dpi = 'print')
})

sapply(1:length(cfg.list), function(i) {
  cfg = cfg.list[[i]]
  ggsave(series.plot('pi2') + 
           geom_hline(yintercept = params$beta[2], lty = 3),
         filename = file.path(out.dir.list[[i]], cfg$sub_paths$figures, 
                              cfg$sub_paths$comparisons, 'pi2_series.png'), 
         dpi = 'print')
})

sapply(trace.theta, function(trace) {
  nrow(trace$state$theta)
})

sapply(trace.theta, function(trace) {
  effectiveSize(mcmc(trace$state$theta))
})


#
# combined plot of all parameter recovery
#

# build plot
pl = do.call(rbind, trace.summaries) %>%
  # munge parameter labels
  mutate(param = recode_factor(param, pi1='pi^(1)', pi2='pi^(3)', 
                                 lambda1='lambda^(1)', lambda2='lambda^(2)',
                                 lambda3='lambda^(3)')) %>%
  # base plot
  ggplot(aes(x = factor(tstep), y = mean, ymin = lower, ymax = upper)) + 
  # parameter recovery
  geom_pointrange() + 
  # truth reference lines
  geom_hline(mapping = aes(yintercept = truth), lty = 3) + 
  # formatting
  facet_grid(param~., scales = 'free_y', 
             labeller = label_parsed, switch = 'both') + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 3.5)) +
  # labels
  xlab('Time between observations (sec.)') + 
  ylab('Posterior estimate') + 
  theme_few() + 
  theme(strip.text.y = element_text(angle = 180, face = 'bold', size = 12), 
        strip.placement = 'outside')

pl

# save plot
sapply(1:length(cfg.list), function(i) {
  cfg = cfg.list[[i]]
  sc = 1.75
  ggsave(pl,
         filename = file.path(out.dir.list[[i]], cfg$sub_paths$figures, 
                              cfg$sub_paths$comparisons, 'param_recovery.png'), 
         dpi = 'print', width = 3*sc, height = 4*sc)
})
