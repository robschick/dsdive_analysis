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


#
# Set up environment
#

# groups.compare = list(
#   exact = list(validation="holdout_half",
#        observation_model="exact_systematic",
#        sampler="prod"),
#   uniform = list(validation="holdout_half",
#        observation_model="uniform_systematic",
#        sampler="prod")
# )

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
  )
)

# build configuration
cfg.list = lapply(groups.compare, function(g) {
  compose_cfg(file = file.path('conf', 'config.yaml'), groups = g)
})
rm(groups.compare)

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
    LoadToEnvironment(file.path(o, cfg$sub_paths$figures, 
                                cfg$sub_paths$posterior_predictions, 
                                filename))
  })
}


#
# compare depth-by-time validations
#

valcmp = LoadValidationData('depths_by_time.RData')

names.group = CapStr(names(cfg.list))

df = rbind(
  valcmp[[2]]$df %>% ungroup() %>% mutate(T0 = 'Validation') %>% 
    filter(Distribution != 'Post. Predictive'),
  valcmp[[1]]$df %>% ungroup() %>% mutate(T0 = names.group[1]) %>% 
    filter(Distribution == 'Post. Predictive') %>% 
    mutate(cdf.lwr = NA, cdf.upr = NA),
  valcmp[[2]]$df %>% ungroup() %>% mutate(T0 = names.group[2]) %>% 
    filter(Distribution == 'Post. Predictive') %>% 
    mutate(cdf.lwr = NA, cdf.upr = NA)
) %>% mutate(T0 = factor(T0))

pl = ggplot(df, aes(x = (depth), y = cdf, ymin = cdf.lwr, ymax = cdf.upr,
                    fill = T0, col = T0, group = T0)) + 
  geom_ribbon(alpha = .05, col = NA) + 
  geom_point() + 
  geom_line(lty = 3, alpha = .6) + 
  scale_color_brewer(type = 'qual', palette = 'Dark2') +
  scale_fill_brewer(type = 'qual', palette = 'Dark2') +
  xlab('Depth (m)') + 
  ylab('CDF') + 
  facet_wrap(~factor(time/60)) +  
  theme_few()

# save comparison plot to each of the cfg directories
sapply(1:length(cfg.list), function(i) {
  cfg = cfg.list[[i]]
  o = out.dir.list[i]
  o2 = file.path(o, cfg$sub_paths$figures, cfg$sub_paths$posterior_predictions, 
                 cfg$sub_paths$comparisons)
  dir.create(o2, recursive = TRUE)
  ggsave(pl, filename = file.path(o2, 'depths_by_time_comparison.png'), 
         dpi = 'print', width = 14, height = 14)
})


#
# compare depth-by-time imse's
#

valcmp = LoadValidationData('depths_by_time.RData')

names.group = CapStr(names(cfg.list))

df = do.call(rbind, mapply(function(nom, val) {
  data.frame(do.call(rbind, lapply(val$r.imse, function(x) x)),
        Model = nom)
}, names.group, valcmp, SIMPLIFY = FALSE))

ggplot(df, aes(x = Model, y = imse)) + 
  geom_point() + 
  facet_wrap(~factor(s/60), scales = 'free') + 
  theme_few() + 
  theme(panel.border = element_blank())

# save comparison summaries to each of the cfg directories
sapply(1:length(cfg.list), function(i) {
  cfg = cfg.list[[i]]
  o = out.dir.list[i]
  o2 = file.path(o, cfg$sub_paths$figures, cfg$sub_paths$posterior_predictions, 
                 cfg$sub_paths$comparisons)
  dir.create(o2, recursive = TRUE)
  
  sink(file.path(o2, 'depths_by_time_imse_comparison.txt'))
  
  cat("==================================================\n")
  cat("Average integrated mean square error criterion\n")
  cat("==================================================\n")
  
  cat("\n")
  
  print(
    df %>%
      filter(s/60 <= 60) %>% 
      group_by(Model) %>% 
      summarise(mean.imse = mean(imse)) %>% 
      arrange(mean.imse)
  )
  
  sink()
  
})


#
# compare duration imse's
#

valcmp = LoadValidationData('observed_duration.RData')

names.group = CapStr(names(cfg.list))

df = do.call(rbind, mapply(function(nom, val) {
  data.frame(imse = val$r.imse, Model = nom)
}, names.group, valcmp, SIMPLIFY = FALSE))

# save comparison summaries to each of the cfg directories
sapply(1:length(cfg.list), function(i) {
  cfg = cfg.list[[i]]
  o = out.dir.list[i]
  o2 = file.path(o, cfg$sub_paths$figures, cfg$sub_paths$posterior_predictions, 
                 cfg$sub_paths$comparisons)
  dir.create(o2, recursive = TRUE)
  
  sink(file.path(o2, 'observed_duration_imse_comparison.txt'))
  
  cat("==================================================\n")
  cat("Integrated mean square error criterion\n")
  cat("==================================================\n")
  
  cat("\n")
  
  print(
    df %>%
      arrange(imse)
  )
  
  sink()
  
})


#
# compare stage duration imse's
#

valcmp = LoadValidationData('stage_duration.RData')

names.group = CapStr(names(cfg.list))

df = do.call(rbind, mapply(function(nom, val) {
  res = data.frame(do.call(rbind, lapply(val$r.imse, function(x) x)),
                   Model = nom)
  res$imse = as.numeric(levels(res$imse)[res$imse])
  res
}, names.group, valcmp, SIMPLIFY = FALSE))

# save comparison summaries to each of the cfg directories
sapply(1:length(cfg.list), function(i) {
  cfg = cfg.list[[i]]
  o = out.dir.list[i]
  o2 = file.path(o, cfg$sub_paths$figures, cfg$sub_paths$posterior_predictions, 
                 cfg$sub_paths$comparisons)
  dir.create(o2, recursive = TRUE)
  
  sink(file.path(o2, 'stage_duration_imse_comparison.txt'))
  
  cat("==================================================\n")
  cat("Integrated mean square error criterion\n")
  cat("==================================================\n")
  
  cat("\n")
  
  print(df %>% filter(s=='Stage 1') %>% arrange(imse)) 
  
  cat('\n')
  
  print(df %>% filter(s=='Stage 2') %>% arrange(imse)) 
  
  cat('\n')
  
  print(df %>% filter(s=='Stage 3') %>% arrange(imse)) 
  
  cat('\nOverall\n')
  
  print(df %>% 
          group_by(Model) %>% 
          summarise(mean.imse = mean(imse)) %>% 
          arrange(mean.imse))
  
  sink()
  
})


#
# compare max observed depth imse's
#

valcmp = LoadValidationData('max_observed_depth.RData')

names.group = CapStr(names(cfg.list))

df = do.call(rbind, mapply(function(nom, val) {
  data.frame(imse = val$r.imse, Model = nom)
}, names.group, valcmp, SIMPLIFY = FALSE))

# save comparison summaries to each of the cfg directories
sapply(1:length(cfg.list), function(i) {
  cfg = cfg.list[[i]]
  o = out.dir.list[i]
  o2 = file.path(o, cfg$sub_paths$figures, cfg$sub_paths$posterior_predictions, 
                 cfg$sub_paths$comparisons)
  dir.create(o2, recursive = TRUE)
  
  sink(file.path(o2, 'max_observed_depth_imse_comparison.txt'))
  
  cat("==================================================\n")
  cat("Integrated mean square error criterion\n")
  cat("==================================================\n")
  
  cat("\n")
  
  print(
    df %>% arrange(imse)
  )
  
  sink()
  
})
