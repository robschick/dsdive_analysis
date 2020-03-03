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

groups.compare = list(
  exact = list(validation="holdout_half",
       observation_model="exact_systematic",
       sampler="prod"),
  uniform = list(validation="holdout_half",
       observation_model="uniform_systematic",
       sampler="prod")
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
# compare validations for depths_by_time
#

valcmp = LoadValidationData('depths_by_time.RData')

df = rbind(
  valcmp[[2]]$df %>% ungroup() %>% mutate(T0 = 'Validation') %>% 
    filter(Distribution != 'Post. Predictive'),
  valcmp[[1]]$df %>% ungroup() %>% mutate(T0 = 'Exact') %>% 
    filter(Distribution == 'Post. Predictive') %>% 
    mutate(cdf.lwr = NA, cdf.upr = NA),
  valcmp[[2]]$df %>% ungroup() %>% mutate(T0 = 'Uniform') %>% 
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
