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
library(stringr)
library(ggpubr)


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

cfg = read_yaml(file = 'output/zc84/all_dives/no_validation/uniform_systematic/standard_priors/cfg.yaml')

# output paths
out.dir = file.path(cfg$base_paths$fit, cfg$data$name, cfg$subset$name, 
                    cfg$validation$name, cfg$observation_model$name, 
                    cfg$priors$name)

# load data and utility functions
source(file.path('scripts', 'utils', 'datafns.R'))
source(file.path('scripts', 'utils', 'LoadToEnvironment.R'))


#
# load data
#

dives.obs = dives.load(path = cfg$data$path, 
                       dive_pattern = cfg$data$file_patterns$dive,
                       depth_pattern = cfg$data$file_patterns$depths)

load(file.path(out.dir, cfg$base_names$fit_inds))


#
# plot imputations
#

if(is.null(cfg$sub_paths$imputations)) {
  defaults = compose_cfg(file = 'conf/config.yaml')
  cfg$sub_paths$imputations = defaults$sub_paths$imputations
}

o = file.path(out.dir, cfg$sub_paths$figures, cfg$sub_paths$imputations)

dir.create(o, recursive = TRUE)

for(dive.ind in 1:length(fit.inds$fit)) {
  
  # extract dive
  dive.id = fit.inds$fit[dive.ind]
  d = dives.obs[[dive.id]]$dive
  depth.bins = dives.obs[[dive.id]]$depth.bins
  
  # select imputations to plot
  impute.inds = round(seq(from = cfg$sampler$burn + 1, 
                          to = cfg$sampler$iterations, 
                          length.out = 100))
  
  # identify imputation files
  imputation.files = dir(path = file.path(out.dir, cfg$sub_paths$imputations, 
                                          paste('dive', dive.id, sep = '')), 
                         full.names = TRUE)
  
  # determine ranges of imputations
  imputation.ranges = sapply(
    str_extract_all(basename(imputation.files), '[0-9]+', simplify = FALSE), 
    function(inds) as.numeric(inds)[-1]
  )
  
  # re-sort file list in numerical (vs. alphabetical) order
  o.num = order(imputation.ranges[1,])
  imputation.files = imputation.files[o.num]
  imputation.ranges = imputation.ranges[,o.num]
  
  # determine which files have the requested imputations
  target.files = findInterval(impute.inds, imputation.ranges[1,])
  
  # extract imputations
  imputations = sapply(1:length(impute.inds), function(i) {
    # load imputation file
    load(imputation.files[target.files[i]])
    # determine offset for imputation
    offset = impute.inds[i] - imputation.ranges[1,target.files[i]] + 1
    # return extracted dive
    imputed[offset]
  })
  
  #
  # compute stage probabilities as an approximate function of time
  #
  
  # determine grid of times to investigate stage at
  t.win = range(d$times)
  t.locs = seq(from = t.win[1], to = t.win[2], length.out = 500) - t.win[1]
  
  stages.samples = sapply(imputations, function(d) {
    d$stages[findInterval(t.locs, d$times)]
  })
  
  df = rbind(
    data.frame(t = t.locs, p = rowMeans(stages.samples==1), Stage = 1),
    data.frame(t = t.locs, p = rowMeans(stages.samples==2), Stage = 2),
    data.frame(t = t.locs, p = rowMeans(stages.samples==3), Stage = 3)
  )
  
  df$t = as.POSIXct(df$t, origin = "1970-01-01", tz = "UTC")
  
  
  #
  # assemble plot
  #
  
  # plot imputed dive
  pl1 = plot(x = d, depth.bins = depth.bins, errorbars = TRUE, 
             imputed.list = imputations, imputed.alpha = .025, 
             time.as.POSIXct = TRUE) + 
    xlab('Dive time (hh:mm)') + 
    ggtitle(paste('Dive', dive.id, sep = ' ')) + 
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank())
  
  # plot stage transition probabilities
  pl2 = ggplot(df, aes(x=t, y = p, color = factor(Stage))) + 
    geom_line(alpha = .7, lwd = .7) + 
    scale_color_brewer(type = 'qual', palette = 'Dark2') + 
    xlim(as.POSIXct(t.win - t.win[1], origin = "1970-01-01", tz = "UTC")) + 
    xlab('Dive time (hh:mm)') + 
    ylab('P(Stage|Y)') + 
    guides(color = 'none') + 
    theme_few() + 
    theme(panel.border = element_blank())
  
  
  #
  # save plot
  #
  
  ggsave(file.path(o, paste('dive', dive.id, '.png', sep = '')), 
         ggarrange(pl1, pl2, nrow = 2, ncol = 1), 
         dpi = 'print', width = 8, height = 6)
}
