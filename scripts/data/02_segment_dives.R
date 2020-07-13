library(anytime, lib.loc = c('singularity/libs', .libPaths()))
library(parallel, lib.loc = c('singularity/libs', .libPaths()))
library(stringr, lib.loc = c('singularity/libs', .libPaths()))
library(dplyr, lib.loc = c('singularity/libs', .libPaths()))
library(tidyr, lib.loc = c('singularity/libs', .libPaths()))


#
# output directories
#

out.dir = file.path('data', 'tag_labels')
dir.create(out.dir, recursive = TRUE)

diagnostic.dir = file.path(out.dir, 'diagnostics')
dir.create(diagnostic.dir, recursive = TRUE)


#
# set up workspace
#

# segmentation script
source(file.path('scripts', 'data', 'segment_fn.R'))

# tag plotting script
source(file.path('scripts', 'data', 'plot_tag.R'))

# load template bins
depth.bins = read.csv(file.path('data', 'imputed_bins', 'template.csv'))

# merge ratios
merge.seq = seq(from = .005, to = .99, by = .01)


#
# set up parallellization
#

cl = makeCluster(spec = detectCores(), type = 'SOCK')

clusterEvalQ(cl, library(anytime))
clusterEvalQ(cl, library(stringr))
clusterExport(cl, c('dive.segmentation', 'depth.bins', 'merge.seq'))


#
# find and process raw data
#

# locate raw sattag records
tag.files = dir(path = file.path('data', 'raw'), pattern = 'series_', 
                  full.names = TRUE)

# process raw files
tag.records = clusterApply(cl, tag.files, function(f) {
  
  # load data
  d = read.csv(file = f)
  d$Date = anytime(paste(d$Day, d$Time))
  
  # extract tag name
  tag.name = str_extract(f, 'Zc[0-9A-Za-z]+')

  # map all depths to standardized bins
  d$depth.bin = sapply(d$Depth, function(depth) {
    which.min(abs(depth - depth.bins$center))
  })
  d$depth.standardized = depth.bins$center[d$depth.bin]

  # label dives using range of ratios
  label.seq = lapply(merge.seq, function(mr) {
    list(
      labels = dive.segmentation(y = d$depth.bin, merge.ratio = mr,
                                 depth.bins = depth.bins),
      merge.ratio = mr
    )
  })

  # package results
  list(
    depths = d,
    label.seq = label.seq,
    name = tag.name
  )
})


stopCluster(cl)


#
# segmentation diagnostics
#

# extract summary features from records
labels.diagnostics = do.call(rbind, lapply(tag.records, function(record) {
  # loop over segmentations within record
  do.call(rbind, lapply(record$label.seq, function(segmentation) {
    # join depth data with segmentation labels and configuration
    cbind(record$depths, 
          mode = segmentation$labels, 
          merge.ratio = segmentation$merge.ratio) %>% 
      # restrict summary to known diving behavior
      filter(mode > 0) %>%
      # compute dive-level summaries
      group_by(mode, merge.ratio) %>% 
      summarise(maxd = max(depth.standardized),
                start = min(Date),
                end = max(Date),
                duration = difftime(time1 = end, time2 = start, units = 'mins'),
                deep = (maxd >= 800)) %>% 
      # aggregate summaries by dive type
      group_by(deep) %>% 
      summarise(n = length(deep), 
                duration.q25 = quantile(duration, probs = .25),
                duration.q5 = quantile(duration, probs = .5),
                duration.mean = mean(duration),
                duration.q75 = quantile(duration, probs = .75),
                merge.ratio = merge.ratio[1]) %>% 
      # add tag label and format for plotting
      mutate(tag = record$name,
             deep = factor(deep, 
                           levels = c(TRUE, FALSE), 
                           labels = c('Deep',' Shallow'))
             )
  }))
}))


# compare number of dives wrt. merge ratios
pl = ggplot(labels.diagnostics %>% 
              group_by(merge.ratio, deep) %>%
              summarise(n.mean = mean(n),
                        n.lwr = quantile(n, probs = .25),
                        n.upr = quantile(n, probs = .75)), 
            aes(x = merge.ratio, y = n.mean)) + 
  geom_line() + 
  geom_line(mapping = aes(x = merge.ratio, y = n.lwr), lty = 3, 
            inherit.aes = FALSE) + 
  geom_line(mapping = aes(x = merge.ratio, y = n.upr), lty = 3, 
            inherit.aes = FALSE) + 
  xlab('Merge ratio') + 
  ylab('Dives per record') + 
  facet_wrap(~deep, ncol = 1, scales = 'free') + 
  theme_few() + 
  theme(panel.background = element_blank(), strip.placement = 'inside') + 
  ggtitle('(Solid: mean across tags; Dotted: 25% and 75% quantiles)')

ggsave(pl, filename = file.path(diagnostic.dir, 'dives_per_ratio.pdf'), 
       type = 'cairo')



# compare dive durations wrt. merge ratios
pl = ggplot(labels.diagnostics %>% 
              group_by(merge.ratio, deep) %>% 
              summarise(duration.mean = mean(duration.mean),
                        duration.q25 = quantile(duration.q25, probs = .25),
                        duration.q75 = quantile(duration.q75, probs = .75)), 
            aes(x = merge.ratio, y = duration.mean)) + 
  geom_line() + 
  geom_line(mapping = aes(x= merge.ratio, y = duration.q25), 
            inherit.aes = FALSE, lty = 3) + 
  geom_line(mapping = aes(x= merge.ratio, y = duration.q75), 
            inherit.aes = FALSE, lty = 3) + 
  xlab('Merge ratio') + 
  ylab('Dive duration (min)') + 
  facet_wrap(~deep, ncol = 1, scales = 'free') + 
  theme_few() + 
  theme(panel.background = element_blank(), strip.placement = 'inside') + 
  ggtitle('(Solid: mean across tags; Dotted: Est. 25% and 75% quantiles)')

ggsave(pl, filename = file.path(diagnostic.dir, 'durations_per_ratio.pdf'),
       type = 'cairo')




#
# tag plots and save labels
#

# merge ratios selected
merge.ratios.ideal = rep(.6, length(tag.records))

mapply(FUN = function(record, ratio) {
  
  # find the labeled depth series that best matches the desired merge ratio
  merge.ind = which.min(abs(
    sapply(record$label.seq, function(lab) lab$merge.ratio) - ratio
  ))
  
  # extract dive labels
  labs = record$label.seq[[merge.ind]]$labels
  
  # plot dive record
  pl = tagplot(depths = record$depths, depth.bins = depth.bins, 
               dives.labeled = labs)
  
  # save plot of dive record 
  ggsave(pl, filename = file.path(out.dir, 
                                  paste(record$name, '.pdf', sep = '')),
         dpi = 'print', height = 12, width = 12*12, type = 'cairo',
         limitsize = FALSE)
  
  # save dive labels
  saveRDS(labs, 
          file = file.path(out.dir, 
                           paste(record$name, 'labels.rds', sep = '_'))
          )
  
}, tag.records, merge.ratios.ideal)
