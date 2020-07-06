library(dplyr)
library(ggplot2)
library(ggthemes)


rm(list = ls())

# directory for imputed bins
out.dir = file.path('data', 'imputed_bins')
dir.create(out.dir, recursive = TRUE)

# directory for diagnostics
diagnostic.dir = file.path(out.dir, 'diagnostics')
dir.create(diagnostic.dir, recursive = TRUE)


#
# Load/Extract depth bins from all 4-hr message blocks from all sattag records
#

# list all message bin and depth files
messages.files = dir(path = file.path('data', 'raw'), pattern = 'seriesrange_', 
                     full.names = TRUE)
depths.files = dir(path = file.path('data', 'raw'), pattern = 'series_', 
                   full.names = TRUE)

# load/extract/aggregate depth bins from each sattag record 
bins_by_sattag = data.frame(
  do.call(rbind, mapply(function(message.file, depth.file) {
  
  # load metadata about each 4-hr message block
  messages = read.csv(message.file)
  
  # load observations
  depths = read.csv(depth.file)
  
  # add message block ids to each depth record
  depths$message.id = findInterval(depths$Date, messages$End) + 1
  
  # many message blocks share the same MaxDepth value; 
  # extract unique depth bins reported for each MaxDepth value
  as.matrix(depths %>% 
    left_join(messages %>% mutate(message.id = 1:nrow(messages)), 
              by = 'message.id') %>% 
    select('message.id', 'MaxDepth', 'Depth', 'DRange') %>% 
    select(-message.id) %>% 
    group_by(MaxDepth) %>% 
    unique() %>% 
    ungroup())
  
}, messages.files, depths.files)))


# merge depth bins from all sattag records; enumerate bin indices
bins_merged = bins_by_sattag %>% 
  group_by(MaxDepth) %>% 
  unique() %>% 
  arrange(Depth) %>% 
  mutate(bin.ind = 1:length(Depth)) %>% 
  ungroup()


#
#  Learn relationship between depth bin midpoint/range and bin number (1...16)
#

# subset where all message bins are available for a max depth
bins_complete = bins_merged %>% 
  group_by(MaxDepth) %>% 
  filter(max(bin.ind) == 16) %>% 
  ungroup()

# filter subset where not all message bins are available
bins_missing = bins_merged %>% 
  group_by(MaxDepth) %>% 
  filter(max(bin.ind) != 16) %>% 
  ungroup()

# verify high correlation between Depth and DRange, which is good for prediction
cor(bins_complete %>% select(Depth, DRange))

# build a predictive model for (Depth, DRange) given MaxDepth and bin index
fit.by.ind = lapply(1:16, function(ind) {
  lm(cbind(Depth,DRange) ~ MaxDepth, bins_complete %>% filter(bin.ind==ind))
})



#
# Output complete bins
#

max.depths = sort(unique(bins_complete$MaxDepth))

for(m in max.depths) {
  
  # extract observed depth bins associated with MaxDepth == m
  obs = bins_complete %>% filter(MaxDepth==m) %>% select(Depth, DRange)
  
  # format and save depth bins
  colnames(obs) = c('center', 'halfwidth')
  write.csv(obs, file = file.path(out.dir, 
                                  paste('bin', m, '.csv', sep='')), 
            row.names = FALSE)
  
  # print(c(rounded_centers = all(obs$center == round(obs$center)),
  #         half_rounded_halfwidths = all(obs$halfwidth == 
  #                                       round(obs$halfwidth/.5)*.5)
  # ))
  
  print(all(
    # centers are all rounded to nearest whole meter
    obs$center == round(obs$center),
    # midpoints are all rounded to nearest half meter
    obs$halfwidth == round(obs$halfwidth/.5)*.5
    
  ))
}

# save deepest observed depth bins 
depths.deepest = obs


#
#  Impute missing depth bins
#


# verify we do not have MaxDepth groups with more than 16 bins
all(bins_missing %>% 
  group_by(MaxDepth) %>%
  summarise(n.inds = max(bin.ind)) %>%
  select(n.inds) < 16)

max.depths = sort(unique(bins_missing$MaxDepth))
for(m in max.depths) {
  
  # predict depth bins associated with MaxDepth == m
  pred.bins = sapply(fit.by.ind, function(fit) {
    matrix(c(1, m), nrow = 1, ncol = 2) %*% coef(fit)
  })
  pred.bins = data.frame(Depth = pred.bins[1,], 
                         DRange = pred.bins[2,], 
                         group = rep(1:4,rep(4,4)))
  
  # get group means for predicted bins
  pred.means = pred.bins %>%  
    group_by(group) %>% 
    summarise(mean.drange = mean(DRange),
              mean.depth = mean(Depth))
  
  # extract observed depth bins associated with MaxDepth == m
  obs = bins_missing %>% filter(MaxDepth==m) %>% select(Depth, DRange)
  
  # determine which depth-group each of the observed bins is closest to;
  # use DRange because there is greater group separation vs. using Depth
  obs$group = apply(outer(obs$DRange, pred.means$mean.drange, '-')^2, 
                     1, function(x) which.min(x))
  
  # bias correct predicted DRanges by projecting predicted DRanges onto an 
  # empirical local, linear subspace determined by observed DRanges
  pred.bins = do.call(rbind, lapply(1:4, function(g) {
    fit.obs = lm(DRange~Depth, obs %>% filter(group==g))
    df = pred.bins %>% filter(group==g)
    df$DRange = predict(fit.obs, df)
    df
  }))
  
  # diagnostic: plot observed depth bins
  pl.truth = rbind(
    obs %>% mutate(Bin = 'Observed')
  ) %>% ggplot(aes(x = Depth, y = DRange, color = Bin, shape = Bin)) + 
    geom_point(alpha = .7) + 
    scale_color_brewer(type = 'qual', palette = 'Dark2') + 
    ggtitle(paste('MaxDepth = ', m, sep = '')) + 
    theme_few() + 
    theme(panel.border = element_blank())
  
  # diagnostic: plot predictions alongside observations
  pl = rbind(
    obs %>% mutate(Bin = 'Observed'),
    pred.bins %>% mutate(Bin = 'Predicted')
  ) %>% ggplot(aes(x = Depth, y = DRange, color = Bin, shape = Bin)) + 
    geom_point(alpha = .7) + 
    scale_color_brewer(type = 'qual', palette = 'Dark2') + 
    ggtitle(paste('MaxDepth = ', m, sep = '')) + 
    theme_few() + 
    theme(panel.border = element_blank())
  
  # output diagnostic plots
  ggsave(pl, file = file.path(diagnostic.dir, 
                              paste('imputed_bins_', m, '.png', sep = '')),
         width = 5, height = 5)
  ggsave(pl.truth, file = file.path(diagnostic.dir, 
                              paste('observed_bins_', m, '.png', sep = '')),
         width = 5, height = 5)
  
  # format and save depth bins
  pred.bins = pred.bins[,1:2]
  colnames(pred.bins) = c('center', 'halfwidth')
  pred.bins$center = round(pred.bins$center)
  pred.bins$halfwidth = round(pred.bins$halfwidth/.5)*.5
  write.csv(pred.bins, file = file.path(out.dir, 
                                        paste('bin', m, '.csv', sep='')), 
            row.names = FALSE)
}


# update absolute deepest set of depth bins
if(max(pred.bins$center) > max(depths.deepest$center)) {
  depths.deepest = pred.bins
}

# save absolute deepest set of depth bins as template
write.csv(depths.deepest, file = file.path(out.dir, 'template.csv'))
