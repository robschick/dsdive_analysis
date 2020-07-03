library(dplyr)
library(ggplot2)
library(ggthemes)


# list all message bin and depth files
messages.files = dir(path = file.path('data', 'raw'), pattern = 'seriesrange_', 
                     full.names = TRUE)
depths.files = dir(path = file.path('data', 'raw'), pattern = 'series_', 
                   full.names = TRUE)

# load message bins
messages = read.csv(messages.files[[1]])

# load observations
depths = read.csv(depths.files[[1]])

# add message ids to each depth record
depths$message.id = findInterval(depths$Date, messages$End) + 1

# extract dataset
dat.all = depths %>% 
  left_join(messages %>% mutate(message.id = 1:nrow(messages)), 
            by = 'message.id') %>% 
  select('message.id', 'MaxDepth', 'Depth', 'DRange') %>% 
  select(-message.id) %>%
  group_by(MaxDepth) %>%
  unique() %>%
  arrange(Depth) %>%
  mutate(bin.ind = 1:length(Depth)) %>%
  ungroup()

# there are only two messages that have complete records
unique(
depths %>% 
  left_join(messages %>% mutate(message.id = 1:nrow(messages)), 
            by = 'message.id') %>% 
  select('message.id', 'MaxDepth', 'Depth', 'DRange') %>% 
  group_by(message.id) %>%
  unique() %>%
  arrange(Depth) %>%
  mutate(bin.ind = 1:length(Depth)) %>%
  filter(max(bin.ind)==16) %>% 
  select(message.id)
  )

# get subset where all message bins are available for a max depth
dat.complete = dat.all %>% 
  group_by(MaxDepth) %>% 
  filter(max(bin.ind)==16) %>%
  ungroup()

# filter subset where not all message bins are available
dat.test = dat.all %>% 
  group_by(MaxDepth) %>% 
  filter(max(bin.ind)!=16) %>%
  ungroup()

# verify high correlation between Depth and DRange, which is good for prediction
cor(dat.complete %>% select(Depth, DRange))

# build a predictive model for (Depth, DRange) given MaxDepth and bin index
fit.by.ind = lapply(1:16, function(ind) {
  lm(cbind(Depth,DRange) ~ MaxDepth, dat.complete %>% filter(bin.ind==ind))
})


#
# test model
#

# verify we do not have a case where there are more than 16 bins
all(dat.test %>% 
  group_by(MaxDepth) %>%
  summarise(n.inds = max(bin.ind)) %>%
  select(n.inds) < 16)

dir.create(file.path('imputed_bins', 'bins'), recursive = TRUE)

max.depths = unique(dat.test$MaxDepth)
for(m in max.depths) {
  
  # predict depth bins for this depth
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
  
  # extract observed depth bins
  obs = dat.test %>% filter(MaxDepth==m) %>% select(Depth, DRange)
  
  # determine which group each of the observations is closest to;
  # use DRange because there is greater group separation vs. using Depth
  obs$group = apply(outer(obs$DRange, pred.means$mean.drange, '-')^2, 
                     1, function(x) which.min(x))
  
  # bias correct DRanges by projecting predicted DRanges onto an empirical 
  # local, linear subspace for DRanges; predicted depths look ok
  pred.bins = do.call(rbind, lapply(1:4, function(g) {
    fit.obs = lm(DRange~Depth, obs %>% filter(group==g))
    df = pred.bins %>% filter(group==g)
    df$DRange = predict(fit.obs, df)
    df
  }))
  
  # plot truth only
  pl.truth = rbind(
    obs %>% mutate(Bin = 'Observed')
  ) %>% ggplot(aes(x = Depth, y = DRange, color = Bin, shape = Bin)) + 
    geom_point(alpha = .7) + 
    scale_color_brewer(type = 'qual', palette = 'Dark2') + 
    ggtitle(paste('MaxDepth = ', m, sep = '')) + 
    theme_few() + 
    theme(panel.border = element_blank())
  
  # plot predictions alongside truth
  pl = rbind(
    obs %>% mutate(Bin = 'Observed'),
    pred.bins %>% mutate(Bin = 'Predicted')
  ) %>% ggplot(aes(x = Depth, y = DRange, color = Bin, shape = Bin)) + 
    geom_point(alpha = .7) + 
    scale_color_brewer(type = 'qual', palette = 'Dark2') + 
    ggtitle(paste('MaxDepth = ', m, sep = '')) + 
    theme_few() + 
    theme(panel.border = element_blank())
  
  ggsave(pl, file = paste('imputed_bins/imputed_bins_', m, '.png', sep = ''), 
         width = 5, height = 5)
  
  ggsave(pl.truth, file = paste('imputed_bins/imputed_bins_observed_', m, 
                                '.png', sep = ''), 
         width = 5, height = 5)
  
  pred.bins = pred.bins[,1:2]
  colnames(pred.bins) = c('center', 'halfwidth')
  
  # save depth bins
  write.csv(pred.bins, file = file.path('imputed_bins', 'bins', 
                                        paste('bin', m, '.csv', sep='')), 
            row.names = FALSE)
}




pl = dat.test %>% select(Depth, DRange) %>% 
  ggplot(aes(x = Depth, y = DRange)) + 
  geom_point() +
  theme_few() + 
  theme(panel.border = element_blank())

ggsave(pl, file = 'imputed_bins/all_depths.png', width = 5, height = 5)
