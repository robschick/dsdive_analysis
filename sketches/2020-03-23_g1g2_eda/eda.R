library(lubridate)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(dsdive)
library(MASS)
library(yaml)


# output path

o = file.path('sketches', '2020-03-23_g1g2_eda')


#
# load and format depth series
#

depths.raw = read.csv('data/raw/zc84/ZcTag084_DUML_series_20200108.csv')

depths = depths.raw %>% mutate(
  Date = as.POSIXct(Date, origin = "1970-01-01", tz = "UTC")
)

rm(depths.raw)


#
# map all depths to standardized depth bins
#

depth.bins = read.csv('data/depth_template.csv')

depths$depth.bin = sapply(depths$Depth, 
                          function(d) which.min(abs(d - depth.bins$center)))

# build non-overlapping depth bins, for display
depth.bins = depth.bins %>% 
  mutate(bin.min = center - halfwidth, 
         bin.max = center + halfwidth,
         bin.ind = 1:n())
depth.bins$bin.min[-1] = depth.bins$bin.max[1:(nrow(depth.bins)-1)]

# format depth bin display labels
depth.bins = depth.bins %>% mutate(
  bin.range = paste('[', round(bin.min), ', ', round(bin.max), ')', 
                    sep =''),
  bin.range = ordered(bin.range, bin.range[bin.ind])
)


#
# load manual labels
#

load('sketches/2020-03-23_g1g2_eda/dive_labels_simple_800_all.RData')


#
# compute dive durations and time between dives
#

# augment dive summaries with start and end times, and surface periods
dive.lengths = do.call(rbind, apply(dive.labels, 1, function(r) {
  # depths-aligned starting index of underwater period
  depth.ind.start = uw.periods[r['uw.period'],1]
  # depths-aligned indices of dive within the underwater period
  dive.inds = depth.ind.start + r[c('d.start', 'd.end')] - 1
  # extract times associated with these indices
  times = depths$Date[unlist(dive.inds)]
  data.frame(cbind(t(r), t.start = times[1], t.end = times[2]))
})) %>% mutate(
  duration = t.end - t.start,
  surface.previous = t.start - lag(t.end)
)


#
# summary statistics for periods between deep dives
#

# identify underwater periods without deep dives
shallow.periods = setdiff(1:nrow(uw.periods), dive.labels$uw.period)

# extract depth bins for shallow dives
shallow.dives = lapply(shallow.periods, function(uw) {
  # indices of the underwater period
  dive.inds = seq(from = uw.periods[uw,1], to = uw.periods[uw,2])
  # depths 
  depth.bins$bin.range[depths$depth.bin[dive.inds]]
})


#
# eda fits
#

# gamma model for time between deep dives
fit.surface_interval = fitdistr(
  dive.lengths %>% 
    mutate(surface.previous = surface.previous/60) %>% 
    filter(is.finite(surface.previous), surface.previous > 0) %>%
    dplyr::select(surface.previous) %>% 
    unlist(), 
  'gamma')


#
# eda plots
#

# we see a low propensity for whales to start a deep dive too soon after a 
# deep dive.  however, after 30 mins (dotted line) we see a relatively small 
# increase in the dive duration with respect to the recovery time, but the 
# biggest thing to note really is that it looks like there is roughly a 
# minimum amount of time between deep dives.  there may be second-order 
# time dependence or other characteristics to explain the deep dives before 
# 30 minutes, but the previous statements hold for the bulk of dives.
#
# the implication is that sequences of deep dives may not be modeled well by a 
# HMM, but rather a semi-markov model.  of course, the gamma-like distribution 
# can result from a HMM where shallow dives are exponentially distributed, and 
# switching probability between classes is low 
lambda = 5
pl = dive.lengths %>% ggplot(aes(x = surface.previous/60, 
                            y = duration/60)) + 
  geom_point(alpha = .15) + 
  geom_quantile(method = 'rqss', quantiles = .5, col = 'black', 
                lambda = lambda) + 
  geom_quantile(method = 'rqss', quantiles = c(.1,.9), col = 'gray80', 
                lambda = lambda) + 
  geom_vline(xintercept = 30, lty = 3) + 
  xlab('Time since last deep dive (min.)') + 
  ylab('Dive duration (min.)') + 
  theme_few() + 
  theme(panel.border = element_blank())

ggsave(pl, filename = file.path(o, 'duration_time.png'), dpi = 'print')

# a gamma fit to the time between last deep dive looks reasonable with the 
# oversimplified HMM assumptions, and seem to pair well with Tyack (2006)
pl = dive.lengths %>% ggplot(aes(x = surface.previous/60)) + 
  stat_density(geom='line') + 
  stat_function(fun = dgamma, args = as.list(fit.surface_interval$estimate),  
                lty = 4) + 
  xlab('Time since last deep dive (min.)') + 
  ylab('Density') + 
  theme_few() + 
  theme(panel.border = element_blank())

ggsave(pl, filename = file.path(o, 'marginal_fit.png'), dpi = 'print')

# distribution of depth bins for shallow dives
pl = data.frame(x = depth.bins$bin.range[do.call(c, shallow.dives)]) %>% 
  mutate(total = n()) %>% 
  group_by(x) %>% 
  summarise(prob = length(x) / total[1]) %>% 
  ggplot(aes(x = x, y = prob)) + 
  xlab('Depth bin (m)') + 
  ylab('Empirical probability') + 
  geom_bar(stat = 'identity')  + 
  ggtitle('Depth bin distribution outside deep dives') + 
  theme_few() + 
  theme(panel.border = element_blank(), 
        text = element_text(size = 10))

ggsave(pl, filename = file.path(o, 'shallow_depth_bin_distribution.pdf'), 
       dpi = 'print', width = 8, height = 4)


#
# long time series plot of depth bins
#

df = data.frame(Date = depths$Date, 
                Depth = depth.bins$center[depths$depth.bin],
                bin.ind = depths$depth.bin) %>% 
  mutate(surf = Depth * (bin.ind==1)) 

pl = ggplot(df, aes(x = Date, y = Depth)) + 
  geom_line(col = 'grey80', lwd = .5) +
  geom_point(mapping = aes(x = Date, y = surf), 
             data = df %>% filter(surf > 0),
             size = .75) + 
  scale_y_reverse('Depth bin center (m)') +
  scale_x_datetime('Time', breaks = dateseq, date_labels = '%B %d') + 
  ggtitle('zc84 (Surface observations marked by solid points)') + 
  theme_few()

ggsave(pl, filename = file.path(o, 'long_tag_record.pdf'), 
       width = 150, height = 5, limitsize = FALSE)
