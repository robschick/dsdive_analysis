library(lubridate)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(dsdive)
library(MASS)
library(yaml)

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

#
# roughly split dive data
#

# determine segments of underwater periods
surface.inds = which(depths$depth.bin == 1)
uw.periods = data.frame(
  start.ind = surface.inds[1:(length(surface.inds)-1)],
  end.ind = surface.inds[-1]
)

rm(surface.inds)


#
# identify long and deep dive segments
#

# dive summaries
dive.summaries = data.frame(t(apply(uw.periods, 1, function(r) {
  c(duration.min = as.numeric(diff(r)*5), 
    max.depth = depth.bins[max(depths$depth.bin[r[1]:r[2]]),1])
})))

# indices for deep dives
deep.inds = which(dive.summaries$max.depth >= 800)


#
# load manual labels
#

load('data/raw/zc84/dive_labels_simple_800_surfaceonly.RData')


#
# manually label remaining deep dives data from rough extraction
#

n = length(deep.inds)

unlabeled.dives = which(!complete.cases(dive.labels))

for(i in unlabeled.dives) {
  
  # get indices of dive observation
  inds.range = uw.periods[deep.inds[i],]
  inds = inds.range$start.ind:inds.range$end.ind
  
  # build dive object
  d = list(
    depths = depths$depth.bin[inds],
    times = as.numeric(depths$Date[inds]),
    label = 1:length(inds)
  )
  class(d) = 'dsobs'
  
  
  
  # build plot of dive
  pl = plot(x = d, depth.bins = depth.bins, time.as.POSIXct = TRUE, 
            errorbars = TRUE) + 
    geom_label(mapping = aes(x = t.start, y = depth.mid, label = ind), 
              data = data.frame(ds.df(depths = d$depths, times = d$times, 
                                      depth.bins = depth.bins, 
                                      time.as.POSIXct = TRUE),
                                ind = 1:length(inds)), inherit.aes = FALSE)
  
  # display dive
  print(pl)
  
  # solicit dive labels
  input = readline('Enter "n" to skip dive, or press enter to continue: ')
  if(input!='n') {
    
    input = readline('Enter start and end indices of dive: ')
    
    inds = as.numeric(unlist(strsplit(input, ' ')))
    
    # save dive labels
    dive.labels[i,c('d.start','d.end')] = inds
    
  }
  
}

# save manual labels
save(uw.periods, deep.inds, dive.labels, 
     file = 'sketches/2020-03-23_g1g2_eda/dive_labels_simple_800_all.RData')