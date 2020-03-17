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
# manually label dive data from rough extraction
#

n = length(deep.inds)


dive.labels = data.frame(dive.summaries[deep.inds,], d.start = NA, d.end = NA)

for(i in 1:n) {
  
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
    
    # save dive labels
    dive.labels[i,c('d.start','d.end')] = c(1,length(inds))
  }
  
}

dive.labels$uw.period = deep.inds

# save manual labels
save(uw.periods, deep.inds, dive.labels, 
     file = 'data/raw/zc84/dive_labels_simple_800_surfaceonly.RData')

load('data/raw/zc84/dive_labels_simple_800_surfaceonly.RData')


#
# extract model parameters from extraction
#

labeled.dives = complete.cases(dive.labels)

dives.complete = dive.labels[labeled.dives,]


#
# save dives to csv formats
#


dives.to.export = which(dives.complete$d.start==1)

o = file.path('data', 'zc84_800')
dir.create(o, recursive = TRUE)

# export each dive to csv formats
for(dive.ind in dives.to.export) {
  
    dive = dives.complete[dive.ind,]
    
    dive.id = dive$uw.period
  
    # get depth record for dive
    d = depths %>% 
      # filter data for appropriate underwater period
      slice(uw.periods$start.ind[dive.id]:uw.periods$end.ind[dive.id]) %>% 
      # only include the diving portions of the uw period
      slice(dive$d.start:dive$d.end)
  
    # save csv with dive profile
    write.csv(
      x = data.frame(
        depths = d$depth.bin,
        # times of observations, in seconds
        times = as.numeric(d$Date)
      ),
      file = file.path(o, paste('dive', dive.id, '.csv', sep = '')),
      row.names = FALSE
    )
    
    # save csv with depth bins
    write.csv(
      x = depth.bins, 
      file = file.path(o, paste('depths', dive.id, '.csv', sep = '')), 
      row.names = FALSE
    )
}


#
# write config file for dataset
#

cfg = list(
  data = list(
    name = 'zc84_800',
    path = 'data/zc84_800',
    tstep = 300,
    file_patterns = list(
      dive = 'dive',
      depths = 'depths'
    ),
    description = paste(
      'Deep dives (dives over 800m) from zc84.  Dives are additionally',
      'filtered such that they begin and end in a surface bin.  Dives are',
      'stored in separate files, but observations are timestamped so that the',
      'dataset could support sequential modeling.  The dives are preprocessed',
      'so that the observations all use the same exact depth bins.  Depth bins',
      'are reported once every 300 seconds.',
      sep = ' '
    )
  )
)

write_yaml(cfg, file = file.path('conf','data','zc84_800.yaml'))
