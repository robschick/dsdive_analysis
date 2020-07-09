library(lubridate, lib.loc = c('singularity/libs', .libPaths()))
library(ggplot2, lib.loc = c('singularity/libs', .libPaths()))
library(ggthemes, lib.loc = c('singularity/libs', .libPaths()))
library(dplyr, lib.loc = c('singularity/libs', .libPaths()))

tagplot = function(depths, depth.bins, dives.labeled, date_range = NULL) {
  
  # build dataframe
  df = cbind(depths, mode = dives.labeled)
  
  # compute date range
  if(is.null(date_range)) {
    date_range = c(floor_date(with_tz(df$Date[1], 'America/New_York'), 'day'),
                   ceiling_date(with_tz(df$Date[nrow(df)], 'America/New_York'), 
                                'day'))
  } else {
    date_range = c(floor_date(with_tz(date_range[1], 'America/New_York'), 'day'),
                   ceiling_date(with_tz(date_range[2], 'America/New_York'), 
                                'day'))
    
    df = df %>% filter(Date >= date_range[1], Date <= date_range[2])
  }
  
  # compute date axis breaks
  dateseq = seq(from = date_range[1], to = date_range[2], by = 'day')
  
  # plot data
  pl = ggplot(df, aes(x = Date, y = depth.standardized)) + 
    # shading dive segments
    geom_rect(mapping = aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax), 
              data = df %>% 
                filter(mode != 0) %>% 
                group_by(mode) %>% 
                summarise(start = min(Date), 
                          end = max(Date),
                          ymin = depth.bins$center[1], 
                          ymax = max(depth.standardized)), 
              inherit.aes = FALSE, alpha = .15) + 
    # horizontal surface line
    geom_hline(yintercept = depth.bins$center[1], lty = 3) + 
    # time series plot of depth bins
    geom_line(col = 'grey80') + 
    geom_point(aes(col = mode > 0)) + 
    # dive/"not dive" coloring of points
    scale_color_brewer(type = 'qual', palette = 'Dark2') + 
    # axis and plot formatting
    scale_x_datetime('Time', breaks = dateseq, date_labels = '%B %d') + 
    scale_y_reverse() + 
    guides(col = 'none') + 
    ylab('Depth (m)') + 
    theme_few()
  
  pl
}