# modeling tools
library(dsdive, lib.loc = c('.', .libPaths()))
# plotting tools
library(ggplot2)
library(ggthemes)
library(ggpubr)


depth.bins = read.csv('data/depth_template.csv')


#
# baseline dive
#

T1 = 10*60
T2 = 15*60

set.seed(1918)

A = dsdive.fwdsample.dive(depth.bins = depth.bins, t0 = 0, steps.max = 1e3,
                          beta = c(.99, .05), 
                          lambda = c(1.25, .6, 1), 
                          T1 = T1, T2 = T2)


#
# faster-motion dive
#

T1 = 10*60
T2 = 15*60

set.seed(1918)

B = dsdive.fwdsample.dive(depth.bins = depth.bins, t0 = 0, steps.max = 1e3,
                          beta = c(.99, .05), 
                          lambda = c(2.5, 1.2, 2), 
                          T1 = T1, T2 = T2)


#
# less-directed dive
#

T1 = 20*60
T2 = 15*60

set.seed(1918)

C = dsdive.fwdsample.dive(depth.bins = depth.bins, t0 = 0, steps.max = 1e3,
                          beta = c(.8, .2), 
                          lambda = c(1.25, .6, 1), 
                          T1 = T1, T2 = T2)


#
# less-directed, faster-motion dive
#

T1 = 10*60
T2 = 15*60

set.seed(1918)

D = dsdive.fwdsample.dive(depth.bins = depth.bins, t0 = 0, steps.max = 1e3,
                          beta = c(.8, .2), 
                          lambda = c(2.5, 1.2, 2), 
                          T1 = T1, T2 = T2)


#
# arrange plots
#

yax = scale_y_reverse(limits = rev(c(0, max(rowSums(depth.bins)))))
xax = scale_x_continuous(limits = c(0, max(A$times, B$times, 
                                           C$times, D$times))/60)

plA = plot(x = A, depth.bins = depth.bins) + 
  xax + yax + 
  xlab('') + 
  guides(fill = 'none')

plB = plot(x = B, depth.bins = depth.bins) + 
  xax + yax + 
  ylab('') + 
  xlab('') + 
  guides(fill = 'none')

plC = plot(x = C, depth.bins = depth.bins) + 
  xax + yax +
  xlab('Time (min)') + 
  guides(fill = 'none')

plD = plot(x = D, depth.bins = depth.bins) + 
  xax + yax +
  xlab('Time (min)') + 
  ylab('') + 
  guides(fill = 'none')

pl = ggarrange(plA, plB, plC, plD, ncol = 2, nrow = 2, 
               labels = paste(LETTERS[1:4], ')', sep=''))

pl
#
# save plot
#

o = file.path('output', 'example_dives')

dir.create(o, recursive = TRUE)

ggsave(pl, filename = file.path(o, 'example_dives.png'), dpi = 'print', 
       width = 8, height = 4)
