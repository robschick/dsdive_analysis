# build and save simulated dive trajectories
#
# Simulation goal is to look at parameter recovery for variable observation 
# windows

library(dsdive)
library(ggplot2)
library(dplyr)
library(ggthemes)
library(yaml)


#
# simulation function
#

sim.gen = function(beta, lambda, T1.params, T2.params, N, t.win, 
                   out.path = NULL, seed = NULL, require.deep = FALSE,
                   known.end = FALSE) {
  # Parameters:
  #  beta, lambda, T1.params, T2.params - model parameters
  #  N - number of dives to simulate
  #  t.win - vector with number of seconds between observing the dive
  #  out.path - directory in which to save simulated trajectories
  #  seed - seed from which to start simulation
  #  require.deep - TRUE to rejection sample simulations until observed depth
  #    is at least 1,000m
  #  known.end - TRUE to set the time of the final observation (e.g., at the 
  #    surface) equal to the true end time of the dive
  
  if(!is.null(seed)) {
    set.seed(seed)
  }
  
  # simulate family of dives 
  x = replicate(N, {
    
    max.attempts = 1e5
    sampled.deep = FALSE
    
    attempts = 0
    while(!sampled.deep) {
      
      # sample stage durations
      T1 = 60 * rgamma(n = 1, shape = T1.params[1], rate = T1.params[2])
      T2 = 60 * rgamma(n = 1, shape = T2.params[1], rate = T2.params[2])
      
      # simulate dive
      d = dsdive.fwdsample.dive(depth.bins = depth.bins, beta = beta, 
                                lambda = lambda, t0 = 0, steps.max = 1e3, 
                                T1 = T1, T2 = T2)
      
      # observe dive at different time intervals
      d.obs = lapply(t.win, function(t.win) {
        # determine observation times
        if(known.end) {
          d.end = max(d$times)
          t.obs = unique(c(seq(from = 0, to = d.end, by = t.win), d.end))
        } else {
          t.obs = seq(from = 0, to = max(d$times) + t.win, by = t.win)
        }
        # observe dive
        dsdive.observe(depths = d$depths, times = d$times, 
                       stages = d$stages, t.obs = t.obs)
      })
      
      # check to see if sampled dive is deep
      if(require.deep) {
        attempts = attempts + 1
        sampled.deep = all(sapply(d.obs, function(dobs) { 
          max(depth.bins[dobs$depths,1]) >= 1e3
        }))
        if(attempts >= max.attempts) {
          break
        }
      } else {
        sampled.deep = TRUE
      }
      
    }
    
    if(!sampled.deep) {
      stop('Deep dive not successfully sampled')
    }
      
    # extract empirical stages times
    t.stages = d$times[c(FALSE, diff(d$stages)==1)]
    
    # sample dive
    list(list(
      dive = d,
      dive.obs = d.obs,
      t.stages = t.stages
    ))
  })

  # package parameters
  params = list(beta = beta, lambda = lambda, T1.params = T1.params, 
                T2.params = T2.params)
  
  # save dives to disk
  if(!is.null(out.path)) {
    
    # build output paths
    truth.path = file.path(out.path, 'truth')
    params.path = file.path(out.path, 'params')
    obs.path = file.path(out.path, paste('observations', t.win, sep ='_'))
    
    # sim_base_30
    
    # create output directories
    dir.create(path = truth.path, recursive = TRUE)
    dir.create(path = params.path, recursive = TRUE)
    for(o in obs.path) {
      dir.create(path = o, recursive = TRUE)
    }
    
    # save parameters
    save(params, seed, file = file.path(params.path, 'params.RData'))
    
    # generate and save yaml config file for each observation dataset
    lapply(1:length(t.win), function(i) {
      cfg = list(
        data = list(
          name = paste(c(strsplit(out.path, '/')[[1]][-1], t.win[i]), 
                       collapse = '_'),
          path = obs.path[i],
          tstep = t.win[i],
          file_patterns = list(
            dive = 'dive',
            depths = 'depths'
          ),
          description = paste(
            'Observations of simulated dives.  Dives are observed every ',
            t.win[i], ' seconds.  Complete dives are stored in ../truth.',
            sep = ''
          )
        )
      )
      
      write_yaml(cfg, file = file.path('conf', 'data', 
                                       paste(cfg$data$name, '.yaml', sep = '')))
    })
    
    # output dive information
    for(i in 1:length(x)) {
      
      write.csv(
        x = data.frame(x[[i]]$dive[c('depths','durations','times','stages')]), 
        file = file.path(truth.path, paste('dive', i, '.csv', sep = '')),
        row.names = FALSE)
      
      write.csv(
        x = depth.bins, 
        file = file.path(truth.path, paste('depths', i, '.csv', sep = '')),
        row.names = FALSE)
      
      for(j in 1:length(obs.path)) {
        write.csv(
          x = data.frame(x[[i]]$dive.obs[[j]][c('depths','times','stages')]), 
          file = file.path(obs.path[j], paste('dive', i, '_obs.csv', sep = '')),
          row.names = FALSE)
        
        write.csv(
          x = depth.bins, 
          file = file.path(obs.path[j], paste('depths', i, '.csv', sep = '')),
          row.names = FALSE)
      }
      
    }
    
  }
  
  list(dives.sim = x, params = params)
}


#
# base simulation configuration
#

n.sim = 100
seed = 2019

# load SATTAG depth template
depth.bins = read.csv('data/depth_template.csv')

# default parameters
beta = c(.9,.1)
lambda = c(1.5, .3, .8)

# mean and sd. for time in stage
T1.mean = 10
T1.sd = 2
T2.mean = 15
T2.sd = 2

# convert mean and sd. for time in stage to gamma shape/rate parameters
T1.params = c(T1.mean^2/T1.sd^2, T1.mean/T1.sd^2)
T2.params = c(T2.mean^2/T2.sd^2, T2.mean/T2.sd^2)

# curve(dgamma(x = x, shape = T1.params[1], rate = T1.params[2]),
#       from = 0, to = 20)
# curve(dgamma(x = x, shape = T2.params[1], rate = T2.params[2]),
#       from = 0, to = 20)

lambda.tyack = c(1.5, .3, .7)

# mean and sd. for time in stage
T1.mean = 17.3
T1.sd = 3
T2.mean = 32.8
T2.sd = 7.6

# convert mean and sd. for time in stage to gamma shape/rate parameters
T1.params.tyack = c(T1.mean^2/T1.sd^2, T1.mean/T1.sd^2)
T2.params.tyack = c(T2.mean^2/T2.sd^2, T2.mean/T2.sd^2)


# simulation based on Tyack paper
tyack.series = sim.gen(
  beta = beta, lambda = lambda.tyack, T1.params = T1.params.tyack, 
  T2.params = T2.params.tyack,
  N = n.sim, out.path = file.path('data', 'sim', 'tyack_alldeep_more'), 
  seed = seed, t.win = c(.5,1,5) * 60, require.deep = TRUE
)

# non-deep simulation based on Tyack paper
tyack.series.free = sim.gen(
  beta = beta, lambda = lambda.tyack, T1.params = T1.params.tyack, 
  T2.params = T2.params.tyack,
  N = n.sim, out.path = file.path('data', 'sim', 'tyack_more_known_end'), 
  seed = seed, t.win = c(.5,1,5) * 60, require.deep = FALSE, known.end = TRUE
)


# test the simulation settings
base.series = sim.gen(
  beta = beta, lambda = lambda, T1.params = T1.params, T2.params = T2.params,
  N = n.sim, out.path = file.path('data', 'sim', 'base_alldeep_more'), 
  seed = seed, t.win = c(.5,1,5) * 60, require.deep = TRUE
)

# test the simulation settings with faster diving rates
deeper.series = sim.gen(
  beta = c(.95, .1), lambda = c(1.8, .3, 1.2), 
  T1.params = T1.params, T2.params = T2.params,
  N = n.sim, out.path = file.path('data', 'sim', 'deeper_alldeep_more'), 
  seed = seed, t.win = c(.5,1,5) * 60, require.deep = TRUE
)


# mean and sd. for time in stage
T1.mean = 10
T1.sd = .5
T2.mean = 15
T2.sd =.5

# convert mean and sd. for time in stage to gamma shape/rate parameters
T1.params = c(T1.mean^2/T1.sd^2, T1.mean/T1.sd^2)
T2.params = c(T2.mean^2/T2.sd^2, T2.mean/T2.sd^2)

# test the simulation settings with faster diving rates
sharper.series = sim.gen(
  beta = c(.95, .1), lambda = c(1.8, .3, 1.2), 
  T1.params = T1.params, T2.params = T2.params,
  N = n.sim, out.path = file.path('data', 'sharper'), seed = seed, 
  t.win = c(.5,1,5) * 60
)


# function to plot dives in series
pl.dives = function(series) {
  for(dive.ind in 1:n.sim) {
    pl = plot(x = series$dives.sim[[dive.ind]]$dive.obs[[3]],
              depth.bins = depth.bins,
              stages = series$dives.sim[[dive.ind]]$dive.obs[[3]]$stages,
              errorbars = TRUE,
              imputed.list = series$dives.sim[[dive.ind]]$dive) +
      ggtitle(dive.ind)
    print(pl)
  }
}

# pl.dives(base.series)
# pl.dives(deeper.series)
pl.dives(sharper.series)

# function to extract max depths and durations
df.depths_durations = function(series) {
  data.frame(
    dive.id = 1:n.sim, 
    max.depth = sapply(series$dives.sim, function(d) {
      depth.bins[max(d$dive$depths),1]}),
    dive.duration =sapply(series$dives.sim, function(d) {
      d$dive$times[length(d$dive$times)]}) 
  )
}

df.base = df.depths_durations(base.series)
df.deeper = df.depths_durations(deeper.series)
df.sharper = df.depths_durations(sharper.series)


# review max depths

plot(density(df.base$max.depth), xlab = 'Max depth (m)', 
     main='Baseline (black) and Deeper (red) dives', ylim = c(0,.0025))
lines(density(df.deeper$max.depth), col = 2)
lines(density(df.sharper$max.depth), col = 3)

plot(density(df.base$dive.duration/60), xlab = 'Dive duration (min)', 
     ylim = c(0,.04), main='Baseline (black) and Deeper (red) dives')
lines(density(df.deeper$dive.duration/60), col = 2)
lines(density(df.sharper$dive.duration/60), col = 3)


# function to extract number of observations per stage
df.stage_obs = function(series) {
  do.call(rbind, lapply(1:length(series$dives.sim), function(ind) {
    d = series$dives.sim[[ind]]
    # process observations
    do.call(rbind, lapply(d$dive.obs, function(o) {
      df = data.frame(
        dive.id = ind,
        tstep = diff(o$times[1:2]),
        table(o$stages)
      ) %>% mutate(prop = Freq / sum(Freq))
      colnames(df)[3] = 'stage'
      df
    }))
  }))
}

stage_obs.base = df.stage_obs(base.series)
stage_obs.deeper = df.stage_obs(deeper.series)

# review distribution of observations per stage

ggplot(stage_obs.base %>% filter(stage==3), 
       aes(x = prop, col = factor(tstep))) + 
  # stat_density(geom = 'line') +
  stat_ecdf(geom = 'line') + 
  theme_few() + 
  theme(panel.border = element_blank())


# 
# #
# # configure and generate simulations
# #
# 
# beta1.seq = c(.8, .9, .95)
# beta2.seq = 1 - beta1.seq
# lambda1.seq = c(1, 1.5, 2)
# lambda2.seq = c(.1, .3, .5, .7)
# lambda3.seq = c(.6, .8, 1, 1.2)
# 
# 
# beta1.series = lapply(beta1.seq, function(theta) {
#   sim.gen(beta = c(theta, beta[2]), lambda = lambda, T1.params = T1.params, 
#           T2.params = T2.params, N = n.sim, 
#           out.path = file.path('data', paste('beta1_', theta, sep ='')), 
#           seed = seed)
# })
# 
# beta2.series = lapply(beta2.seq, function(theta) {
#   sim.gen(beta = c(beta[1], theta), lambda = lambda, T1.params = T1.params, 
#           T2.params = T2.params, N = n.sim, 
#           out.path = file.path('data', paste('beta2_', theta, sep ='')), 
#           seed = seed)
# })
# 
# lambda1.series = lapply(lambda1.seq, function(theta) {
#   sim.gen(beta = beta, lambda = c(theta, lambda[-1]), T1.params = T1.params, 
#           T2.params = T2.params, N = n.sim, 
#           out.path = file.path('data', paste('lambda1_', theta, sep ='')), 
#           seed = seed)
# })
# 
# lambda2.series = lapply(lambda2.seq, function(theta) {
#   sim.gen(beta = beta, lambda = c(lambda[1], theta, lambda[3]), 
#           T1.params = T1.params, T2.params = T2.params, N = n.sim, 
#           out.path = file.path('data', paste('lambda2_', theta, sep ='')), 
#           seed = seed)
# })
# 
# lambda3.series = lapply(lambda3.seq, function(theta) {
#   sim.gen(beta = beta, lambda = c(lambda[1:2], theta), 
#           T1.params = T1.params, T2.params = T2.params, N = n.sim, 
#           out.path = file.path('data', paste('lambda3_', theta, sep ='')), 
#           seed = seed)
# })
# 
