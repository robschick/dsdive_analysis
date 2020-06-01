library(coda)
library(dplyr)
library(composr)

# clear workspace
rm(list = ls())


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

groups = list(
  data = 'zc84_800_covariates',
  observation_model = 'uniform_systematic',
  priors = 'tyack_cov_duration_priors',
  sampler = 'prod',
  subset = 'all_dives',
  validation= 'holdout_half'
)

# build configuration
cfg = compose_cfg(file = file.path('conf', 'config.yaml'), groups = groups)
rm(groups)

# output paths
out.dir = file.path(cfg$base_paths$fit, cfg$data$name, cfg$subset$name, 
                    cfg$validation$name, cfg$observation_model$name, 
                    cfg$priors$name)
dir.create(out.dir, recursive = TRUE)

fig.dir = file.path(out.dir, cfg$sub_paths$figures, cfg$sub_paths$posteriors)
dir.create(fig.dir, recursive = TRUE)

# load data and utility functions
source(file.path('scripts', 'utils', 'datafns.R'))
source(file.path('scripts', 'utils', '85pct_rule.R'))


#
# load data
#

dives.obs = dives.load(path = cfg$data$path, 
                       dive_pattern = cfg$data$file_patterns$dive,
                       depth_pattern = cfg$data$file_patterns$depths, 
                       covariates = cfg$data$file_patterns$covariates)


#
# extract lists
#

depth.bins = dives.obs$dives[[1]]$depth.bins
dives.obs.list = lapply(dives.obs$dives, function(d) d$dive)

t.stages.list = lapply(dives.obs.list, function(d) {
  seq(from = d$times[1], to = d$times[length(d$times)], length.out = 4)[2:3]
})


#
# load posterior samples and priors
#

load(file.path(out.dir, cfg$base_names$fit))
load(file.path(out.dir, cfg$base_names$ctmc_priors))


#
# posterior diagnostics
#

attach(state)
names(state)

it = 1:sum(!is.na(state$theta$beta1[,1]))
burn = 1:1e3
post = it[-burn]

plot.coda = function(series, priors) {
  
  priors = switch(series, 
                  beta1 = beta.priors[[1]],
                  beta2 = beta.priors[[2]],
                  alpha1 = alpha.priors[[1]],
                  alpha2 = alpha.priors[[2]], 
                  alpha3 = alpha.priors[[3]])
  
  par(mfrow = c(ncol(theta[[series]]),2), mar = c(3,3,2,1) + .1)
  
  for(i in 1:ncol(theta[[series]])) {
    # print(plot(mcmc(theta[[series]][post,])))
    m = mcmc(theta[[series]][post,i])
    hpd = HPDinterval(m)
    print(plot(m, trace = TRUE, density = FALSE, 
               main = paste(series, i), auto.layout = FALSE))
    print(plot(m, trace = FALSE, density = TRUE, auto.layout = FALSE, 
               main = 'Density'))
    curve(dnorm(x = x, mean = priors$mu[i], sd = priors$sd[i]), add = TRUE, 
                col = 2)
    abline(v = priors$mu[i], lty = 3, col = 2)
    abline(v = hpd, col = 'grey60')
  }
  
  print(summary(mcmc(theta[[series]][post,])))
  print(effectiveSize(mcmc(theta[[series]][post,])))
  cat('\n95% HPD\n')
  print(HPDinterval(mcmc(theta[[series]][post,]), prob = .95))
  cat('\n90% HPD\n')
  print(HPDinterval(mcmc(theta[[series]][post,]), prob = .9))
}


tgt.seq = c('beta1', 'beta2', 'alpha1', 'alpha2', 'alpha3')

for(tgt in tgt.seq) {
  sink(file.path(fig.dir, paste(tgt, '.txt', sep = '')))
  # png(file.path(fig.dir, paste(tgt, '.png', sep = '')), 
  #     height = 480*2, width = 480*4)
  pdf(file.path(fig.dir, paste(tgt, '.pdf', sep = '')), 
      width = 7*1.5, height = 7)
  plot.coda('beta1')
  dev.off()
  sink()
}


detach(state)
