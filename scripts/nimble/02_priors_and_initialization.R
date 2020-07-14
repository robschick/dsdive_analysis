# density fitting
library(MASS, lib.loc = c('singularity/libs', .libPaths()))

# 85% rule script, for dive priors
source(file.path('scripts', 'utils', '85pct_rule.R'))

# load flattened data
nim_pkg = readRDS(file.path('data', 'tag_endpoints', 'flattened_endpoints.rds'))

# output directory
out.dir = file.path('output', 'multitag')
dir.create(out.dir, recursive = TRUE)

# template bins
depth.bins = read.csv(file.path('data', 'imputed_bins', 'template.csv'))


#
# priors for dive durations
#

# reformat dive data for 85% rule function
dives.obs = apply(nim_pkg$consts$dive_relations, 1, function(r) {
  list(
    depth.bins = depth.bins,
    dive = list(
      depths = nim_pkg$data$depths[r['depth_first']:r['depth_last']],
      times = nim_pkg$data$times[r['depth_first']:r['depth_last']] - 
        nim_pkg$data$times[r['depth_first']]
    )
  )
})

# estimate stage durations by tag
durations.est = cbind(
  times.stages(dives.obs) * 60,
  tag = nim_pkg$consts$dive_relations[,'tag']
)

colnames(durations.est) = c('sub.time.sec', 'bottom.time.sec', 'tag')


duration_priors = do.call(
  rbind, lapply(sort(unique(durations.est$tag)), function(i) {
  # identify estimates for ith tag
  inds = durations.est$tag == i
  # estimate dive durations and correlation on log scale
  res = c(
    fitdistr(durations.est$sub.time.sec[inds], 'lognormal')$estimate,
    fitdistr(durations.est$bottom.time.sec[inds], 'lognormal')$estimate,
    cov(log(durations.est$sub.time.sec[inds]), 
        log(durations.est$bottom.time.sec[inds]))
  )
  # format and return
  names(res) = c('G1_mean', 'G1_sd', 'G2_mean', 'G2_sd', 'cov_log')
  res
}))

# reformat duration priors for nimble
nim_pkg$consts$xi_prior_means = duration_priors[,c('G1_mean', 'G2_mean')]
nim_pkg$consts$xi_prior_covs = array(data = NA, 
                                     dim = c(nrow(duration_priors), 2, 2))
for(i in 1:nrow(duration_priors)) {
  nim_pkg$consts$xi_prior_covs[i,1,1] = duration_priors[i, 'G1_sd']^2
  nim_pkg$consts$xi_prior_covs[i,1,2] = duration_priors[i, 'cov_log']
  nim_pkg$consts$xi_prior_covs[i,2,1] = duration_priors[i, 'cov_log']
  nim_pkg$consts$xi_prior_covs[i,2,2] = duration_priors[i, 'G2_sd']^2
}



#
# random effect and covariate priors
#


# prior means for logit-pi model coefficients
nim_pkg$consts$beta_pi_prior_mean = rbind(
  c(intercept = qlogis(.95), sex = 0),
  c(intercept = qlogis(.5), sex = 0),
  c(intercept = qlogis(.05), sex = 0)
)

# prior uncertainty for logit-pi model coefficients
nim_pkg$consts$beta_pi_prior_sd = rbind(
  c(intercept = 1, sex = 1),
  c(intercept = 0, sex = 0),
  c(intercept = 1, sex = 1)
)

# prior distribution for logit-pi model random effect scale
nim_pkg$consts$sigma_pi_priors = rbind(
  c(shape = 2, rate = 1),
  c(shape = 0, rate = 0),
  c(shape = 2, rate = 1)
)

# prior means for log-lambda model coefficient
nim_pkg$consts$beta_lambda_prior_mean = rbind(
  c(intercept = log(1.5), sex = 0),
  c(intercept = log(.3), sex = 0),
  c(intercept = log(1), sex = 0)
)

# prior uncertainty for log-lambda model coefficients
nim_pkg$consts$beta_lambda_prior_sd = rbind(
  c(intercept = 1, sex = 1),
  c(intercept = 1, sex = 1),
  c(intercept = 1, sex = 1)
)

# prior distribution for log-lambda model random effect scale
nim_pkg$consts$sigma_lambda_priors = rbind(
  c(shape = 2, rate = 1),
  c(shape = 2, rate = 1),
  c(shape = 2, rate = 1)
)


#
# set initial values for random variables
#
  
nim_pkg$inits$endpoints = runif(
  n = nim_pkg$consts$N_endpoints, 
  min = nim_pkg$consts$endpoint_priors[,'t_lwr'],
  max = nim_pkg$consts$endpoint_priors[,'t_upr']
)

nim_pkg$inits$xi = as.matrix(
  durations.est[,c('sub.time.sec', 'bottom.time.sec')]
)
colnames(nim_pkg$inits$xi) = c('sub_time_sec', 'bottom_time_sec')

nim_pkg$inits$log_xi = log(nim_pkg$inits$xi)

nim_pkg$inits$T = cbind(
  T0 = nim_pkg$inits$endpoints[nim_pkg$consts$dive_relations[,'T0_endpoint']],
  T1 = NA,
  T2 = NA,
  T3 = nim_pkg$inits$endpoints[nim_pkg$consts$dive_relations[,'T3_endpoint']]
)

nim_pkg$inits$T[,'T1'] = nim_pkg$inits$T[,'T0'] + 
  nim_pkg$inits$xi[,'sub_time_sec']

nim_pkg$inits$T[,'T2'] = nim_pkg$inits$T[,'T1'] + 
  nim_pkg$inits$xi[,'bottom_time_sec']

# verify all time random effects are well ordered
all(
  nim_pkg$inits$T[,'T0'] < nim_pkg$inits$T[,'T1'],
  nim_pkg$inits$T[,'T1'] < nim_pkg$inits$T[,'T2'],
  nim_pkg$inits$T[,'T2'] < nim_pkg$inits$T[,'T3']
)


nim_pkg$inits$beta_pi = nim_pkg$consts$beta_pi_prior_mean

nim_pkg$inits$logit_pi = matrix(
  rep(nim_pkg$consts$beta_pi_prior_mean[,'intercept'], nim_pkg$consts$N_tags), 
  ncol = 3, byrow = TRUE
)

nim_pkg$inits$pi = plogis(nim_pkg$inits$logit_pi)


nim_pkg$inits$beta_lambda = nim_pkg$consts$beta_lambda_prior_mean

nim_pkg$inits$log_lambda = matrix(
  rep(nim_pkg$consts$beta_lambda_prior_mean[,'intercept'], 
      nim_pkg$consts$N_tags),
  ncol = 3, byrow = TRUE
)
  
nim_pkg$inits$lambda = exp(nim_pkg$inits$log_lambda)

nim_pkg$inits$sigma_pi = apply(nim_pkg$consts$sigma_pi_priors, 1, function(r) {
  rgamma(n = 1, shape = r['shape'], rate = r['rate'])
})

nim_pkg$inits$sigma_lambda = apply(nim_pkg$consts$sigma_lambda_priors, 1, 
                                   function(r) {
  rgamma(n = 1, shape = r['shape'], rate = r['rate'])
})

# save nimble-formatted data and inits
saveRDS(nim_pkg, file = file.path(out.dir, 'nim_pkg.rds'))
