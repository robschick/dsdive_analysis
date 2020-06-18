
#
# flatten data structures
#

depths.flat = do.call(c, apply(dives.processed$dive.ranges, 1, function(r) { 
  dives.processed$depths[r['start.ind']:r['end.ind']]
}))

times.flat = do.call(c, apply(dives.processed$dive.ranges, 1, function(r) { 
  dives.processed$times[r['start.ind']:r['end.ind']]
}))

dive.lengths = apply(dives.processed$dive.ranges, 1, function(r) {
  r['end.ind'] - r['start.ind'] + 1
})

start.inds = as.numeric(c(1, 1 + cumsum(dive.lengths)))

dives.processed$dive.ranges$length = dive.lengths
dives.processed$dive.ranges$start.ind.flat = start.inds[
  1:(length(start.inds)-1)
]
  

# reverse-map to associate T with E

E_map = matrix(nrow = nrow(Tmat), ncol = 2)
for(i in 1:nrow(E_map)) {
  start_ind = which(dives.processed$endpoint.inds$dive_start == i)
  if(length(start_ind) == 0) {
    E_map[i,1] = 1
  } else {
    E_map[i,1] = start_ind
  }
  
  end_ind = which(dives.processed$endpoint.inds$dive_end == i)
  if(length(end_ind) == 0) {
    E_map[i,2] = 1
  } else {
    E_map[i,2] = end_ind
  }
}
  
#
# associate data and priors with model
#

deep.inds = dives.processed$dive.ranges$type == 1
shallow.inds = dives.processed$dive.ranges$type == 2

consts = list(
  M = nrow(depth.bins),
  N_endpoints = nrow(dives.processed$endpoint.inds),
  N_ranges_deep = sum(dives.processed$dive.ranges$type == 1),
  N_ranges_shallow = sum(dives.processed$dive.ranges$type == 2),
  tstep = 300,
  widths = depth.bins$halfwidth * 2,
  pi_priors = rbind(pi1.prior, pi2.prior, pi1.prior, pi2.prior),
  lambda_priors = rbind(lambda1.prior, lambda2.prior, lambda3.prior, 
                        lambda1.prior.shallow, lambda2.prior.shallow),
  dive_priors_deep = as.matrix(dives.processed$dive.ranges[deep.inds,]),
  dive_priors_shallow = as.matrix(dives.processed$dive.ranges[shallow.inds,]),
  stage_duration_priors = rbind(T1.prior, T2.prior, T1.prior.shallow),
  delta = 1e-10,
  E_priors = as.matrix(dives.processed$endpoint.inds),
  E_map = E_map,
  training_dive = dives.processed$dive.flags
)

dat = list(
  depths = depths.flat,
  times = times.flat
)

inits = list(
  pi = c(params.deep$pi, params.shallow$pi),
  lambda = c(params.deep$lambda, params.shallow$lambda),
  xi = xi,
  T = Tmat,
  xi = xi,
  E = apply(consts$E_priors, 1, function(r) { mean(r[c('t_lwr', 't_upr')]) })
)

inits$logit_pi = qlogis(inits$pi)
inits$log_lambda = log(inits$lambda)
inits$logit_E = qlogis((inits$E - consts$E_priors[,'t_lwr']) / 
  (consts$E_priors[,'t_upr'] - consts$E_priors[,'t_lwr']))

model = nimbleModel(
  code = modelCode,
  constants = consts,
  data = dat,
  inits = inits,
  name = 'ctdsDives'
)

cmodel = compileNimble(model, projectName = 'ctdsDives', resetFunctions = TRUE)

cmodel$calculate()