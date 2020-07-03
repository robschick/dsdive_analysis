
#
# associate data and priors with model
#

consts = list(
  M = nrow(depth.bins),
  N = length(fit.inds$fit),
  tstep = 300,
  inds = start.inds,
  widths = depth.bins$halfwidth * 2,
  pi_priors = rbind(pi1.prior, pi2.prior),
  lambda_priors = rbind(lambda1.prior, lambda2.prior, lambda3.prior),
  stage_duration_priors = rbind(T1.prior.params, T2.prior.params),
  dive_start_priors = matrix(c(-300, 300), nrow = length(fit.inds$fit), 
                             ncol = 2, byrow = TRUE),
  dive_end_priors = do.call(rbind, 
                            lapply(dives.obs.list[fit.inds$fit], 
                                   function(d) d$times[length(d$times)] + 
                                     c(-300, 300))),
  delta = 1e-10
)

dat = list(
  depths = depths,
  times = times
)

inits = list(
  pi = c(params$beta[1], .5, params$beta[2]),
  lambda = params$lambda,
  xi = as.matrix(times.stages.est * 60),
  T = cbind(0, 
            t(apply(times.stages.est, 1, function(x) c(x[1], x[1]+x[2]))* 60), 
            apply(consts$dive_end_priors, 1, mean))
)

inits$logit_pi = qlogis(inits$pi)
inits$log_lambda = log(inits$lambda)

model = nimbleModel(
  code = modelCode,
  constants = consts,
  data = dat,
  inits = inits,
  name = 'ctdsDives'
)

cmodel = compileNimble(model, projectName = 'ctdsDives', resetFunctions = TRUE)

cmodel$calculate()



#
# begin to configure sampler
#

cfg_mcmc = configureMCMC(model)

cfg_mcmc$removeSampler(c('logit_pi', 'log_lambda', 'xi', 'T'))

cfg_mcmc$addMonitors(c('pi', 'lambda'))


#
# semi-warm start for model parameters, and use covariances for RW proposals
#


# stage 1 and 3 parameters
o_joint = lapply(c(1,3), function(s) {
  
  # approximate posterior mode
  o = optim(par = c(0,0), function(theta) {
    cmodel$logit_pi[s] = theta[1]
    cmodel$log_lambda[s] = theta[2]
    cmodel$calculate()
  }, control = list(fnscale = -1), hessian = TRUE)
  
  # add sampler
  cfg_mcmc$addSampler(
    target = paste(c('logit_pi[', 'log_lambda['), s, ']', sep = ''), 
    type = 'RW_block', 
    control = list(propCov = solve(-o$hessian))
  )
  
  o
})

# stage 2 parameter
o = optim(par = 0, function(theta) {
  cmodel$log_lambda[2] = theta
  cmodel$calculate()
}, control = list(fnscale = -1), hessian = TRUE)

cmodel$pi
cmodel$lambda

# stage 2 sampler
cfg_mcmc$addSampler(
  target = 'log_lambda[2]', 
  type = 'RW',
  control = list(scale = as.numeric(sqrt(solve(-o$hessian))))
)

for(i in 1:consts$N) {
  
  print(i)
  
  t1.tgt = paste('T[', i, ',1]', sep = '')
  t4.tgt = paste('T[', i, ',4]', sep = '')
  xi1.tgt = paste('xi[', i, ',1]', sep = '')
  xi2.tgt = paste('xi[', i, ',2]', sep = '')
    
  # jointly optimize dive-specific random effects
  o = optim(par = c(inits$T[i,c(1,4)], inits$xi[i,]), function(u) {
    cmodel[[t1.tgt]] = u[1]
    cmodel[[t4.tgt]] = u[2]
    cmodel[[xi1.tgt]] = u[3]
    cmodel[[xi2.tgt]] = u[4]
    cmodel$calculate()
  }, control = list(fnscale = -1), method = 'BFGS', hessian = FALSE)
  
  #
  # get marginal sd's for initial proposal distributions
  #
  
  tgt.ranges = list()
  tgt.ranges[[t1.tgt]] = consts$dive_start_priors[i,]
  tgt.ranges[[t4.tgt]] = consts$dive_end_priors[i,]
  
  for(tgt in c(t1.tgt, t4.tgt)) {
    # optimize parameter
    o = optim(par = cmodel[[tgt]], function(u) {
      cmodel[[tgt]] = u
      cmodel$calculate()
    }, control = list(fnscale = -1), hessian = TRUE, method = 'Brent',
    lower = tgt.ranges[[tgt]][1], upper = tgt.ranges[[tgt]][2])
    
    # add sampler
    cfg_mcmc$addSampler(
      target = tgt, 
      type = 'RW', 
      control = list(scale = as.numeric(sqrt(solve(abs(o$hessian))))))
  }
  
  for(tgt in c(xi1.tgt, xi2.tgt)) {
    # optimize parameter
    o = optim(par = cmodel[[tgt]], function(u) {
      cmodel[[tgt]] = u
      cmodel$calculate()
    }, control = list(fnscale = -1), hessian = TRUE)
    
    # add sampler
    cfg_mcmc$addSampler(
      target = tgt, 
      type = 'RW', 
      control = list(scale = as.numeric(sqrt(solve(abs(o$hessian))))))
  }

}


#
# construct sampler
#


model_mcmc = buildMCMC(cfg_mcmc)

cmcmc = compileNimble(model_mcmc, resetFunctions = TRUE)


#
# posterior samples
#

niter = 1e4

ncheckpoints = 20

for(i in 1:ncheckpoints) {
  cmcmc$run(niter = ceiling(niter/ncheckpoints), reset = FALSE, 
            progressBar = TRUE)
  samples = as.matrix(cmcmc$mvSamples)
  save.time = Sys.time()
  save(samples, save.time,
       file = 'sketches/2020-06-04_nimble_implementation/methods/samples.RData')
}


# library(coda)
# 
# plot(mcmc((samples[,'pi[1]'])))
# plot(mcmc((samples[,'pi[3]'])))
# plot(mcmc((samples[,'lambda[1]'])))
# plot(mcmc((samples[,'lambda[2]'])))
# plot(mcmc((samples[,'lambda[3]'])))
# 
# plot(mcmc((samples[,'xi[3, 1]'])))
# plot(mcmc((samples[,'xi[3, 2]'])))
# plot(mcmc((samples[,'T[3, 1]'])))
# plot(mcmc((samples[,'T[3, 4]'])))
# 
# 
# model$getDependencies('pi[1]')