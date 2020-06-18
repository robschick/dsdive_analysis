

#
# begin to configure sampler
#

cfg_mcmc = configureMCMC(model)

cfg_mcmc$removeSampler(c('logit_pi', 'log_lambda'))

cfg_mcmc$addMonitors(c('pi', 'lambda'))


#
# semi-warm start for model parameters, and use covariances for RW proposals
#


# stage 1 and 3 parameters
o_joint = lapply(c(1,3,4,5), function(s) {

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


# for(i in 1:consts$N_endpoints) {
#   e.tgt = paste('E[', i, ']', sep = '')
#   
#   o = optim(par = inits$E[i], function(u) {
#     cmodel[[e.tgt]] = u
#     cmodel$calculate()
#   }, control = list(fnscale = -1), method = 'Brent', hessian = TRUE,
#   lower = consts$E_priors[i,'t_lwr'],
#   upper = consts$E_priors[i,'t_upr'])
#   
#   cfg_mcmc$addSampler(
#     target = e.tgt,
#     type = 'RW',
#     control = list(scale = as.numeric(sqrt(solve(-o$hessian))))
#   )
#   
# }





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
       file = 'sketches/2020-06-17_nimble_fulltag/methods/samples.RData')
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