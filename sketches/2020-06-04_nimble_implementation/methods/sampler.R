
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
  dive_end_priors = do.call(rbind, lapply(dives.obs.list[fit.inds$fit], 
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
  T = cbind(0, t(apply(times.stages.est, 1, 
                       function(x) c(x[1], x[1]+x[2]))* 60) , 
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
# warm start for model parameters; get initial proposal covariances
#

for(i in 1:nrow(inits$T)) {
  print(i)
  o = optim(par = c(inits$T[i,c(1,4)], inits$xi[i,]), function(u) {
    cmodel[[paste('T[', i, ',1]', sep = '')]] = u[1]
    cmodel[[paste('T[', i, ',4]', sep = '')]] = u[2]
    cmodel[[paste('xi[', i, ',1]', sep = '')]] = u[3]
    cmodel[[paste('xi[', i, ',2]', sep = '')]] = u[4]
    cmodel$calculate()
  }, control = list(fnscale = -1), method = 'BFGS')
  print(o$par)
}


o_init = lapply(c(1,3), function(s) {
  optim(par = c(0,0), function(theta) {
    cmodel$logit_pi[s] = theta[1]
    cmodel$log_lambda[s] = theta[2]
    cmodel$calculate()
  }, control = list(fnscale = -1), hessian = TRUE)
})

o_last = optim(par = 0, function(theta) {
  cmodel$log_lambda[2] = theta
  cmodel$calculate()
}, control = list(fnscale = -1), hessian = TRUE)

cmodel$pi
cmodel$lambda


#
# construct sampler
#

cfg_mcmc = configureMCMC(model)

cfg_mcmc$removeSampler(c('logit_pi', 'log_lambda'))

cfg_mcmc$addSampler(target = c('logit_pi[1]', 'log_lambda[1]'), 
                    type = 'RW_block', 
                    control = list(propCov = solve(-o_init[[1]]$hessian)))

cfg_mcmc$addSampler(target = c('logit_pi[3]', 'log_lambda[3]'), 
                    type = 'RW_block', 
                    control = list(propCov = solve(-o_init[[2]]$hessian)))

cfg_mcmc$addSampler(target = 'log_lambda[2]',
                    type = 'RW',
                    control = list(scale = as.numeric(sqrt(solve(-o_last$hessian)))))

cfg_mcmc$addMonitors(c('pi', 'lambda'))


cfg_mcmc
model_mcmc = buildMCMC(cfg_mcmc)


cmcmc = compileNimble(model_mcmc, resetFunctions = TRUE)


#
# posterior samples
#

samples = runMCMC(cmcmc, niter = 1e2)


library(coda)

plot(mcmc((samples[,'pi[1]'])))
plot(mcmc((samples[,'pi[3]'])))
plot(mcmc((samples[,'lambda[1]'])))
plot(mcmc((samples[,'lambda[2]'])))
plot(mcmc((samples[,'lambda[3]'])))

plot(mcmc((samples[,'xi[3, 1]'])))
plot(mcmc((samples[,'xi[3, 2]'])))
plot(mcmc((samples[,'T[3, 1]'])))
plot(mcmc((samples[,'T[3, 4]'])))


model$getDependencies('pi[1]')