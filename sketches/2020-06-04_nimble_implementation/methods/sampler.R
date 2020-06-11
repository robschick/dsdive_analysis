
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
            apply(consts$dive_end_priors, 1, mean))#,
  # expm = array(0, dim = c(3, consts$M, consts$M)),
  # evecs = array(0, dim = c(3, consts$M, consts$M)),
  # evals = matrix(0, nrow = 3, ncol = consts$M),
  # d = matrix(0, nrow = 3, ncol = consts$M),
  # dInv = matrix(0, nrow = 3, ncol = consts$M)
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
#
#

for(i in 1:1) {
  print(i)
  o = optim(par = c(inits$T[i,c(1,4)], inits$xi[i,]), function(u) {
    cmodel[[paste('T[', i, ',1]', sep = '')]] = u[1]
    cmodel[[paste('T[', i, ',4]', sep = '')]] = u[2]
    cmodel[[paste('xi[', i, ',1]', sep = '')]] = u[3]
    cmodel[[paste('xi[', i, ',2]', sep = '')]] = u[4]
    cmodel$calculate()
  }, control = list(fnscale = -1))
  print(o$par)
}


#
# construct sampler
#

cfg_mcmc = configureMCMC(model)

# find initial covariance matrix for block proposals
# o = optim(par = c(jjff))


for(s in c(1,3)) {
  o = optim(par = c(0,0), function(theta) {
    cmodel$logit_pi[s] = theta[1]
    cmodel$log_lambda[s] = theta[2]
    cmodel$calculate()
  }, control = list(fnscale = -1), hessian = FALSE)
}

cmodel$pi
cmodel$lambda

library(ggplot2)

for(i in fit.inds$fit) {
  svec = stagevec(length.out = length(dives.obs.list[[i]]$depths), 
                  breaks = which(diff(findInterval(dives.obs.list[[i]]$times, 
                                                   cmodel$T[i,2:3]))==1))
  print(plot(x = dives.obs.list[[i]], depth.bins = depth.bins, 
             stages = svec, errorbars = TRUE) + ggtitle(i))
}

curve(sapply(x, function(x) {
  cmodel$xi[1,1] = x
  cmodel$calculate()
}), from = 0, to = 1200)


curve(sapply(x, function(x) {
  cmodel[['logit_pi[1]']] = logit(x)
  cmodel$calculate()
}))

cfg_mcmc$removeSampler(c('logit_pi[1]', 'logit_pi[3]', 
                         'log_lambda[1]', 'log_lambda[3]',
                         'xi'))

cfg_mcmc$addSampler(target = c('logit_pi[1]', 'log_lambda[1]'), 
                    type = 'RW_block')

cfg_mcmc$addSampler(target = c('logit_pi[3]', 'log_lambda[3]'), 
                    type = 'RW_block')

cfg_mcmc$addMonitors(c('pi', 'lambda'))

xi_nodes = model$getNodeNames()[grep('xi.*', model$getNodeNames())]
for(x in xi_nodes) {
  cfg_mcmc$addSampler(target = x, type = 'slice')
}


cfg_mcmc
model_mcmc = buildMCMC(cfg_mcmc)


cmcmc = compileNimble(model_mcmc, resetFunctions = TRUE)


#
# posterior samples
#

samples = runMCMC(cmcmc, niter = 1e2)


library(coda)

plot(mcmc((samples[,'pi[1]'])))
(mcmc(plogis(samples[(1:50),380:382])))


effectiveSize(mcmc(samples[-(1:50),]))



model$getDependencies(model$getNodeNames()[264])

samplerFn = cmcmc$samplerFunctions[[96]][[1]]
samplerFn$memberData(name = 'control')

edit(sampler_RW_block)
