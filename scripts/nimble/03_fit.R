library(nimble)

# output directory
out.dir = file.path('output', 'multitag')
dir.create(out.dir, recursive = TRUE)

# load data
nim_pkg = readRDS(file.path(out.dir, 'nim_pkg.rds'))

# add timestep and approximation, for likelihood computation
nim_pkg$consts$tstep = 300
nim_pkg$consts$delta = 1e-10

# load model and distributions
source(file.path('scripts', 'nimble', 'methods', 'distributions.R'))
source(file.path('scripts', 'nimble', 'methods', 'model_multitag.R'))


#
# begin to configure sampler
#

dive_model = nimbleModel(code = modelCode, constants = nim_pkg$consts, 
                         data = nim_pkg$data, inits = nim_pkg$inits, 
                         name = 'dives')

dive_model$initializeInfo()

cmodel = compileNimble(dive_model)

if(is.infinite(cmodel$calculate())) {
  cmodel.names = cmodel$getNodeNames()
  ll.nodes = sapply(cmodel.names, function(x) cmodel$calculate(x))
  names.inf = cmodel.names[which(is.infinite(ll.nodes))]
  print(names.inf)
}

cfg_mcmc = configureMCMC(cmodel, print = TRUE)


##
## semi-warm start for model parameters, and use covariances for RW proposals
##

#
# movement parameter samplers
#

cfg_mcmc$removeSamplers(c('logit_pi', 'log_lambda'))

o_joint = lapply(1:nim_pkg$consts$N_tags, function(tagId) {

  # stage 1 and 3 parameters
  lapply(c(1,3), function(s) {

    # status update
    message(paste('tag:', tagId, 'stage:', s, sep = ' '))
    
    # specify target nodes, and extract dependencies
    tgt = paste(c('logit_pi[', 'log_lambda['), tagId, ', ', s, ']', sep = '')
    deps = cmodel$getDependencies(tgt)
  
    
    # find posterior mode and hessian
    o = optim(par = c(cmodel$logit_pi[tagId, s], cmodel$log_lambda[tagId, s]), 
              fn = function(theta) {
      cmodel$logit_pi[tagId, s] = theta[1]
      cmodel$log_lambda[tagId, s] = theta[2]
      cmodel$calculate(deps)
    }, method = 'BFGS', control = list(fnscale = -1), hessian = TRUE)
    
    cov = solve(-o$hessian)
    if(any(eigen(cov)$values < 0)) {
      cov = diag(2)
    }

    # add sampler
    cfg_mcmc$addSampler(target = tgt, type = 'RW_block', silent = TRUE,
                        control = list(propCov = cov)
    )
  })

  message(paste('tag:', tagId, 'stage:', 2, sep = ' '))

  # stage 2 parameter
  tgt = paste('log_lambda[', tagId, ', ', 2, ']', sep = '')
  deps = cmodel$getDependencies(tgt)
  
  o = optim(par = cmodel$log_lambda[tagId, 2], 
            fn = function(theta) {
    cmodel$log_lambda[tagId, 2] = theta
    cmodel$calculate(deps)
  }, method = 'BFGS', control = list(fnscale = -1), hessian = TRUE)
  
  cfg_mcmc$addSampler(target = tgt, type = 'RW', silent = TRUE,
                      control = list(scale = as.numeric(sqrt(-1/o$hessian)))
  )
  
  o
  
})


#
# dive endpoint samplers
#

cfg_mcmc$removeSamplers('endpoints')

for(i in 1:nim_pkg$consts$N_endpoints) {
  
  tgt = paste('endpoints[', i, ']', sep = '')
  deps = cmodel$getDependencies(tgt)
  
  o = optim(par = cmodel[[tgt]], fn = function(theta) {
    cmodel[[tgt]] = theta
    cmodel$calculate(deps)
  }, method = 'Brent', control = list(fnscale = -1), hessian = TRUE,
  lower = nim_pkg$consts$endpoint_priors[i, 't_lwr'],
  upper = nim_pkg$consts$endpoint_priors[i, 't_upr'])
  
  sd = as.numeric(sqrt(-1/o$hessian))
  sd = ifelse(is.na(sd), 100, sd)
  
  cfg_mcmc$addSampler(target = tgt, type = 'RW', silent = TRUE,
                      control = list(scale = sd))
  
}


#
# stage duration samplers
#

cfg_mcmc$removeSamplers('log_xi')

for(i in 1:nim_pkg$consts$N_dives) {
  
  if(i %% 50 == 0) {
    message(paste('dive:', i))
  }
  
  tgt = paste('log_xi[', i, ', 1:2]', sep = '')
  deps = cmodel$getDependencies(tgt)
  
  o = optim(par = cmodel[[tgt]], fn = function(theta) {
    cmodel[[tgt]] = theta
    cmodel$calculate(deps)
  }, method = 'BFGS', control = list(fnscale = -1), hessian = TRUE)
  
  cov = solve(-o$hessian)
  if(any(eigen(cov)$values < 0)) {
    cov = diag(2)
  }
  
  cfg_mcmc$addSampler(target = tgt, type = 'RW_block', silent = TRUE, 
                      control = list(propCov = cov))
  
}


#
# construct sampler
#

model_mcmc = buildMCMC(cfg_mcmc)

cmcmc = compileNimble(model_mcmc, projectName = 'ctdsDives')


#
# posterior samples
#


niter = 1e4

ncheckpoints = 20

for(i in 1:ncheckpoints) {
  chunk_iter = ceiling(niter/ncheckpoints)
  cmcmc$run(niter = chunk_iter, reset = FALSE, progressBar = TRUE)
  samples = as.matrix(cmcmc$mvSamples)
  save.time = Sys.time()
  save(samples, save.time, file = file.path(out.dir, 'samples.RData'))
}


# library(coda)
# 
# plot(mcmc((plogis(samples[,'logit_pi[1, 1]']))))
# plot(mcmc((plogis(samples[,'logit_pi[1, 3]']))))
# 
# plot(mcmc((exp(samples[,'log_lambda[6, 1]']))))
# plot(mcmc((exp(samples[,'log_lambda[6, 2]']))))
# plot(mcmc((exp(samples[,'log_lambda[6, 3]']))))
# 
# plot(mcmc((exp(samples[,'log_xi[300, 2]']))))
# 
# plot(mcmc(((samples[,'endpoints[1]']))))
