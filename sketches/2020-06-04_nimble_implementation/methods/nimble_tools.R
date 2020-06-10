library(nimble)


#
# support for CTMC model across depth bins
#

# nimbleList for storing entries of a tridiagonal matrix
tridiagonalEntries = nimbleList(
  diag = double(1),
  dsuper = double(1),
  dsub = double(1)
)

# build components for infinitesimal generator matrix for depth bin transitions
buildInfinitesimalGeneratorEntries = nimbleFunction(
  run = function(pi = double(0), lambda = double(0), M = integer(0),
                 stage = integer(0), widths = double(1)) {
    # Parameters:
    #  pi - probability of making a downward descent
    #  lambda - speed of moving through bin
    #  M - total number of depth bins
    #  stage - build generator matrix for "stage"
    #  widths - vector of depth bin widths
    
    returnType(tridiagonalEntries())
    
    entries <- tridiagonalEntries$new(diag = numeric(M, init = FALSE), 
                                      dsub = numeric(M-1, init = FALSE),
                                      dsuper = numeric(M-1, init = FALSE))
    
    # only allow (downward) transitions from shallowest bin during stage 1
    if(stage == 1) {
      rate <- lambda / widths[1]
      entries$diag[1] <- -rate
      entries$dsuper[1] <- rate
    }
    
    # only allow surfacing from second bin to occur in stage 3
    rate <- lambda / widths[2]
    if(stage == 1 | stage == 2) {
      entries$diag[2] <- -rate
      entries$dsuper[2] <- rate
    } else {
      entries$diag[2] <- rate
      entries$dsuper[2] <- rate * pi
      entries$dsub[1] <- rate * (1-pi)
    }
    
    # intermediate bins may always transition up or down
    for(i in 3:(M-1)) {
      rate <- lambda / widths[i]
      entries$diag[i] <- -rate
      entries$dsuper[i] <- rate * pi
      entries$dsub[i-1] <- rate * (1-pi)
    }
    
    # deepest bin can only transition to next shallowest depth
    rate <- lambda / widths[M]
    entries$diag[M] <- -rate
    entries$dsub[M-1] <- rate
    
    return(entries)
  }
)

# build infinitesimal generator matrix for depth bin transitions
buildInfinitesimalGenerator = nimbleFunction(
  run = function(pi = double(0), lambda = double(0), M = integer(0),
                 stage = integer(0), widths = double(1)) {
    # Parameters:
    #  pi - probability of making a downward descent
    #  lambda - speed of moving through bin
    #  M - total number of depth bins
    #  stage - build generator matrix for "stage"
    #  widths - vector of depth bin widths

    returnType(double(2))
    
    # compute matrix's tridiagonal entries
    entries <- buildInfinitesimalGeneratorEntries(pi = pi, lambda = lambda, 
                                                  M = M, stage = stage, 
                                                  widths = widths)

    # initialize matrix
    A <- matrix(0, nrow = M, ncol = M)
    
    #
    # populate matrix
    #
    
    for(i in 1:(M-1)) {
      A[i,i] <- entries$diag[i]
      A[i,i+1] <- entries$dsuper[i]
      A[i+1,i] <- entries$dsub[i]
    }
    
    A[M,M] <- entries$diag[M]
    
    return(A)
  }
)


#
# support for matrix exponentials
#

# compile C++ implementation of LAPACK-based decomposition of generators
system(paste(
  'R CMD COMPILE',
  paste('CXXFLAGS=', '"',
        paste('-I', 
              file.path(find.package(c('RcppEigen', 'Rcpp')), 'include'), 
              sep = '', collapse = ' '), 
        '"', sep = ''),
  'sketches/2020-06-04_nimble_implementation/exp_symtools.cpp',
  sep = ' '
))

# link NIMBLE to external C++ (decomposition)
decomposeGeneratorCpp = nimbleExternalCall(
  prototype = function(diag = double(1), dsuper = double(1), dsub = double(1), 
                       N = integer(0), expm = double(2), evecs = double(2), 
                       d = double(1), dInv = double(1), delta = double(0),
                       t = double(0)){},  
  returnType = void(), 
  Cfun = 'expm_cpp', 
  headerFile = file.path(getwd(), 'sketches/2020-06-04_nimble_implementation', 
                         'exp_symtools.h'), 
  oFile = file.path(getwd(), 'sketches/2020-06-04_nimble_implementation', 
                    'exp_symtools.o')
)

# nimbleList for eigendecomposition of matrix exponential
eigenExpm = nimbleList(
  expm = double(2),
  evecs = double(2),
  evals = double(1),
  d = double(1),
  dInv = double(1),
  t = double(0),
  N = integer(0)
)

# construct eigendecomposition of matrix exponential
decomposeGenerator = nimbleFunction(
  run = function(Aentries = tridiagonalEntries(), delta = double(0), 
                 t = double(0), N = integer(0)) {

    returnType(eigenExpm())
    
    expm <- matrix(0, nrow = N, ncol = N)
    evecs <- matrix(0, nrow = N, ncol = N)
    d <- numeric(0, length = N)
    dInv <- numeric(0, length = N)
    
    evals <- Aentries$diag

    decomposeGeneratorCpp(diag = evals, dsuper = Aentries$dsuper, 
                          dsub = Aentries$dsub, N = N, expm = expm, 
                          evecs = evecs, d = d, dInv = dInv, delta = delta, 
                          t = t)
    
    res <- eigenExpm$new(
      expm = expm,
      evecs = evecs,
      evals = evals,
      d = d,
      dInv = dInv,
      t = t,
      N = N
    )
    
    return(res)
  }
)

# compile nimble code, so we can use it
cdecomposeGenerator = compileNimble(decomposeGenerator)
  
# link NIMBLE to external C++ (multiplication)
expmAtvCpp = nimbleExternalCall(
  prototype = function(evecs = double(2), evals = double(1), M = integer(0), 
                       v = double(1), d = double(1), dInv = double(1), 
                       t = double(0), x = double(1), 
                       preMultiply = logical(0)){},
  returnType = void(), 
  Cfun = 'expmAtv_cpp', 
  headerFile = file.path(getwd(), 'sketches/2020-06-04_nimble_implementation', 
                         'exp_symtools.h'), 
  oFile = file.path(getwd(), 'sketches/2020-06-04_nimble_implementation', 
                    'exp_symtools.o')
)

# use generator matrix decomposition to compute x = exp(At)*v
expmAtv = nimbleFunction(
  run = function(Adecomp = eigenExpm(0), v = double(1), t = double(0), 
                 preMultiply = logical(0)) {
    
    returnType(double(1))
    
    x <- numeric(0, length = Adecomp$N)
    
    if(t == Adecomp$t) {
      
      N <- Adecomp$N
      
      if(preMultiply == TRUE) {
        x[1:N] <- t(v[1:N]) %*% Adecomp$expm[1:N,1:N]
      } else {
        x[1:N] <- Adecomp$expm[1:N,1:N] %*% v[1:N]
      }
    } else {
      expmAtvCpp(evecs = Adecomp$evecs, evals = Adecomp$evals, M = Adecomp$N, 
                 v = v, d = Adecomp$d, dInv = Adecomp$dInv, t = t, x = x, 
                 preMultiply = preMultiply)
    }
    
    return(x)
  }
)

# compile nimble code, so we can use it
cexpmAtv = compileNimble(expmAtv)


# #
# # Test some code
# #
# 
# Aentries = buildInfinitesimalGeneratorEntries(pi = .95, lambda = 1.4,
#                                               M = nrow(depth.bins), stage = 1,
#                                               widths = 2 * depth.bins$halfwidth)
# 
# Adecomp = cdecomposeGenerator(Aentries = Aentries, delta = 1e-15, t = 300,
#                               N = nrow(depth.bins))
# 
# t0 = 300
# 
# x0 = c(1, rep(0, 15))
# 
# prob = cexpmAtv(Adecomp = Adecomp, v = x0, t = t0, preMultiply = TRUE)
# 
# A = buildInfinitesimalGenerator(pi = .95, lambda = 1.4,
#                                 M = nrow(depth.bins), stage = 1,
#                                 widths = 2 * depth.bins$halfwidth)
# 
# max(abs(x0 %*% Matrix::expm(A * t0) - prob))


#
# prior densities
#

# density for logit(X) where X ~ Beta(shape1, shape2)
dlogitBeta = nimbleFunction(
  run = function(x = double(0), shape1 = double(0), shape2 = double(0), 
                 log = logical(0, default = 0)) {
    
    returnType(double(0))
    
    res <- dbeta(x = ilogit(x), shape1 = shape1, shape2 = shape2, log = TRUE) + 
      x - 2 * log(exp(x) + 1)
    
    if(log) { return(res) } else { return(exp(res)) }
  }
)

# density for log(X) where X ~ Gamma(shape, rate)
dlogGamma = nimbleFunction(
  run = function(x = double(0), shape = double(0), rate = double(0), 
                 log = logical(0, default = 0)) {
    
    returnType(double(0))
    
    res <- dgamma(x = exp(x), shape = shape, rate = rate, log = TRUE) + x
    
    if(log) { return(res) } else { return(exp(res)) }
  }
)


#
# likelihood function
#

# transition kernel for a single dive transition
ddiveTx = nimbleFunction(
  run = function(x = integer(0), x_prev = integer(0), Ad1 = eigenExpm(),
                 Ad2 = eigenExpm(), Ad3 = eigenExpm(),
                 M = integer(0), t_prev = double(0), t = double(0),
                 T = double(1), log = logical(0, default = 0)) {
    # Parameters:
    #   x - depth bin transitioned to
    #   x_prev - depth bin transitioned from
    #   Ad1 - decomposition object for stage 1 infinitesimal generator matrix
    #   Ad2 - decomposition object for stage 2 infinitesimal generator matrix
    #   Ad3 - decomposition object for stage 3 infinitesimal generator matrix
    #   M - number of depth bins
    #   t_prev - time at which x_last was observed
    #   t - time at which x was observed
    #   T - vector of dive-specific random effects;  c(T0, T1, T2, T3) s.t. 
    #       T0 is dive start time, T1 and T2 are respecitively end of stage 1/2,
    #       and T3 is dive end time.  note that the indexing implies T[1] is 
    #       start of dive, T[2] is start of stage 2, T[3] is start of stage 3, 
    #       and T[4] is start of surface period (i.e., end of dive).
    #   log - If TRUE, then log-likelihood is returned
    
    returnType(double(0))
    
    
    if(t < T[1]) {
      # transition occurs before dive begins
      prob <- 1
    } else if( t_prev > T[4]) {
      # transition occurs after dive ends
      prob <- 1
    } else {
      # transition occurs during dive
      
      # (x_prev, t_prev) occur before dive, so align transition with dive start
      if(t_prev < T[1]) {
        x_prev_local <- 1
        t_prev_local <- T[1]
      } else {
        x_prev_local <- x_prev
        t_prev_local <- t_prev
      }
      
      # (x, t) occur after dive, so align transition with dive end
      if(t > T[4]) {
        x_local <- 1
        t_local <- T[4]
      } else {
        x_local <- x
        t_local <- t
      }
      
      # initialize transition distribution
      u <- numeric(0, length = M)
      u[x_prev_local] <- 1
      
      # determine dive stage at x_prev_local
      if(t_prev_local >= T[3]) {
        stage_prev <- 3
      } else if(t_prev_local >= T[2]) {
        stage_prev <- 2
      } else {
        stage_prev <- 1
      }
      
      # time at which x_prev_local "entered" stage_prev
      t_stage <- t_prev_local
      
      # diffuse transition mass across stages
      diffusing <- TRUE
      for(s in stage_prev:3) {
        if(diffusing) {
          # transition ends after current stage
          if(t_local >= T[s+1]) {
            # time spent in stage s
            d_t <- T[s+1] - t_stage
            # time at which x_prev "enters" next stage
            t_stage <- T[s+1]
          }
          # transition ends in current stage
          else {
            # time spent in stage s
            d_t <- t_local - t_stage
            # signal end of loop
            diffusing <- FALSE
          }
          
          # diffuse mass across stage s
          if(s == 1) {
            u[1:M] <- expmAtv(Adecomp = Ad1, v = u[1:M], t = d_t, 
                              preMultiply = TRUE)[1:M]
          } else if(s == 2) {
            u[1:M] <- expmAtv(Adecomp = Ad2, v = u[1:M], t = d_t, 
                              preMultiply = TRUE)[1:M]
          } else if(s == 3) {
            u[1:M] <- expmAtv(Adecomp = Ad3, v = u[1:M], t = d_t, 
                              preMultiply = TRUE)[1:M]
          }
          
        }
      }
      
      # extract transition probability
      prob <- u[x_local]
      
    }
    
    if(log) { return(log(prob)) } else { return(prob) }
  }
)

registerDistributions(list(
  ddiveTx = list(
    BUGSdist = "ddiveTx(x_prev, Ad1, Ad2, Ad3, M, t_prev, t, T)",
    discrete = TRUE,
    range = c(0, Inf),
    types = c('value = integer()', 'x_prev = integer()', 'Ad1 = eigenExpm()',
              'Ad2 = eigenExpm()', 'Ad3 = eigenExpm()', 'M = integer()', 
              't_prev = double()', 't = double()', 'T = double()')
)))


# #
# # Test likelihood
# #
# 
# pi = c(.95, .5, .05)
# lambda = c(1.5, .3, 1)
# 
# Adecomps = lapply(1:3, function(s) {
#   Aentries = buildInfinitesimalGeneratorEntries(
#     pi = pi[s], lambda = lambda[s], M = 16, stage = s, 
#     widths = 2 * depth.bins$halfwidth)
#   
#   cdecomposeGenerator(Aentries = Aentries, delta = 1e-15, t = 300, N = 16)
# })
# 
# 
# cddiveTx = compileNimble(ddiveTx)
# 
# cddiveTx(x = 5, x_prev = 1, Ad1 = Adecomps[[1]], Ad2 = Adecomps[[2]], 
#          Ad3 = Adecomps[[3]], M = 16, t_prev = 0, t = 300, 
#          T = c(500, 600, 800, Inf), log = TRUE)




#
# dive model
#

modelCode = nimbleCode({

  #
  # priors
  #

  # directional preferences
  logit_pi[1] ~ dlogitBeta(shape1 = pi_priors[1,1], shape2 = pi_priors[1,2])
  logit_pi[3] ~ dlogitBeta(shape1 = pi_priors[2,1], shape2 = pi_priors[2,2])
  
  # transformed parameters, for sampling
  pi[1] <- ilogit(logit_pi[1])
  pi[3] <- ilogit(logit_pi[3])

  # speeds
  for(i in 1:3) {
    log_lambda[i] ~ dlogGamma(shape = lambda_priors[i,1], 
                              rate = lambda_priors[i,2])
    # transformed parameters, for sampling
    lambda[i] <- exp(log_lambda[i])
  }

  # dive-specific random effects
  for(i in 1:N) {
    # stage durations
    for(j in 1:2) {
      xi[i,j] ~ dgamma(shape = stage_duration_priors[j,1],
                       rate = stage_duration_priors[j,2])
    }
    # dive start and end times
    T[i,1] ~ dunif(dive_start_priors[i,1], dive_start_priors[i,2])
    T[i,4] ~ dunif(dive_end_priors[i,1], dive_end_priors[i,2])
    # stage 2 and 3 start times
    T[i,2] <- T[i,1] + xi[i,1]
    T[i,3] <- T[i,2] + xi[i,2]
  }

  #
  # CTMC generator matrix decompositions
  #
  
  A1entries <- buildInfinitesimalGeneratorEntries(
    pi = pi[1], lambda = lambda[1], M = M, stage = 1, widths = widths[1:M])
  Ad1 <- decomposeGenerator(Aentries = A1entries, delta = delta, t = tstep, 
                            N = M)
  
  A2entries <- buildInfinitesimalGeneratorEntries(
    pi = pi[2], lambda = lambda[2], M = M, stage = 2, widths = widths[1:M])
  Ad2 <- decomposeGenerator(Aentries = A1entries, delta = delta, t = tstep, 
                            N = M)
  
  A3entries <- buildInfinitesimalGeneratorEntries(
    pi = pi[3], lambda = lambda[3], M = M, stage = 3, widths = widths[1:M])
  Ad3 <- decomposeGenerator(Aentries = A1entries, delta = delta, t = tstep, 
                            N = M)
  
  #
  # likelihood
  #

  # loop over dives
  for(i in 1:N) {
    
    # initial depth observation
    depths[inds[i]] ~ ddiveTx(x_prev = 1, Ad1 = Ad1, Ad2 = Ad2, Ad3 = Ad3, 
                              M = M, t_prev = -Inf, t = times[inds[i]], 
                              T = T[i,1:4])
    # subsequent depth observations
    for(ind in (inds[i] + 1):(inds[i+1] - 1)) {
      depths[ind] ~ ddiveTx(x_prev = depths[ind-1], Ad1 = Ad1, Ad2 = Ad2, 
                            Ad3 = Ad3, M = M, t_prev = times[ind-1], 
                            t = times[ind], T = T[i,1:4])
    }
  }

})


#
# associate data and priors with model
#

consts = list(
  M = nrow(depth.bins),
  N = n.dives,
  tstep = 300,
  inds = start.inds,
  widths = depth.bins$halfwidth * 2,
  pi_priors = rbind(pi1.prior, pi2.prior),
  lambda_priors = rbind(lambda1.prior, lambda2.prior, lambda3.prior),
  stage_duration_priors = rbind(T1.prior.params, T2.prior.params),
  dive_start_priors = matrix(c(-300, 300), nrow = n.dives, ncol = 2, 
                             byrow = TRUE),
  dive_end_priors = do.call(rbind, lapply(dives.obs.list[fit.inds$fit], 
                                          function(d) d$times[length(d$times)] + 
                                            c(-300, 300)))
)

dat = list(
  depths = depths,
  times = times
)

inits = list(
  pi = c(params$beta[1], .5, params$beta[2]),
  lambda = params$lambda,
  xi = as.matrix(times.stages.est * 60),
  T = cbind(0, do.call(rbind, t.stages.list), 
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


#
# construct sampler
#

cfg_mcmc = configureMCMC(model)

# find initial covariance matrix for block proposals
# o = optim(par = c(jjff))

cfg_mcmc$removeSampler(c('logit_pi[1]', 'logit_pi[3]', 
                         'log_lambda[1]', 'log_lambda[3]',
                         'xi'))

cfg_mcmc$addSampler(target = c('logit_pi[1]', 'log_lambda[1]'), 
                    type = 'RW_block')

cfg_mcmc$addSampler(target = c('logit_pi[3]', 'log_lambda[3]'), 
                    type = 'RW_block')

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

plot(mcmc(exp(samples[-(1:50),1:3])))
plot(mcmc(plogis(samples[,c(4,6)])))


effectiveSize(mcmc(samples[-(1:50),]))

              