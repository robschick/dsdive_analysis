
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
  
  for(s in 1:3) {
    expm_decomp[s,1:(2*M+3),1:M] <- buildAndDecomposeGenerator(
      pi = pi[s], lambda = lambda[s], M = M, stage = s, widths = widths[1:M], 
      delta = delta, t = tstep)
  }
  
  #
  # likelihood
  #
  
  # loop over dives
  for(i in 1:N) {
    
    depths[inds[i]:(inds[i+1]-1)] ~ ddive(
      times = times[inds[i]:(inds[i+1]-1)], N = inds[i+1] - inds[i] + 1, 
      expm = expm_decomp[1:3,1:M,1:M],
      evecs = expm_decomp[1:3,(M+1):(2*M),1:M],
      evals = expm_decomp[1:3,2*M+1,1:M],
      d = expm_decomp[1:3,2*M+2,1:M],
      dInv = expm_decomp[1:3,2*M+3,1:M],
      tstep = tstep, M = M, T = T[i,1:4])
  }
  
})

