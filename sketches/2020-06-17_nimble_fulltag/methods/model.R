
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
  logit_pi[4] ~ dlogitBeta(shape1 = pi_priors[3,1], shape2 = pi_priors[3,2])
  logit_pi[5] ~ dlogitBeta(shape1 = pi_priors[4,1], shape2 = pi_priors[4,2])
  
  # transformed parameters, for sampling
  pi[1] <- ilogit(logit_pi[1])
  pi[3] <- ilogit(logit_pi[3])
  pi[4] <- ilogit(logit_pi[4])
  pi[5] <- ilogit(logit_pi[5])
  
  # speeds
  for(i in 1:5) {
    log_lambda[i] ~ dlogGamma(shape = lambda_priors[i,1], 
                              rate = lambda_priors[i,2])
    
    # transformed parameters, for sampling
    lambda[i] <- exp(log_lambda[i])
  }
  
  # dive endpoint random effects
  for(i in 1:N_endpoints) {
    E[i] ~ dunif(E_priors[i,2], E_priors[i,3])
  }
  
  #
  # CTMC generator matrix decompositions
  #
  
  # deep dive generators
  for(s in 1:3) {
    expm_decomp[s,1:(2*M+3),1:M] <- buildAndDecomposeGenerator(
      pi = pi[s], lambda = lambda[s], M = M, stage = s, widths = widths[1:M], 
      delta = delta, t = tstep)
  }
  
  # shallow dive descent generator
  expm_decomp[4,1:(2*M+3),1:M] <- buildAndDecomposeGenerator(
    pi = pi[4], lambda = lambda[4], M = M, stage = 1, widths = widths[1:M], 
    delta = delta, t = tstep)
  
  # shallow dive ascent generator
  expm_decomp[5,1:(2*M+3),1:M] <- buildAndDecomposeGenerator(
    pi = pi[5], lambda = lambda[5], M = M, stage = 3, widths = widths[1:M], 
    delta = delta, t = tstep)
  
  #
  # likelihood, and additional priors
  #

  # deep dive stage transition times, and likelihood
  for(i in 1:N_ranges_deep) {

    # stage duration priors
    for(j in 1:2) {
      xi[dive_priors_deep[i,1],j] ~ dgamma(
        shape = stage_duration_priors[j,1],
        rate = stage_duration_priors[j,2])
    }

    # dive start and end times
    T[dive_priors_deep[i,1],1] <- E[E_map[dive_priors_deep[i,1],1]]
    T[dive_priors_deep[i,1],4] <- E[E_map[dive_priors_deep[i,1],2]]

    # stage transition times
    T[dive_priors_deep[i,1],2] <- T[dive_priors_deep[i,1],1] +
                                  xi[dive_priors_deep[i,1],1]
    T[dive_priors_deep[i,1],3] <- T[dive_priors_deep[i,1],2] +
                                  xi[dive_priors_deep[i,1],2]

    # stage 3 duration
    xi[dive_priors_deep[i,1],3] <- T[dive_priors_deep[i,1],4] -
                                   T[dive_priors_deep[i,1],3]

    # likelihood
    depths[dive_priors_deep[i,6]:(dive_priors_deep[i,6] + dive_priors_deep[i,5] - 1)] ~ ddive(
      times = times[dive_priors_deep[i,6]:(dive_priors_deep[i,6] + dive_priors_deep[i,5] - 1)],
      N = dive_priors_deep[i,3] - dive_priors_deep[i,2] + 1,
      expm = expm_decomp[1:3,1:M,1:M],
      evecs = expm_decomp[1:3,(M+1):(2*M),1:M],
      evals = expm_decomp[1:3,2*M+1,1:M],
      d = expm_decomp[1:3,2*M+2,1:M],
      dInv = expm_decomp[1:3,2*M+3,1:M],
      tstep = tstep, M = M, T = T[dive_priors_deep[i,1],1:4],
      include = training_dive[dive_priors_deep[i,1]]
    )
  }

  # shallow dive stage transition times, and likelihood
  for(i in 1:N_ranges_shallow) {

    # stage duration priors
    xi[dive_priors_shallow[i,1],1] ~ dgamma(
      shape = stage_duration_priors[3,1],
      rate = stage_duration_priors[3,2])

    # dive start and end times
    T[dive_priors_shallow[i,1],1] <- E[E_map[dive_priors_shallow[i,1],1]]
    T[dive_priors_shallow[i,1],4] <- E[E_map[dive_priors_shallow[i,1],2]]

    # stage transition times
    T[dive_priors_shallow[i,1],2] <- T[dive_priors_shallow[i,1],1] +
                                     xi[dive_priors_shallow[i,1],1]

    # stage 2 duration
    xi[dive_priors_shallow[i,1],2] <- T[dive_priors_shallow[i,1],4] -
                                      T[dive_priors_shallow[i,1],2]

    # likelihood
    depths[dive_priors_shallow[i,6]:(dive_priors_shallow[i,6] + dive_priors_shallow[i,5] - 1)] ~ ddiveShallow(
      times = times[dive_priors_shallow[i,6]:(dive_priors_shallow[i,6] + dive_priors_shallow[i,5] - 1)],
      N = dive_priors_shallow[i,3] - dive_priors_shallow[i,2] + 1,
      expm = expm_decomp[4:5,1:M,1:M],
      evecs = expm_decomp[4:5,(M+1):(2*M),1:M],
      evals = expm_decomp[4:5,2*M+1,1:M],
      d = expm_decomp[4:5,2*M+2,1:M],
      dInv = expm_decomp[4:5,2*M+3,1:M],
      tstep = tstep, M = M, T = T[dive_priors_shallow[i,1],1:4],
      include = training_dive[dive_priors_shallow[i,1]]
    )
  }
  
})
