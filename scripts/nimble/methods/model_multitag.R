
#
# dive model
#

modelCode = nimbleCode({
  
  #
  # priors
  #
  
  for(tagInd in 1:N_tags) {
    
    # directional preferences
    logit_pi[tagInd, 1] ~ dlogitBeta(shape1 = pi_priors[1, 1], 
                                     shape2 = pi_priors[1, 2])
    logit_pi[tagInd, 3] ~ dlogitBeta(shape1 = pi_priors[2, 1], 
                                     shape2 = pi_priors[2, 2])
    
    # transformed parameters, for sampling
    pi[tagInd, 1] <- ilogit(logit_pi[tagInd, 1])
    pi[tagInd, 3] <- ilogit(logit_pi[tagInd, 3])
    
    # speeds and CTMC generator decompositions
    for(i in 1:3) {
      log_lambda[tagInd, i] ~ dlogGamma(shape = lambda_priors[i, 1], 
                                        rate = lambda_priors[i, 2])
      
      # transformed parameters, for sampling
      lambda[tagInd, i] <- exp(log_lambda[tagInd, i])
      
      # generators
      expm_decomp[tagInd, i, 1:(2*N_bins+3), 1:N_bins] <- 
        buildAndDecomposeGenerator(
          pi = pi[tagInd, i], lambda = lambda[tagInd, i], M = N_bins, stage = i,
          widths = widths[1:N_bins], delta = delta, t = tstep
        )
      
    }
    
  }
  
  # dive endpoint random effects
  for(i in 1:N_endpoints) {
    endpoints[i] ~ dunif(endpoint_priors[i, 1], endpoint_priors[i, 2])
  }
  
  #
  # likelihood, and additional priors
  #

  for(diveId in 1:N_dives) {

    # stage duration priors
    log_xi[diveId, 1:2] ~ dmnorm(
      mean = xi_prior_means[dive_relations[diveId, 3], 1:2], 
      cov = xi_prior_covs[dive_relations[diveId, 3], 1:2, 1:2]
    )
    
    # back-transform stage durations
    xi[diveId, 1] <- exp(log_xi[diveId, 1])
    xi[diveId, 2] <- exp(log_xi[diveId, 2])
    
    # dive start and end times
    T[diveId, 1] <- endpoints[dive_relations[diveId, 1]]
    T[diveId, 4] <- endpoints[dive_relations[diveId, 2]]
    
    # stage transition times
    T[diveId, 2] <- T[diveId, 1] + xi[diveId, 1]
    T[diveId, 3] <- T[diveId, 2] + xi[diveId, 2]

    # likelihood
    depths[dive_relations[diveId, 4]:dive_relations[diveId, 5]] ~ ddive(
      times = times[dive_relations[diveId, 4]:dive_relations[diveId, 5]],
      N = dive_relations[diveId, 5] - dive_relations[diveId, 4] + 1,
      expm = expm_decomp[dive_relations[diveId, 3], 1:3, 1:N_bins, 1:N_bins],
      evecs = expm_decomp[dive_relations[diveId, 3], 1:3, 
                          (N_bins + 1):(2 * N_bins), 1:N_bins],
      evals = expm_decomp[dive_relations[diveId, 3], 1:3, 2 * N_bins + 1, 
                          1:N_bins],
      d = expm_decomp[dive_relations[diveId, 3], 1:3, 2 * N_bins + 2, 1:N_bins],
      dInv = expm_decomp[dive_relations[diveId, 3], 1:3, 2* N_bins + 3, 
                         1:N_bins],
      tstep = tstep, M = N_bins, T = T[diveId, 1:4], include = 1
    )
  }
  
})
