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
    
    # only allow (downward) transitions from shallowest bin during stages 1 & 2
    if(stage == 1 | stage == 2) {
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
      entries$diag[2] <- -rate
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
    
    evals <- numeric(length = N, init = FALSE)
    evals[1:N] <- Aentries$diag[1:N]

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


# # compile nimble code, so we can use it
# cdecomposeGenerator = compileNimble(decomposeGenerator)
  
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
  run = function(expm = double(2), evecs = double(2), evals = double(1), 
                 d = double(1), dInv = double(1), tstep = double(0), 
                 N = integer(0), v = double(1), t = double(0), 
                 preMultiply = logical(0)) {
    
    returnType(double(1))
    
    x <- numeric(0, length = N)
    
    if(t == tstep) {
      
      if(preMultiply == TRUE) {
        x[1:N] <- t(v[1:N]) %*% expm[1:N,1:N]
      } else {
        x[1:N] <- expm[1:N,1:N] %*% v[1:N]
      }
    } else {
      expmAtvCpp(evecs = evecs, evals = evals, M = N, v = v, d = d, dInv = dInv, 
                 t = t, x = x, preMultiply = preMultiply)
    }
    
    return(x)
  }
)

# # compile nimble code, so we can use it
# cexpmAtv = compileNimble(expmAtv)

buildAndDecomposeGenerator = nimbleFunction(
  run = function(pi = double(0), lambda = double(0), M = integer(0),
                 stage = integer(0), widths = double(1), delta = double(0), 
                 t = double(0)) {
    
    returnType(double(2))
    
    res <- matrix(nrow = 2 * M + 3, ncol = M, init = FALSE)
    
    # build generator entries
    entries <- buildInfinitesimalGeneratorEntries(
      pi = pi, lambda = lambda, M = M, stage = stage, widths = widths
    )
    
    # decompose generator
    decomposition <- decomposeGenerator(Aentries = entries, delta = delta, 
                                        t = t, N = M)
    
    # assign output
    res[1:M,1:M] <- decomposition$expm[1:M,1:M]
    res[(M+1):(2*M),1:M] <- decomposition$evecs[1:M,1:M]
    res[2*M+1, 1:M] <- decomposition$evals[1:M]
    res[2*M+2,1:M] <- decomposition$d[1:M]
    res[2*M+3,1:M] <- decomposition$dInv[1:M]
    
    return(res)
  }
)

# cbuildAndDecomposeGenerator = compileNimble(buildAndDecomposeGenerator)
# 
# cbuildAndDecomposeGenerator(pi = .95, lambda = 1.5, M = 16, stage = 1, 
#                            widths = 2*depth.bins$halfwidth, delta = 1e-15, 
#                            t = 300)

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
# prob = cexpmAtv(expm = Adecomp$expm, evecs = Adecomp$evecs, 
#                 evals = Adecomp$evals, d = Adecomp$d, dInv = Adecomp$dInv, 
#                 tstep = Adecomp$t, N = Adecomp$N, v = x0, t = t0, 
#                 preMultiply = TRUE)
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

# likelihood for an entire dive
ddive = nimbleFunction(
  run = function(x = double(1), times = double(1), expm = double(3), 
                 N = integer(0),
                 evecs = double(3), evals = double(2), d = double(2), 
                 dInv = double(2), tstep = double(0), 
                 M = integer(0), T = double(1), log = logical(0, default = 0)) {
    # Parameters:
    #   x - sequence of observed depth bins 
    #   times - observation times
    #   N - number of observations
    #   M - number of depth bins
    #   T - vector of dive-specific random effects;  c(T0, T1, T2, T3) s.t. 
    #       T0 is dive start time, T1 and T2 are respecitively end of stage 1/2,
    #       and T3 is dive end time.  note that the indexing implies T[1] is 
    #       start of dive, T[2] is start of stage 2, T[3] is start of stage 3, 
    #       and T[4] is start of surface period (i.e., end of dive).
    #   log - If TRUE, then log-likelihood is returned
    
    returnType(double(0))
    
    # density is not defined for poorly-ordered dive-specific random effects
    if(T[1] >= T[2]) { return(-Inf) }
    if(T[2] >= T[3]) { return(-Inf) }
    if(T[3] >= T[4]) { return(-Inf) }
    
    # initialize log-likelihood
    ll <- 0
    
    # initialize Markov transitions from beginning of dive
    x_prev <- 1
    time_prev <- T[1]
    stage_prev <- 1
    
    # TODO: finish stage 3, account for tx. from last obs. to the surface!

    # loop over observations
    more_observations <- TRUE
    for(i in 1:N) {
      # ad-hoc flag to "break" loop
      if(more_observations == TRUE) {
        # only process in-sync observations (e.g., during dive window)
        if(times[i] > time_prev) {
          
          # extract current transition: special handling for last observation
          if(times[i] >= T[4]) {
            x_loc <- 1
            time_loc <- T[4]
            more_observations <- FALSE
          } else {
            x_loc <- x[i]
            time_loc <- times[i]
          }
          
          # determine dive stage
          if(time_loc >= T[3]) {
            stage_loc <- 3
          } else if(time_loc >= T[2]) {
            stage_loc <- 2
          } else {
            stage_loc <- 1
          }
          
          # initialize transition distribution
          u <- numeric(0, length = M)
          u[x_prev] <- 1
          
          # diffuse transition mass across stages
          time_stage <- time_prev
          for(s in stage_prev:stage_loc) {
            # time spent in stage
            d_t <- min(T[s+1], time_loc) - time_stage
            # update time-in-stage counter
            time_stage <- T[s+1]
            # diffuse mass across stage s
            u[1:M] <- expmAtv(expm = expm[s,1:M,1:M], evecs = evecs[s,1:M,1:M], 
                              evals = evals[s,1:M], d = d[s,1:M], 
                              dInv = dInv[s,1:M], tstep = tstep, N = M, 
                              v = u[1:M], t = d_t, preMultiply = TRUE)[1:M]
          }
          
          # aggregate likelihood
          ll <- ll + log(u[x_loc])
          
          # update last state, for computing Markov probabilities
          x_prev <- x_loc
          time_prev <- time_loc
          stage_prev <- stage_loc
        }
      }
      
    }
    
    if(log) { return(ll) } else { return(exp(ll)) }
  }
)

# transition kernel for a single dive transition
ddiveTx = nimbleFunction(
  run = function(x = integer(0), x_prev = integer(0), expm = double(3), 
                 evecs = double(3), evals = double(2), d = double(2), 
                 dInv = double(2), tstep = double(0), 
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
          u[1:M] <- expmAtv(expm = expm[s,1:M,1:M], evecs = evecs[s,1:M,1:M], 
                            evals = evals[s,1:M], d = d[s,1:M], 
                            dInv = dInv[s,1:M], tstep = tstep, N = M, 
                            v = u[1:M], t = d_t, preMultiply = TRUE)[1:M]
        }
      }
      
      # extract transition probability
      prob <- u[x_local]
      
    }
    
    if(log) { return(log(prob)) } else { return(prob) }
  }
)


# #
# # Test likelihood
# #
# 
# pi = c(.95, .5, .05)
# lambda = c(1.5, .3, 1)
# 
# cdecomposeGenerator = compileNimble(decomposeGenerator)
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
# cddive = compileNimble(ddive)
# 
# expm = array(dim = c(3,16,16))
# evecs = array(dim = c(3,16,16))
# evals = matrix(nrow = 3, ncol = 16)
# d = matrix(nrow = 3, ncol = 16)
# dInv = matrix(nrow = 3, ncol = 16)
# 
# for(s in 1:3) {
#   expm[s,,] = Adecomps[[s]]$expm
#   evecs[s,,] = Adecomps[[s]]$evecs
#   evals[s,] = Adecomps[[s]]$evals
#   d[s,] = Adecomps[[s]]$d
#   dInv[s,] = Adecomps[[s]]$dInv
# }
# 
# 
# cddiveTx(x = 5, x_prev = 1, expm = expm, evecs = evecs, evals = evals, d = d,
#          dInv = dInv, tstep = 300, M = 16, t_prev = 0, t = 300,
#          T = c(0, 600, 800, Inf), log = TRUE)
# 
# cddive(x = dives.obs[[1]]$dive$depths, times = dives.obs[[1]]$dive$times, 
#        expm = expm, N = length(dives.obs[[1]]$dive$times), evecs = evecs, 
#        evals = evals, d = d, dInv = dInv, tstep = 300, M = 16, 
#        T = c(-10, 1000, 2000, 3000), log = TRUE)