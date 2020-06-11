library(dsdive)

#
# dsdive likelihood, etc.
#

data('dive.sim')
attach(dive.sim)
attach(dive.sim$params)

# extract time between observations
tstep = diff(sim.obs$times[1:2])

# get true stage transition times
t.stages = sim$times[c(FALSE, diff(sim$stages)==1)]

# build probability matrices for observations
obstx.mat = lapply(1:3, function(s) {
  dsdive.obstx.matrix(depth.bins = depth.bins, beta = dive.sim$params$beta, 
                      lambda = lambda, 
                      s0 = s, tstep = tstep, include.raw = TRUE, delta = 1e-10)
})

# compute likelihood of observations, given model parameters
ld = dsdive.obsld(dsobs.list = list(sim.obs),
                  t.stages.list = list(t.stages), 
                  P.raw = obstx.mat, s0 = 1, sf = 3)


#
# nimble likelihood, etc.
#

pi = c(dive.sim$params$beta[1], .5, dive.sim$params$beta[2])

cbuildAndDecomposeGenerator = compileNimble(buildAndDecomposeGenerator)

obstx.mat_nimble = lapply(1:3, function(s) {
  cbuildAndDecomposeGenerator(pi = pi[s], lambda = lambda[s], M = 16, stage = s, 
                              widths = 2 * depth.bins$halfwidth, delta = 1e-10,
                              t = tstep)
})

A_nimble = lapply(1:3, function(s) {
  buildInfinitesimalGenerator(pi = pi[s], lambda = lambda[s], M = 16, stage = s, 
                              widths = 2 * depth.bins$halfwidth)
})


expm_decomp = array(dim = c(3,2*16+3,16))
M = 16
widths = 2 * depth.bins$halfwidth
delta = 1e-10
for(s in 1:3) {
  expm_decomp[s,1:(2*M+3),1:M] <- cbuildAndDecomposeGenerator(
    pi = pi[s], lambda = lambda[s], M = M, stage = s, widths = widths[1:M], 
    delta = delta, t = tstep)
}


for(s in 1:3) {
  expm_decomp[s,1:M,1:M] = obstx.mat[[s]]$obstx.mat
  expm_decomp[s,(M+1):(2*M),1:M] = obstx.mat[[s]]$vectors
  expm_decomp[s,2*M+1,1:M] = obstx.mat[[s]]$values
  expm_decomp[s,2*M+2,1:M] = obstx.mat[[s]]$d
  expm_decomp[s,2*M+3,1:M] = obstx.mat[[s]]$dInv
}

cddive = compileNimble(ddive)

cexpmAtv = compileNimble(expmAtv)

# compute likelihood of observations, given model parameters
ld_nimble = ddive(x = sim.obs$depths, times = sim.obs$times, 
                   expm = expm_decomp[1:3,1:M,1:M],
                   evecs = expm_decomp[1:3,(M+1):(2*M),1:M],
                   evals = expm_decomp[1:3,2*M+1,1:M],
                   d = expm_decomp[1:3,2*M+2,1:M],
                   dInv = expm_decomp[1:3,2*M+3,1:M],
                   N = length(sim.obs$depths), tstep = 300, M = 16, 
                   T = c(0,t.stages,sim.obs$times[length(sim.obs$times)]),
                   log = TRUE)

ld = dsdive.obsld(dsobs.list = list(sim.obs),
                  t.stages.list = list(t.stages), 
                  P.raw = obstx.mat, s0 = 1, sf = 3)



#
# comparisons
#

sapply(1:3, function(s) {
  c(stage = s,
    tx_mat = max(abs(obstx.mat[[s]]$obstx.mat - obstx.mat_nimble[[s]][1:16,])),
    vecs = max(abs(obstx.mat[[s]]$vectors - obstx.mat_nimble[[s]][16 + 1:16,])),
    vals = max(abs(obstx.mat[[s]]$values - obstx.mat_nimble[[s]][33,])),
    d = max(abs(obstx.mat[[s]]$d - obstx.mat_nimble[[s]][34,]))
  )
})

sapply(1:3, function(s) {
  c(stage = s,
    tx_mat = max(abs(obstx.mat[[s]]$obstx.mat - expm_decomp[s,1:16,])),
    vecs = max(abs(obstx.mat[[s]]$vectors - expm_decomp[s,16 + 1:16,])),
    vals = max(abs(obstx.mat[[s]]$values - expm_decomp[s,33,])),
    d = max(abs(obstx.mat[[s]]$d - expm_decomp[s,34,]))
  )
})


max(abs(A_nimble[[1]] - obstx.mat[[1]]$A))
max(abs(A_nimble[[2]] - obstx.mat[[2]]$A))
max(abs(A_nimble[[3]] - obstx.mat[[3]]$A))

which.max(abs(obstx.mat[[3]]$vectors - obstx.mat_nimble[[3]][16 + 1:16,]))

obstx.mat[[3]]$values

detach(dive.sim$params)
detach(dive.sim)