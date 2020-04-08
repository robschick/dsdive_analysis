# demonstrating use of Rdsm package for distributed likelihood computation
#
# Note that we use shared memory to aggregate the likelihood components, but the 
# shared memory approach is only marginally faster than simply distributing the 
# likelihood calculation across the cluster using standard R parallel methods.
#
# Shared memory systems, then, may be more useful for problems where the 
# parallel components are tightly-coupled within the parallel call.
#  
# we may still want to use Rdsm to enable parallelization in our workflow 
# because forking is the only other shared-memory tool available to us, and 
# there may be instances where the increased memory overhead of forking is 
# undesirable.  the Rdsm tool additionally allows us to potentially use lock 
# and barrier operations to implement the entire gibbs sampler update with a 
# single call to clusterEvalQ.  
#
# we may also prefer to use Rdsm because the idea behind forking on R is that 
# the nodes have a read-only copy of the main memory, whereas we need to modify 
# memory when we update the random effects in parallel.  However, forking is 
# sufficient for evaluating likelihoods in parallel.

library(parallel)
library(Rdsm)
library(synchronicity)


#
# simulate data
#

n = 1e7
x = rnorm(n = n)


#
# initialize cluster
#

# create Rdsm-enabled cluster 
cl = makeCluster(spec = 3, type = 'FORK')
mgrinit(cls = cl, boost = TRUE, barrback = FALSE)

# initialize shared-memory for batched likelihood contributions
mgrmakevar(cls = cl, varname = 'll.contrib', nr = 1, nc = length(cl), 
           vartype = 'double', savedesc = FALSE)

# distribute data across cluster (not strictly necessary with FORK)
clusterExport(cl, c('x', 'n'))

# pre-compute the indices that each node will compute the likelihood for
clusterEvalQ(cl, { myidxs <- getidxs(n) })

#
# evaluate likelihood in parallel
#

gc()

clusterEvalQ(cl, gc())

tick.parallel = proc.time()[3]

clusterEvalQ(cl = cl, expr = {
    # batch-compute likelihood contributions
    ll.contrib[myinfo$id] = sum(
      dnorm(x = x[myidxs], mean = 0, sd = 1, log = TRUE)
    )
    # null return
    invisible(0)
})
  
# aggregate likelihood contributions
ll.parallel = sum(ll.contrib[])

tock.parallel = proc.time()[3]


#
# evaluate likelihood serially, for comparison
#

gc()

tick.serial = proc.time()[3]

ll.serial = sum(dnorm(x, log = TRUE))

tock.serial = proc.time()[3]


#
# compare results
#

# relative error
abs(ll.serial-ll.parallel)/abs(ll.serial) * 100

# the results show near-linear speedup
(tock.parallel - tick.parallel) / (tock.serial - tick.serial)

stopCluster(cl)
