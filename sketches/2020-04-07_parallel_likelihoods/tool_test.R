# demonstrating use of Rdsm package for distributed likelihood computation

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

# distribute data across cluster (not strictly with FORK)
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
