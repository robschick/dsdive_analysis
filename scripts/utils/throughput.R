# throughput analysis; check memory and time limits for imputations.
# goal is to determine computational needs to post-process dives

# formatting of bytes objects
library(pryr)


#
# costs
#

# average size of a single imputed dive (bytes object)
imputed.size = 1700
class(imputed.size) = 'bytes'

# average time to impute a single dive (seconds)
imputed.time = 6.5/50


#
# needs
#

# number of dives to impute
n.dives = 132

# number of posterior samples to work over
mcit = 10e3

# total size of dive imputations
total.size = imputed.size * n.dives * mcit

# total time required for computation
total.time = imputed.time * n.dives * mcit

  
#
# constraints
#

# maximum core memory size (bytes)
max.coremem = 1024^3
class(max.coremem) = 'bytes'

# desired wall time to finish imputation
desired.time = 3600


#
# computing requirements
#

# minimum number of cores required based on memory
min.cores.mem = ceiling(total.size / max.coremem)

# minimum number of cores required based on desired time
min.cores.time = ceiling(total.time / desired.time)

# minimum number of cores required to satisfy all constraints
min.cores = max(min.cores.mem, min.cores.time)

