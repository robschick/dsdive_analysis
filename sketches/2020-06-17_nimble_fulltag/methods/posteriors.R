library(coda)

load('sketches/2020-06-04_nimble_implementation/methods/samples.RData')

dim(samples)

burn = 1:1e3


post.samples = samples[-burn, setdiff(colnames(samples), 
                                      c('logit_pi[2]', 'pi[2]'))]

min(effectiveSize(mcmc(post.samples)))

min(effectiveSize(mcmc(samples[-burn, c('pi[1]', 'pi[3]')])))
min(effectiveSize(mcmc(samples[-burn, c('lambda[1]', 'lambda[2]', 'lambda[3]')])))

nrow(samples)/900 * 1e3


plot(mcmc(samples[-burn, c('pi[1]', 'pi[3]')]))
plot(mcmc(samples[-burn, c('lambda[1]', 'lambda[2]', 'lambda[3]')]))

1-rejectionRate(mcmc(samples[-burn, c('pi[1]', 'lambda[1]')]))
1-rejectionRate(mcmc(samples[-burn, c('pi[3]', 'lambda[3]')]))


plot(density(as.numeric(samples[-burn, paste('T[', 1:47, ', 1]', sep ='')])),
     xlab = 'Start offet (s)', main = 'Posterior density across all dives')

# NOTE: End offsets require comparing the posterior distribution of T[,4] to 
# the last observation time

plot(density(as.numeric(samples[-burn, paste('xi[', 1:47, ', 1]', sep ='')])/60),
     xlab = 'Stage 1 duration (min)', main = 'Posterior density across all dives')

plot(density(as.numeric(samples[-burn, paste('xi[', 1:47, ', 2]', sep ='')])/60),
     xlab = 'Stage 2 duration (min)', main = 'Posterior density across all dives')


plot(density((as.numeric(samples[-burn, paste('T[', 1:47, ', 4]', sep ='')]) - 
             as.numeric(samples[-burn, paste('T[', 1:47, ', 1]', sep ='')]))/60),
     xlab = 'Dive duration (min)', main = 'Posterior density across all dives')
