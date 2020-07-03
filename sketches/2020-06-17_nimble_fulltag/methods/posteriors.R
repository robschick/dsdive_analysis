library(coda)

load('sketches/2020-06-17_nimble_fulltag/methods/samples.RData')

dim(samples)

burn = 1:1e3

varnames = colnames(samples)


post.samples = samples[-burn, setdiff(varnames, c('logit_pi[2]', 'pi[2]'))]

min(effectiveSize(mcmc(post.samples)))

min(effectiveSize(mcmc(samples[-burn, c('pi[1]', 'pi[3]', 'pi[4]', 'pi[5]')])))
min(effectiveSize(mcmc(samples[-burn, c('lambda[1]', 'lambda[2]', 'lambda[3]',
                                        'lambda[4]', 'lambda[5]')])))


plot(mcmc(samples[-burn, c('pi[1]', 'pi[3]', 'pi[4]', 'pi[5]')]))
plot(mcmc(samples[-burn, c('lambda[1]', 'lambda[2]', 'lambda[3]', 
                           'lambda[4]', 'lambda[5]')]))

1-rejectionRate(mcmc(samples[-burn, c('pi[1]', 'lambda[1]')]))
1-rejectionRate(mcmc(samples[-burn, c('pi[3]', 'lambda[3]')]))


#
# slow mixing of E terms, and some don't even move; need better sampler inits, 
#   or maybe use slice samplers
#


samples.e = samples[-burn, grep(pattern = 'E\\[[0-9]+\\]', x = varnames)]

min(effectiveSize(mcmc(samples.e[, -which(apply(samples.e, 2, var) == 0)])))
