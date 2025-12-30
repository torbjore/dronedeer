library(nimble)
library(coda)

# DEFINING THE MODEL
code <- nimbleCode({
  Y ~ dpois(lambda)
  lambda <- lambda0*area
  lambda0 <- exp(mu0)
  mu0 ~ dnorm(mu0mean, sd = 5)
  mu0mean ~ dnorm(0, sd=5)
})

# Raa April
data <- list(Y = 0)
inits <- function() list(mu0 = rnorm(1))
mcmc.output <- nimbleMCMC(code, 
                          constants = list(area = 5.356925),
                          data = data, inits = inits,
                          monitors = c("lambda0"), thin = 1,
                          niter = 50000, nburnin = 5000, nchains = 3,
                          samplesAsCodaMCMC = TRUE)

plot(mcmc.output)
gelman.diag(mcmc.output)
summary(mcmc.output)
# lambda0 = 6.621e-03 = 0.0066
# with sd = 2.35: lambda0 = 5.634e-03
# with sd = 5: lambda0 = 2.096e-03

# Raa March
data <- list(Y = 0)
inits <- function() list(mu0 = rnorm(1))
mcmc.output <- nimbleMCMC(code, 
                          constants = list(area = 10.580784),
                          data = data, inits = inits,
                          monitors = c("lambda0"), thin = 1,
                          niter = 50000, nburnin = 5000, nchains = 3,
                          samplesAsCodaMCMC = TRUE)

plot(mcmc.output)
gelman.diag(mcmc.output)
summary(mcmc.output)
# lambda0 = 4.607e-03 = 0.0046
# with sd = 2.35: 3.515e-03
# with sd = 5: lambda0 = 1.174e-03

# Sprakehaug (April)
mcmc.output <- nimbleMCMC(code, 
                          constants = list(area = 10.198860),
                          data = data, inits = inits,
                          monitors = c("lambda0"), thin = 1,
                          niter = 50000, nburnin = 5000, nchains = 3,
                          samplesAsCodaMCMC = TRUE)
plot(mcmc.output)
gelman.diag(mcmc.output)
summary(mcmc.output)
# lambda0 = 4.734e-03 = 0.0047
# with sd = 2.35: 3.756e-03
# with sd = 5: lambda0 = 1.268e-03
