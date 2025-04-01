# R-packages used: 
library(nimble)
library(coda)

# Making standardize function (NB! Do not use scale() as this function scales each column of a matrix)
standardize <- function(x) (x-mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE)

# Loading prior-parameters for p
load(file = "data/prior_parameters_for_p.rda")

# Loading data
load(file = "data/nimbleData.rda")

# READING THE MODEL CODE
source("R/nimble_models/nimbleCode_gamma_2_noPPC.q")

load("posterior_samples/gamma_2.RData")
postsamp <- out$samples
#best_postsamp <- as.mcmc.list(best_postsamp)
samp <- as.matrix(postsamp)[,-c(1:4, 9:16)]
p <- apply(samp, 2, median)
truepar <- list(
  beta = p[1:2],
  eta1 = p[3],
  eta2 = p[4],
  mu0 = p[5:12],
  sigma = p[13],
  sigma_p = p[14]
)

# REDUCING sigma
alpha <- 1/20
truepar$mu0 <- truepar$mu0 + truepar$sigma^2*(1 - alpha^2)/2 # NB! using old sigma, so must come first
truepar$sigma <- truepar$sigma*alpha #NB! Run only once


### remaking plot from Supp Inf with modified sigma #####
qs <- function(mean, sigma, p){
  a <- 1/(exp(sigma^2) - 1)
  b <- a/mean
  qgamma(p, shape = a, rate = b)
}

mean <- seq(0.00001, 2, 0.1)
sigma <- truepar$sigma
a <- 1/(exp(sigma^2) - 1)
b <- a/mean
marginal_mean_density <- 0.2621269

qQ.50 <- qs(mean, sigma, 0.5)
qQ.90 <- qs(mean, sigma, 0.9)
qQ.975 <- qs(mean, sigma, 0.975)

par(mfrow = c(1,2))
plot(mean, qQ.975, type = "n", col = "blue", ylim = c(0, max(qQ.975)),
     xlab = "Mean deer density (E[X])",
     ylab = "Deer density (x)")

lines(mean, qQ.975, type = "l", col = "blue")
lines(mean, qQ.90, type = "l", col = "darkgreen")
lines(mean, qQ.50, type = "l", col = "red")
abline(v = marginal_mean_density, lty = 2)
legend("topleft",
       lty = 1,
       col = c("blue", "darkgreen", "red"),
       legend = c("97.5% quantile", "90.0% quantile", "50.0% quantile"),
       cex = 0.8
)

b <- a/marginal_mean_density
Median <- qs(marginal_mean_density, sigma, p = 0.5)
Q90 <- qs(marginal_mean_density, sigma, p = 0.9)
Q975 <- qs(marginal_mean_density, sigma, p = 0.975)

x <- seq(0.0001, 0.45, length.out = 100)
plot(dgamma(x, a, b) ~ x, type = "n", xlab = "Deer density (x)", ylab = "f(x)", xaxs = "i", yaxs = "i", xlim = c(0,max(x)), ylim = c(0, 1.1*max(dgamma(x, a, b))))
box(col = "grey")
lines(dgamma(x, a, b) ~ x, lwd = 2)
abline(v = Median, col = "red")
abline(v = Q90, col = "darkgreen")
abline(v = Q975, col = "blue")
abline(v = a/b, lty = 2)
# legend("topright",
#        lty = 1,
#        col = c("blue", "darkgreen", "red"),
#        legend = c("97.5% quantile", "90.0% quantile", "50.0% quantile"),
#        cex = 0.8
# )

########

# Setting up constants to be used for both simulation and model fitting
N_surv <- length(truepar$mu0)
n_sim_sites <- 500 # Number of sites in each survey
area_val <- as.vector(nimbleData$const$area)
area_val <- area_val[!is.na(area_val)]
x_val <- as.vector(standardize(nimbleData$const$mean_field_dist))
x_val <- x_val[!is.na(x_val)]
constants <- list(
  N_surv = N_surv,
  N_sites = rep(n_sim_sites, N_surv),
  prior_mean_eta = prior_parameters_for_p$mean_eta,
  prior_sd_eta = prior_parameters_for_p$sd_eta,
  sam = 1:N_surv,
  N_sam = N_surv,
  area = array(sample(area_val, N_surv*n_sim_sites, replace = TRUE), dim = c(N_surv, n_sim_sites)),
  x = array(sample(x_val, N_surv*n_sim_sites, replace = TRUE), dim = c(N_surv, n_sim_sites))
)

# Building model
DoubleObsMultisiteModel <- nimbleModel(
  code = nimbleCode_DOMM_gamma_2,
  constants = constants,
  inits = truepar
)

# Compile model
CDoubleObsMultisiteModel <- compileNimble(DoubleObsMultisiteModel)

# Simulation
nodesToSim <- CDoubleObsMultisiteModel$getDependencies(names(p), self = F, downstream = T)
CDoubleObsMultisiteModel$simulate(nodesToSim)
sim_y <- CDoubleObsMultisiteModel$y
sim_Y <- CDoubleObsMultisiteModel$Y

# Fitting model with mcmc
CDoubleObsMultisiteModel$setData(list(y = sim_y, Y = sim_Y))
simMCMC <- buildMCMC(CDoubleObsMultisiteModel)
CsimMCMC <- compileNimble(simMCMC, project = DoubleObsMultisiteModel)
#samples <- runMCMC(CsimMCMC, niter = 10000)

# Initial values (using the original procedure)
Y_per_sam <- apply(sim_Y, 1, sum, na.rm=TRUE) |> tapply(INDEX = constants$sam, sum)

# Simple estimate of detection probability
colsumy <- apply(sim_y, 3, sum, na.rm=TRUE)
p1hat <- colsumy[3]/(colsumy[2]+colsumy[3])
p2hat <- colsumy[3]/(colsumy[1]+colsumy[3])
phat_simple <- 1 - (1-p1hat)*(1-p2hat)

# Area per sam
area_per_sam <- apply(constants$area, 1, sum, na.rm=TRUE) |> tapply(INDEX = constants$sam, sum)

lambdahat <- nimbleData$const$lambdahat
# # Simple estimates of density
# lambdahat <- (Y_per_sam/phat_simple)/area_per_sam
# 
# # For the sam's with zero detentions (see prior_zero_sites.R)
# lambdahat[3] <- 0.0066 # Raa (April)
# lambdahat[4] <- 0.0046 # Raa (March)
# lambdahat[7] <- 0.0047 # Sprakehaug (April)

N <- round(sim_Y/(1-(1-p1hat)*(1-p2hat)), 0)
N[is.na(N)] <- 0
nrowY <- nrow(sim_Y)
ncolY <- ncol(sim_Y)

Inits <- function(){
  sigma_p <- runif(1, 0.1, 0.4)
  p1 <- exp(log(p1hat)*runif(1, 0.9, 1.1))
  p2 <- exp(log(p2hat)*runif(1, 0.9, 1.1))
  eta1 <- log(p1/(1-p1))
  eta2 <- log(p2/(1-p2))
  list(
    mu0 = log(lambdahat*runif(length(lambdahat), 0.9, 1.1)),
    eta1 = eta1,
    eta2 = eta2,
    logit_p1 = matrix(rnorm(nrowY*ncolY, eta1, sigma_p/2) , nrow = nrowY, ncol = ncolY),
    logit_p2 = matrix(rnorm(nrowY*ncolY, eta2, sigma_p/2) , nrow = nrowY, ncol = ncolY),
    N = N,
    lambda = N + 0.01,
    beta = runif(2, -0.5, 0.5),
    sigma = runif(1, 0.5, 1),
    sigma_p = sigma_p
  )
}


init.values <- list(Inits(), Inits(), Inits())
settings <- list(
  niter = 200000,
  nburnin = 20000,
  nchain = 3,
  thin = 6
)

out <- runMCMC(
  CsimMCMC,
  niter = settings$niter,
  nburnin = settings$nburnin,
  nchain = settings$nchain,
  thin = settings$thin,
  inits = init.values,
  samplesAsCodaMCMC = TRUE
  )

# Saving workspace
save(settings, out, file = "posterior_samples/sim/sim_gamma_2_500_sites_low_sigma.RData")
