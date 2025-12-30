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
source("R/nimble_models/nimbleCode_gamma_2.q")

# CONSTANTS USED FOR INITIAL INITIAL VALUES AND FOR PRIORS
colsumy <- apply(nimbleData$data$y, 3, sum, na.rm=TRUE)
p1hat <- colsumy[3]/(colsumy[2]+colsumy[3])
p2hat <- colsumy[3]/(colsumy[1]+colsumy[3])
Nhat <- apply(nimbleData$data$Y, 1, sum, na.rm=TRUE)/(1-(1-p1hat)*(1-p2hat)) # For all sites combined
lambdahat <- nimbleData$const$lambdahat

N <- round(nimbleData$data$Y/(1-(1-p1hat)*(1-p2hat)), 0)
N[is.na(N)] <- 0
nrowY <- nrow(nimbleData$data$Y)
ncolY <- ncol(nimbleData$data$Y)

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

# SETTING UP THE MCMC
DoubleObsMultisiteModel <- nimbleModel(
  nimbleCode_DOMM_gamma_2,
  constants = list(
                   N_surv = length(nimbleData$const$N_sites),
                   N_sites = nimbleData$const$N_sites,
                   prior_mean_eta = prior_parameters_for_p$mean_eta,
                   prior_sd_eta = prior_parameters_for_p$sd_eta,
                   sam = nimbleData$const$sam,
                   N_sam = length(unique(nimbleData$const$sam)),
                   S = sapply(1:8, function(i) range(which(nimbleData$const$sam==i))),
                   area = nimbleData$const$area,
                   x = standardize(nimbleData$const$mean_field_dist)
  ),
  data = nimbleData$data,
  inits = Inits()
)

t1 <- Sys.time()
CDoubleObsMultisiteModel <- compileNimble(DoubleObsMultisiteModel) # Needs to be compiled for the last step
DoubleObsMultisiteConf <- configureMCMC(DoubleObsMultisiteModel, 
                                        monitors = c("Disc_New_Y", "Disc_Y", "Disc_New_y", "Disc_y", "eta1", "eta2", "mu0", "sigma", "sigma_p", "beta", "mean_ds"),
                                        enableWAIC = TRUE)

DoubleObsMultisiteMCMC <- buildMCMC(DoubleObsMultisiteConf)
CDoubleObsMultisiteMCMC <- compileNimble(DoubleObsMultisiteMCMC)

t2 <- Sys.time()
cat("Compilation time:")
t2-t1

t3 <- Sys.time()
init.values <- list(Inits(), Inits(), Inits())
settings <- list(
  niter = 200000,
  nburnin = 20000,
  nchain = 3,
  thin = 6
)
out <- runMCMC(
  CDoubleObsMultisiteMCMC,
  niter = settings$niter,
  nburnin = settings$nburnin,
  nchain = settings$nchain,
  thin = settings$thin,
  inits = init.values,
  samplesAsCodaMCMC = TRUE,
  WAIC = TRUE)
t4 <- Sys.time()
cat("Run time:")
t4-t3

# Saving workspace
save(settings, out, file = "posterior_samples/gamma_2.RData")

#plot(out$samples) # 1 = black, 2 = red, 3 = green
# summary(out$samples)
# gelman.diag(out$samples)

# out$WAIC
# 
# # Posterior predictive checks
#
samp <- as.matrix(out$samples)
plot(samp[,"Disc_New_Y"] ~ samp[,"Disc_Y"], xlim = c(0, 5000), ylim = c(0, 5000), pch = ".")
abline(0, 1, col="red")
mean(samp[, "Disc_New_Y"] > samp[, "Disc_Y"])

lim <- 6000
New_Y_trunc <- samp[, "Disc_New_Y"][samp[, "Disc_New_Y"] < lim]
Y_trunc <- samp[, "Disc_Y"][samp[, "Disc_Y"] < lim]

plot(density(Y_trunc), main = "", col = "red", xlab = "Sum of squared\n Pearson residuals")
lines(density(New_Y_trunc))
legend("topright",
       lty = 1,
       col = c("black", "red"),
       legend = c("Simulated data", "Observed data"),
       cex = 0.8
)

# Wrt y
plot(samp[,"Disc_New_y"] ~ samp[,"Disc_y"])
abline(0, 1, col="red")
mean(samp[, "Disc_New_y"] > samp[, "Disc_y"])
# 
