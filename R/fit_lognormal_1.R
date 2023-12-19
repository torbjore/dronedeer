# R-packages used: 
library(nimble)
library(coda)

# Making standardize function (NB! Do not use scale() as this function scales each column of a matrix)
standardize <- function(x) (x-mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE)

# Loading prior-parameters for p
load(file = "data/prior_parameters_for_p.rda")

# Loading data
load(file = "data/LDDdata.rda")

# READING THE MODEL CODE
source("R/nimble_models/nimbleCode_DOMM_lognormal_1.q")

# CONSTANTS USED FOR INITIAL INITUAL VALUES AND FOR PRIORS
colsumy <- apply(LDDdata$data$y, 3, sum, na.rm=TRUE)
p1hat <- colsumy[3]/(colsumy[2]+colsumy[3])
p2hat <- colsumy[3]/(colsumy[1]+colsumy[3])
Nhat <- apply(LDDdata$data$Y, 1, sum, na.rm=TRUE)/(1-(1-p1hat)*(1-p2hat)) # For all sites combined
lambdahat_surv <- (Nhat/LDDdata$const$N_sites)/apply(LDDdata$const$area, 1, mean, na.rm=TRUE) + 0.01 # Adding a small value since we get -Inf from log(lambdahat=0)
lambdahat <- tapply(lambdahat_surv, list(LDDdata$const$sam), mean)
N <- round(LDDdata$data$Y/(1-(1-p1hat)*(1-p2hat)), 0)
N[is.na(N)] <- 0
nrowY <- nrow(LDDdata$data$Y)
ncolY <- ncol(LDDdata$data$Y)

# FUNCTION FOR INITIAL VALUES
Inits <- function(){
  sigma <- runif(1, 0.5, 1)
  #  sigma <- runif(1, 0.1, 0.3)
  sigma_p <- runif(1, 0.1, 0.4)
  p1 <- exp(log(p1hat)*runif(1, 0.9, 1.1))
  p2 <- exp(log(p2hat)*runif(1, 0.9, 1.1))
  mu_p1 <- log(p1/(1-p1))
  mu_p2 <- log(p2/(1-p2))
  list(
    mu0 = log(lambdahat*runif(length(lambdahat), 0.9, 1.1)),
    epsilon = matrix(rnorm(nrowY*ncolY, 0, sigma) , nrow = nrowY, ncol = ncolY),
    mu_p1 = mu_p1,
    mu_p2 = mu_p2,
    logit_p1 = matrix(rnorm(nrowY*ncolY, mu_p1, sigma_p/2) , nrow = nrowY, ncol = ncolY),
    logit_p2 = matrix(rnorm(nrowY*ncolY, mu_p2, sigma_p/2) , nrow = nrowY, ncol = ncolY),
    N = N,
    sigma = sigma,
    sigma_p = sigma_p,
    beta = runif(1, -0.5, 0.5),
    New_Y = N, # Warning message if not included
    New_y = LDDdata$data$y
  )
}

# SETTING UP THE MCMC
DoubleObsMultisiteModel <- nimbleModel(
  nimbleCode_DOMM_lognormal_1,
  constants = list(lamblow = 0.1*lambdahat,  # 0.1 to 10 times point estimate
                   lambupp = 10*lambdahat,
                   N_surv = length(LDDdata$const$N_sites),
                   N_sites = LDDdata$const$N_sites,
                   prior_mu_logit_p = prior_parameters_for_p$mu_logit_p,
                   prior_sigma_logit_p = prior_parameters_for_p$sigma_logit_p,
                   sam = LDDdata$const$sam,
                   N_sam = length(unique(LDDdata$const$sam)),
                   area = LDDdata$const$area,
                   x = standardize(LDDdata$const$mean_field_dist)
  ),
  data = LDDdata$data,
  inits = Inits()
)

t1 <- Sys.time()
CDoubleObsMultisiteModel <- compileNimble(DoubleObsMultisiteModel) # Needs to be compiled for the last step
DoubleObsMultisiteConf <- configureMCMC(DoubleObsMultisiteModel, 
#                                        monitors = c("p1", "p2", "mu0", "sigma", "beta"), enableWAIC = TRUE)
#                                        monitors = c("median_lambda", "mean_lambda" , "p1", "p2", "mu0", "sigma", "beta"), enableWAIC = TRUE)
#                                        monitors = c("Disc_New_Y", "Disc_Y", "Disc_New_y", "Disc_y", "median_lambda", "mean_lambda" , "mu_p1", "mu_p2", "mu0", "sigma", "sigma_p", "beta"), enableWAIC = TRUE)
                                        monitors = c("Disc_New_Y", "Disc_Y", "Disc_New_y", "Disc_y", "mu_p1", "mu_p2", "mu0", "sigma", "sigma_p", "beta"), enableWAIC = TRUE)


DoubleObsMultisiteMCMC <- buildMCMC(DoubleObsMultisiteConf)
CDoubleObsMultisiteMCMC <- compileNimble(DoubleObsMultisiteMCMC)
t2 <- Sys.time()

cat("Compilation time:")
t2-t1

t3 <- Sys.time()
posterior_lognormal_1 <- runMCMC(
  CDoubleObsMultisiteMCMC,
  niter=50000,
  nburnin=10000,
  nchain=3,
  thin= 4,
  inits = Inits,
  samplesAsCodaMCMC = TRUE,
  WAIC = TRUE)
t4 <- Sys.time()

cat("Run time:")
t4-t3

#load(file = "data/posterior_samples/posterior_lognormal_1.RData")

plot(posterior_lognormal_1$samples)
crosscorr.plot(posterior_lognormal_1$samples)
summary(posterior_lognormal_1$samples)
crosscorr(posterior_lognormal_1$samples)

posterior_lognormal_1$WAIC

gelman.diag(posterior_lognormal_1$samples)

#save(posterior_lognormal_1, file = "data/posterior_samples/posterior_lognormal_1.RData")

# Posterior predictive checks

# Wrt Y
samp <- as.matrix(posterior_lognormal_1$samples)
plot(samp[,"Disc_New_Y"] ~ samp[,"Disc_Y"])
abline(0, 1, col="red")
mean(samp[, "Disc_New_Y"] > samp[, "Disc_Y"])
# OK!

# Wrt y
plot(samp[,"Disc_New_y"] ~ samp[,"Disc_y"])
abline(0, 1, col="red")
mean(samp[, "Disc_New_y"] > samp[, "Disc_y"])
# Now OK
# Before adding random effects: [1] 0.03423333. This means that there are more 
# variation in y[] than expected - probably because observers vary in focus from
# picture to picture
# 
# 2. Quantiles for each variable:
# 
# 2.5%     25%     50%     75%    97.5%
# Disc_New_Y  0.14707  0.7933  1.4121  2.1963   4.1576
# Disc_New_y 25.60289 38.2987 47.8034 60.2280 103.6541
# Disc_Y      0.14868  1.0016  1.6883  2.6257   4.9057
# Disc_y     28.16323 41.1363 50.7363 63.4592 107.5061
# beta        0.04675  0.3554  0.5168  0.6756   0.9764
# mu0[1]     -2.38229 -1.9036 -1.5404 -1.1618  -0.4973
# mu0[2]     -2.81964 -2.4402 -2.1251 -1.7892  -1.1899
# mu0[3]     -6.80932 -6.0006 -5.0581 -3.9832  -2.5831
# mu0[4]     -6.83190 -6.1345 -5.2787 -4.2994  -2.7585
# mu0[5]     -4.21993 -3.8896 -3.5496 -3.1393  -2.3198
# mu0[6]     -3.31674 -2.9519 -2.6638 -2.3719  -1.8336
# mu0[7]     -6.84318 -6.2259 -5.4742 -4.5561  -2.9776
# mu0[8]     -4.54373 -4.2160 -3.8433 -3.3806  -2.4791
# mu_p1       0.63152  1.0096  1.2063  1.3976   1.7900
# mu_p2       1.08607  1.4835  1.6970  1.9202   2.3410
# sigma       1.31758  1.5985  1.7541  1.9110   2.2320
# sigma_p     0.69672  1.0833  1.2651  1.3916   1.4892
