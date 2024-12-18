# RUN MODEL FOR FENCED DEER

library(nimble)

# SOURCE THE MODEL CODE
source("R/nimble_models/enclosure_model_N.q")

# Loading the data
load("data/nimbleData_Fence.rda")

# adding the total area covered
nimbleData_Fence$constants$sum_area <- apply(nimbleData_Fence$constants$area, 1, sum, na.rm = TRUE)

# Function for initial values
colsumy = apply(nimbleData_Fence$data$y, 3, sum, na.rm=TRUE)
p1hat = colsumy[3]/(colsumy[2]+colsumy[3])
p2hat = colsumy[3]/(colsumy[1]+colsumy[3])
phat = (p1hat+p2hat)/2
#Nhat = apply(Y, 1, sum, na.rm=TRUE)/(1-(1-p1hat)*(1-p2hat)) # For all sites combined
#lambdahat = (Nhat/N_sites)/apply(area, 1, mean, na.rm=TRUE) + 0.01 # Adding a small value since we get -Inf from log(lambdahat=0)

Inits = function(){
  N = round(nimbleData_Fence$data$Y/(1-(1-p1hat)*(1-p2hat)), 0)
  #N[is.na(N)] = 0
  #sigma = runif(2, 0, 0.1)
  list(
    #mean_lambda = mean(lambdahat)*runif(length(mean(lambdahat)), 0.9, 1.1),
    #epsilon = matrix(rnorm(nrow(Y)*ncol(Y), 0, sigma) , nrow = nrow(Y), ncol = ncol(Y)),
    #p = exp(log(phat)*runif(2, 0.9, 1.1)),
    mu_p = runif(2, -1, 1),
    sigma_p = runif(2, 0.1, 0.2),
    N = N,
    New_Y = nimbleData_Fence$data$Y, # Warning message if not included
    New_y = nimbleData_Fence$data$y,
    lambda_N = N #mean(N, na.rm = TRUE)
  )
}

# Setting up the MCMC
DoubleObsMultisiteModel <- nimbleModel(
  DoubleObsMultisiteCode_fence,
  constants = nimbleData_Fence$constants, #list(N_surv = length(N_sites), N_sites = N_sites, area = area), # 0.1 to 10 times point estimate
  data = nimbleData_Fence$data, #list(y=y, Y=Y), # area = area),
  inits = Inits()
)

CDoubleObsMultisiteModel <- compileNimble(DoubleObsMultisiteModel)
DoubleObsMultisiteConf <- configureMCMC(DoubleObsMultisiteModel, monitors = c("mu_p", "sigma_p", "N_sum", "Disc_New_Y", "Disc_Y", "Disc_New_y", "Disc_y", "N_tot", "var_N_tot", "N_lower", "N_upper"))
DoubleObsMultisiteMCMC <- buildMCMC(DoubleObsMultisiteConf)
CDoubleObsMultisiteMCMC <- compileNimble(DoubleObsMultisiteMCMC)

#inits.list <- list(Inits(), Inits(), Inits())

posterior_samples <- runMCMC(
  CDoubleObsMultisiteMCMC,
  niter = 20000,
  nburnin = 2000,
  nchain=3,
  thin = 2,
  inits = Inits, #inits.list,
  samplesAsCodaMCMC = TRUE)

library(coda)
#plot(posterior_samples)
#gelman.diag(posterior_samples)
summary(posterior_samples)

# Posterior predictive checks
# Wrt Y:
samp <- as.matrix(posterior_samples)
plot(samp[,"Disc_New_Y"] ~ samp[,"Disc_Y"])
abline(0, 1, col="red")
mean(samp[, "Disc_New_Y"] > samp[, "Disc_Y"])
# OK!

# Wrt y:
plot(samp[,"Disc_New_y"] ~ samp[,"Disc_y"])
abline(0, 1, col="red")
mean(samp[, "Disc_New_y"] > samp[, "Disc_y"])

save(posterior_samples, file = "data/posterior_samples/posterior_samples_enclosure_N.RData")