# RUN MODEL FOR FENCED DEER

library(nimble)

# SOURCE THE MODEL CODE
source("R/nimble_models/enclosure_model.q")

# Loading the data
load("data/nimbleData_Fence.rda")

# adding the total area covered
nimbleData_Fence$constants$sum_area <- apply(nimbleData_Fence$constants$area, 1, sum, na.rm = TRUE)

# Function for initial values
colsumy = apply(nimbleData_Fence$data$y, 3, sum, na.rm=TRUE)
p1hat = colsumy[3]/(colsumy[2]+colsumy[3])
p2hat = colsumy[3]/(colsumy[1]+colsumy[3])
phat = (p1hat+p2hat)/2

Inits = function(){
  N = round(nimbleData_Fence$data$Y/(1-(1-p1hat)*(1-p2hat)), 0)
  list(
    eta = runif(2, -1, 1),
    sigma_p = runif(1, 0.1, 0.2),
    N = N,
    New_Y = nimbleData_Fence$data$Y, # Warning message if not included
    New_y = nimbleData_Fence$data$y,
    lambda_N = N 
  )
}

# Setting up the MCMC
DoubleObsMultisiteModel <- nimbleModel(
  DoubleObsMultisiteCode_fence,
  constants = nimbleData_Fence$constants, 
  data = nimbleData_Fence$data,
  inits = Inits()
)

CDoubleObsMultisiteModel <- compileNimble(DoubleObsMultisiteModel)
DoubleObsMultisiteConf <- configureMCMC(DoubleObsMultisiteModel, monitors = c("eta", "sigma_p", "N_sum", "Disc_New_Y", "Disc_Y", "Disc_New_y", "Disc_y", "N_tot", "var_N_tot", "N_lower", "N_upper"))
DoubleObsMultisiteMCMC <- buildMCMC(DoubleObsMultisiteConf)
CDoubleObsMultisiteMCMC <- compileNimble(DoubleObsMultisiteMCMC)

#inits.list <- list(Inits(), Inits(), Inits())

posterior_samples <- runMCMC(
  CDoubleObsMultisiteMCMC,
  niter = 85000,
  nburnin = 5000,
  nchain=3,
  thin = 4,
  inits = Inits, #inits.list,
  samplesAsCodaMCMC = TRUE)

library(coda)
#plot(posterior_samples)
#gelman.diag(posterior_samples)
#summary(posterior_samples)

# Posterior predictive checks
# Wrt Y:
# samp <- as.matrix(posterior_samples)
# plot(samp[,"Disc_New_Y"] ~ samp[,"Disc_Y"])
# abline(0, 1, col="red")
# mean(samp[, "Disc_New_Y"] > samp[, "Disc_Y"])
# # OK!
# 
# # Wrt y:
# plot(samp[,"Disc_New_y"] ~ samp[,"Disc_y"])
# abline(0, 1, col="red")
# mean(samp[, "Disc_New_y"] > samp[, "Disc_y"])

save(posterior_samples, file = "data/posterior_samples/posterior_samples_enclosure_N.RData")

# Saving prior parameters for detection probabilities
#load("data/posterior_samples/posterior_samples_enclosure_N.RData")
eta <- as.matrix(posterior_samples)[,"eta[2]"]
prior_parameters_for_p <- list(
  mean_eta = mean(eta),
  sd_eta = sd(eta)
)

save(prior_parameters_for_p, file = "data/prior_parameters_for_p.rda")
