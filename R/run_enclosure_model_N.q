# RUN MODEL FOR FENCED DEER

library(nimble)

# SOURCE THE MODEL CODE
source("R/nimble_models/enclosure_model_N.q")

# Loading the data
load("data/CountData_Fence.rda")

UseData_fence = CountData_Fence[!is.na(CountData_Fence$Area),]
UseData_fence$Survey = factor(UseData_fence$Survey, levels = unique(UseData_fence$Survey)) # sorted as in data

N_sites = table(UseData_fence$Survey)
UseData_fence$Site = unlist(lapply(N_sites, function(i) 1:i))

Counts_fence = rbind(
  cbind(UseData_fence[,1:4], Site = UseData_fence$Site, Cat = "only_ab", Count = UseData_fence$n_only_ab),
  cbind(UseData_fence[,1:4], Site = UseData_fence$Site, Cat = "only_jb", Count = UseData_fence$n_only_jb),
  cbind(UseData_fence[,1:4], Site = UseData_fence$Site, Cat = "both",    Count = UseData_fence$n_both)
)
Counts_fence$Cat = factor(Counts_fence$Cat, levels = c("only_ab", "only_jb", "both"))

y = tapply(Counts_fence$Count, list(Counts_fence$Survey, Counts_fence$Site, Counts_fence$Cat), sum)
Y = apply(y, c(1,2), sum)
area = tapply(UseData_fence$Area, list(UseData_fence$Survey, UseData_fence$Site), sum)/10000 # Area in hectare


# Function for initial values
colsumy = apply(y, 3, sum, na.rm=TRUE)
p1hat = colsumy[3]/(colsumy[2]+colsumy[3])
p2hat = colsumy[3]/(colsumy[1]+colsumy[3])
phat = (p1hat+p2hat)/2
#Nhat = apply(Y, 1, sum, na.rm=TRUE)/(1-(1-p1hat)*(1-p2hat)) # For all sites combined
#lambdahat = (Nhat/N_sites)/apply(area, 1, mean, na.rm=TRUE) + 0.01 # Adding a small value since we get -Inf from log(lambdahat=0)

Inits = function(){
  N = round(Y/(1-(1-p1hat)*(1-p2hat)), 0)
  #N[is.na(N)] = 0
  #sigma = runif(2, 0, 0.1)
  list(
    #mean_lambda = mean(lambdahat)*runif(length(mean(lambdahat)), 0.9, 1.1),
    #epsilon = matrix(rnorm(nrow(Y)*ncol(Y), 0, sigma) , nrow = nrow(Y), ncol = ncol(Y)),
    p1 = exp(log(phat)*runif(2, 0.9, 1.1)),
    p2 = exp(log(phat)*runif(2, 0.9, 1.1)),
    N = N,
    lambda_N = N #mean(N, na.rm = TRUE)
  )
}

# Setting up the MCMC
DoubleObsMultisiteModel <- nimbleModel(
  DoubleObsMultisiteCode_fence,
  constants = list(N_surv = length(N_sites), N_sites = N_sites, area = area), # 0.1 to 10 times point estimate
  data = list(y=y, Y=Y), # area = area),
  inits = Inits()
)

CDoubleObsMultisiteModel <- compileNimble(DoubleObsMultisiteModel)
DoubleObsMultisiteConf <- configureMCMC(DoubleObsMultisiteModel, monitors = c("p1", "p2", "N_tot", "Dens", "Dens_tot"))
DoubleObsMultisiteMCMC <- buildMCMC(DoubleObsMultisiteConf)
CDoubleObsMultisiteMCMC <- compileNimble(DoubleObsMultisiteMCMC)

posterior_samples_fence_p <- runMCMC(
  CDoubleObsMultisiteMCMC,
  niter=5000, #00,
  nburnin=1000, #00,
  nchain=3,
  thin = 2,
  inits = Inits,
  samplesAsCodaMCMC = TRUE)

library(coda)
plot(posterior_samples_fence_p)

summary(posterior_samples_fence_p)

#save(posterior_samples_fence_p, file = "posterior_samples_fence_p.RData")