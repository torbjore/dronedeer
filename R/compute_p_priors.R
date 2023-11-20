load(file = "data/posterior_samples/posterior_samples_fence_p.RData") # Fra 'enclosure_model.R'
summary(posterior_samples_fence_p)

# Creating a vector with p's from all 3 chains
posterior_samples_p_fence <- c(posterior_samples_fence_p$chain1[,"p[2]"],posterior_samples_fence_p$chain2[,"p[2]"], posterior_samples_fence_p$chain3[,"p[2]"])

# Taking logit of p's
logit <- function(p) log(p/(1-p)) 

logit_posterior_samples_p_fence <- logit(posterior_samples_p_fence)

# Finding logit means and SDs
mean_logit_posterior_samples_p_fence <- mean(logit_posterior_samples_p_fence)

sd_logit_posterior_samples_p_fence <- sd(logit_posterior_samples_p_fence)

# Changing prior names
mu_logit_p <- mean_logit_posterior_samples_p_fence

sigma_logit_p <- 1.2*sd_logit_posterior_samples_p_fence

# Saving prior parameters
prior_parameters_for_p <- list(
  mu_logit_p = mu_logit_p,
  sigma_logit_p = sigma_logit_p
)

save(prior_parameters_for_p, file = "data/prior_parameters_for_p.rda")
