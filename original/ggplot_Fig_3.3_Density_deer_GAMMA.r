# GGPLOT FIG 3.3 DENSITY OF DEER AT SAMPLING AREAS AND MONTH, GAMMA DISTIBUTION

setwd("C:/Users/Bruker/Dropbox (UiO)/Deer/Analysis/Plots")

library(coda)
library(ggplot2)
library(hrbrthemes) 

load("../Make_CountData/CountData.RData")
sam_names = c("Haugen (April)", "Haugen (March)", "Raa (April)", "Raa (March)", "Soere Bjoerkum (April)", "Soere Bjoerkum (March)", "Sprakehaug (April)", "Sprakehaug (March)")

# Field
load("../Model/Manuscript_models/posterior_gamma_2.RData")
predictor_variable = CountData$mean_field_dist
Xlab = "Distance from field (m)"
Ylab = expression(Density~of~deer~(km^-2))
modelcov = "field_plot_"

#summary(posterior_gamma_2$samples)
post_samp_mat = as.matrix(posterior_gamma_2$samples)

x = seq(0, max(predictor_variable, na.rm = TRUE), length.out=50)
x_st = (x - mean(predictor_variable, na.rm=TRUE))/sd(predictor_variable, na.rm=TRUE)

line = mean(predictor_variable, na.rm = TRUE)


FUN_post_samp_mean_lambda_x = function(x, PS = post_samp_mat, sam){
  mu0 = PS[, paste0("mu0[", sam, "]")]
  beta1 = PS[, "beta[1]"]
  beta2 = PS[, "beta[2]"]
  sigma = PS[, "sigma"]
  
  mean_lambda = exp(mu0)
  mean_lnorm = log(mean_lambda) - 0.5*sigma*sigma
  var_lognormal = (exp(sigma*sigma) - 1)*exp(2*mean_lnorm + sigma*sigma)  
  rate = mean_lambda/var_lognormal  
  
  mu = mu0 + x*beta1 + (x^2)*beta2
  shape = exp(mu)*rate
  
  mean_lambda_km2 = (shape/rate)*100 
  return(mean_lambda_km2 = mean_lambda_km2)
}

plots = list() # Empty list to store plots in

for(SAM in 1:length(sam_names)){
  Title = sam_names[SAM]
  cat(Title, "\n")
  df = data.frame(x = x, x_st = x_st)
  for(i in 1:length(x_st)){
    y = FUN_post_samp_mean_lambda_x(x = x_st[i], sam=SAM)
    Q = quantile(y, probs = c(0.025, 0.5, 0.975))
    df$lwr[i] = Q[1]
    df$median[i] = Q[2]
    df$upr[i] = Q[3]
    df$mean[i] = mean(y)
  }
  
  plots[[SAM]] = ggplot(data = df, aes(x=x)) + 
    
    # Adding CI, mean and median
    geom_ribbon(aes(ymin=lwr, ymax=upr, x=x, fill ="Lower/Upper 95% CI", alpha = 0.1), show.legend = TRUE) +
    geom_line(aes(y=mean, color="Posterior mean"), color="#fe2323", size=1.2) + 
    geom_line(aes(y=median, color="Posterior median"), color="#369bbe", size=1.2) +
    geom_vline(aes(xintercept=line), color="#68754D", linetype="dashed", size=1) +
    geom_vline(aes(xintercept = 277.5447), color="#68754D", linetype="dotted", size=1 ) +
    
    # Adding legend
    scale_fill_manual(name = "",
                      values= c("Posterior mean" = "#fe2323",
                                "Posterior median" = "#369bbe",
                                "Lower/Upper 95% CI"= "#d9f5e6")) +
    guides(alpha = FALSE) + # Remove alpha from legend'
    
    # Add title and axis labels
    ggtitle(Title) +
    xlab(Xlab) +
    ylab(Ylab) +
    scale_y_continuous(trans='log10') +
    theme_ipsum(axis_title_size = 11.5, axis_title_just = "cc",)
  
  print(plots[[SAM]])
  png(paste(modelcov, SAM, ".png", sep = ""), width=600, height=500, res=120) # Start export
  print(plots[[SAM]])
  dev.off() # End export
}
