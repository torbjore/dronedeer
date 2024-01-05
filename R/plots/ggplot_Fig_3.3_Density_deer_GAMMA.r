# GGPLOT FIG 3.3 DENSITY OF DEER AT SAMPLING AREAS AND MONTH, GAMMA DISTIBUTION

library(coda)
library(ggplot2)
library(hrbrthemes) 

library(cowplot)

# Loading site-data
load(file = "data/derived/UseData.rda")

# Loading counts
load(file = "data/derived/Counts.rda")

#load("data/CountData.RData")
sam_names = c("Haugen (April)", "Haugen (March)", "Raa (April)", "Raa (March)", "Soere Bjoerkum (April)", "Soere Bjoerkum (March)", "Sprakehaug (April)", "Sprakehaug (March)")

# Field
load("original/posterior_gamma_2.RData")
predictor_variable = UseData$mean_field_dist
Xlab = "Distance from field (m)"
Ylab = expression(Density~of~deer~(km^-2))
modelcov = "field_plot_"

#summary(posterior_gamma_2$samples)
post_samp_mat = as.matrix(posterior_gamma_2$samples)

X = seq(0, max(predictor_variable, na.rm = TRUE), length.out=50)
X_mean <-  mean(predictor_variable, na.rm=TRUE)
X_sd = sd(predictor_variable, na.rm=TRUE)
X_st = (X - X_mean)/X_sd

mpv <- mean(predictor_variable, na.rm = TRUE) # mean predictor variable
pap <- mean(post_samp_mat[,"x_at_peak"])*X_sd + X_mean# predictor variable at peak


pred_gamma2_orig = function(x, PS = post_samp_mat, sam){
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
  return(mean_lambda_km2)
}

pred_plots <- function(predfun, x=X, x_st=X_st, PS = post_samp_mat, names = sam_names){
  plots = list() # Empty list to store plots in
  
  for(SAM in 1:length(sam_names)){
    Title = names[SAM]
    cat(Title, "\n")
    df = data.frame(x = x, x_st = x_st)
    for(i in 1:length(x_st)){
      y = predfun(x = x_st[i], sam=SAM)
      Q = quantile(y, probs = c(0.025, 0.5, 0.975))
      df$lwr[i] = Q[1]
      df$median[i] = Q[2]
      df$upr[i] = Q[3]
      df$mean[i] = mean(y)
    }
    
    plots[[SAM]] = ggplot(data = df, aes(x=x)) + 
      # Adding CI, mean and median
      #geom_ribbon(aes(ymin=lwr, ymax=upr, x=x, fill ="Lower/Upper 95% CI", alpha = 0.1), show.legend = TRUE) +
      geom_ribbon(aes(ymin=lwr, ymax=upr, x=x, fill ="Lower/Upper 95% CI", alpha = 0.1), show.legend = FALSE) +
      geom_line(aes(y=mean, color="Posterior mean"), color="#fe2323", linewidth=1.2) + 
      geom_line(aes(y=median, color="Posterior median"), color="#369bbe", linewidth=1.2) +
      geom_vline(aes(xintercept=mpv), color="#68754D", linetype="dashed", linewidth=1) +
      geom_vline(aes(xintercept = pap), color="#68754D", linetype="dotted", linewidth=1 ) +
      
      # Adding legend
      # scale_fill_manual(name = "",
      #                   values= c("Posterior mean" = "#fe2323",
      #                             "Posterior median" = "#369bbe",
      #                             "Lower/Upper 95% CI"= "#d9f5e6")) +
      #guides(alpha = FALSE) + # Remove alpha from legend'
      
      # Add title and axis labels
      ggtitle(Title) +
      xlab(Xlab) +
      ylab(Ylab) +
      scale_y_continuous(trans='log10') +
      theme_ipsum(axis_title_size = 11.5, axis_title_just = "cc",)
      #theme_ipsum(axis_title_size = 8, axis_title_just = "cc",)
    
    # print(plots[[SAM]])
    # png(paste(modelcov, SAM, ".png", sep = ""), width=600, height=500, res=120) # Start export
    # print(plots[[SAM]])
    # dev.off() # End export
  }
  return(plots)
}

plots <- pred_plots(pred_gamma2_orig)

plot_grid(plots[[1]],
          plots[[2]],
          plots[[3]],
          plots[[4]])
plot_grid(plots[[5]],
          plots[[6]],
          plots[[7]],
          plots[[8]])
