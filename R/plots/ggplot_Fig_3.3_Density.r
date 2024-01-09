# GGPLOT FIG 3.3 DENSITY OF DEER AT SAMPLING AREAS AND MONTH

library(coda)
library(ggplot2)
library(hrbrthemes) 
library(cowplot)

#############
# Functions #
#############

pred <- function(x, PS = post_samp_mat, sam){ #general prediction function
  
  parnames <- dimnames(PS)[[2]]
  
  mu0 <- PS[, paste0("mu0[", sam, "]")]
  
  if("beta[2]" %in% parnames){
    beta1 <- PS[, "beta[1]"]
    beta2 <- PS[, "beta[2]"]
  } else {
    beta1 <- PS[, "beta"]
    beta2 <- 0
  }
  
  if("sigma" %in% parnames){
    sigma <- PS[, "sigma"]
  } else {
    sigma <- 0
  }
  
  mu <- mu0 + x*beta1 + (x^2)*beta2

  mean_lambda_km2 <- exp(mu + 0.5*sigma^2)*100 
  return(mean_lambda_km2)
}


pred_plots <- function(predfun = pred, x=X, x_st=X_st, PS = post_samp_mat, names = sam_names){
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
      #geom_line(aes(y=mean, color="Posterior mean"), color="#fe2323", linewidth=1.2) + 
      geom_line(aes(y=median, color="Posterior median"), color="#369bbe", linewidth=1.2) +
      geom_vline(aes(xintercept=mpv), color="#68754D", linetype="dashed", linewidth=1) +
      #geom_vline(aes(xintercept = pap), color="#68754D", linetype="dotted", linewidth=1 ) +
      
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

#############
# Preparing #
#############

# Loading data and posterior samples
load("data/posterior_samples/gamma_2.RData")
load(file = "data/UseData.rda")
load(file = "data/Counts.rda")

sam_names <- c("Haugen (April)", "Haugen (March)", "Raa (April)", "Raa (March)", "Soere Bjoerkum (April)", "Soere Bjoerkum (March)", "Sprakehaug (April)", "Sprakehaug (March)")

predictor_variable <- UseData$mean_field_dist
Xlab <- "Distance from field (m)"
Ylab <- expression(Density~of~deer~(km^-2))

post_samp_mat <- as.matrix(out$samples)

X <- seq(0, max(predictor_variable, na.rm = TRUE), length.out=50)
X_mean <-  mean(predictor_variable, na.rm=TRUE)
X_sd <- sd(predictor_variable, na.rm=TRUE)
X_st <- (X - X_mean)/X_sd

mpv <- mean(predictor_variable, na.rm = TRUE) # mean predictor variable
#pap <- mean(post_samp_mat[,"x_at_peak"])*X_sd + X_mean # predictor variable at peak

############
# Plotting #
############

plots <- pred_plots()

plot_grid(plots[[1]],
          plots[[2]],
          plots[[3]],
          plots[[4]])
plot_grid(plots[[5]],
          plots[[6]],
          plots[[7]],
          plots[[8]])

###########################
# Convergence diagnostics #
###########################

#plot(out$samples) # 1 = svart, 2 = rød, 3 = grønn
gelman.diag(out$samples)

#######################
# Posterior summaries #
#######################
summary(out$samples)
out$WAIC

###############################
# Posterior predictive checks #
###############################

samp <- as.matrix(out$samples)
#plot(samp[,"Disc_New_Y"] ~ samp[,"Disc_Y"])
abline(0, 1, col="red")
mean(samp[, "Disc_New_Y"] > samp[, "Disc_Y"])
# OK!

# Wrt y
#plot(samp[,"Disc_New_y"] ~ samp[,"Disc_y"])
abline(0, 1, col="red")
mean(samp[, "Disc_New_y"] > samp[, "Disc_y"])
