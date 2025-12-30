###########################
# Posterior-prior overlap #
###########################

ddensity <- function(x, dens) approx(x=dens$x, y=dens$y, xout=x)$y

post_prior_overlap <- function(post_samp, prior_density_fun, plt = FALSE, plt.main = NULL, full = FALSE, ...){
  post_density <- density(post_samp, cut = 0) #, bw = "sj") # cut = 0 for truncated distrubutions
  prior_density <- prior_density_fun(post_density$x, ...)
  if(plt){
    plot(post_density, main = plt.main)
    lines(prior_density ~ post_density$x, col = "red")
  }
  dens <- list(
    x = post_density$x,
    y = pmin(prior_density, post_density$y)
  )
  int <- integrate(ddensity, min(post_density$x), max(post_density$x), dens = dens, stop.on.error = FALSE)
  if(full){
    return(int)
  } else {
    return(int$value)
  }
}

#####################################################
# Functions for predictions and plot of predictions #
#####################################################

pred <- function(x, PS, sam){ #general prediction function
  
  parnames <- dimnames(PS)[[2]]
  
  mu0 <- PS[, paste0("mu0[", sam, "]")]
  
  if("beta[2]" %in% parnames){
    beta1 <- PS[, "beta[1]"]
    beta2 <- PS[, "beta[2]"]
  } else if("beta[1]" %in% parnames){
    beta1 <- PS[, "beta"]
    beta2 <- 0
  } else {
    beta1 <- 0
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


pred_plots <- function(predfun = pred, x=X, x_st=X_st, PS, names = sam_names){
  plots = list() # Empty list to store plots in
  
  for(SAM in 1:length(sam_names)){
    Title = names[SAM]
    #cat(Title, "\n")
    df = data.frame(x = x, x_st = x_st)
    for(i in 1:length(x_st)){
      y = predfun(x = x_st[i], PS=PS, sam=SAM)
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
      # geom_vline(aes(xintercept = pap), color="#68754D", linetype="dotted", linewidth=1 ) +
      
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
      scale_y_continuous(trans='log10')
  }
  return(plots)
}