# GGPLOT FIG 3.3 DENSITY OF DEER AT SAMPLING AREAS AND MONTH

library(coda)
library(ggplot2)
library(hrbrthemes) 
library(cowplot)

#############
# Functions #
#############

source("R/pred_plot_functions.r")

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
pap <- median(-post_samp_mat[,"beta[1]"]/(2*post_samp_mat[,"beta[2]"]))*X_sd + X_mean # predictor variable at peak

############
# Plotting #
############

plots <- pred_plots(PS = post_samp_mat)

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

#plot(out$samples) # 1 = black, 2 = red, 3 = green
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
#abline(0, 1, col="red")
mean(samp[, "Disc_New_Y"] > samp[, "Disc_Y"])
# OK!

# Wrt y
#plot(samp[,"Disc_New_y"] ~ samp[,"Disc_y"])
#abline(0, 1, col="red")
mean(samp[, "Disc_New_y"] > samp[, "Disc_y"])
