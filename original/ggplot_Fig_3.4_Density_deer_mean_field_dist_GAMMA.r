# NOT MEAN DISTANCE FROM FIELD BUT DISTANCE FROM FIELD WITH PEAK DENSITY


setwd("C:/Users/Bruker/Dropbox (UiO)/Deer/Analysis/Plots")

library(coda)
library(ggplot2)
library(tidyverse)

load("../Make_CountData/CountData.RData")
sam_names = c("Haugen April", "Haugen March", "Raa April", "Raa March", "Soere Bjoerkum April", "Soere Bjoerkum March", "Sprakehaug April", "Sprakehaug March")
CountData$StudyArea_Month = paste(CountData$Study_area, CountData$Month, sep="_")

# Field
load("../Model/Manuscript_models/posterior_gamma_2.RData")
predictor_variable = CountData$mean_field_dist
Title = "Peak density of deer"
Xlab = "Mean elevation (MASL)"
Ylab = expression(Density~of~deer~(km^-2))
modelcov = "field_plot_"

#summary(posterior_samples)
post_samp_mat = as.matrix(posterior_gamma_2$samples)

x = mean(predictor_variable, na.rm = TRUE)
x_at_peak = mean(- unlist(posterior_gamma_2$samples[,"beta[1]"])/(2*unlist(posterior_gamma_2$samples[,"beta[2]"])))

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


df = data.frame(x=rep(x, length(sam_names)), 
                x_at_peak=rep(x_at_peak, length(sam_names)), 
                lwr=NA, 
                upr=NA, 
                median=NA, 
                mean=NA, 
                mean_elevation = tapply(CountData$mean_elevation, CountData$StudyArea_Month, mean)[-1],
                Legend = sam_names)

for(SAM in 1:length(sam_names)){
  y = FUN_post_samp_mean_lambda_x(x = x_at_peak, sam=SAM)
  Q = quantile(y, probs = c(0.025, 0.5, 0.975))
  df[SAM,"lwr"] = Q[1]
  df[SAM, "median"] = Q[2]
  df[SAM, "upr"] = Q[3]
  df[SAM, "mean"] = mean(y)
}


df = rownames_to_column(df, var = "Sampling area & month")
df

write.table(df)

#require(scales)
all <- ggplot(df, aes(x=mean_elevation, y=mean, color = Legend, shape = Legend)) +
  
  # CI:
  geom_pointrange(aes(ymin=lwr, ymax=upr), size = 0.6) +
  
  # Title and axis:
  ggtitle(Title) +
  xlab(Xlab) +
  ylab(Ylab) +
  scale_y_continuous(trans='log10') +
  theme_ipsum(axis_title_size = 11.5, axis_title_just = "cc",) +
  
  scale_shape_manual(values = c(16, 17, 16, 17, 16, 17, 16, 17)) +
  
  #Legend:
  
  scale_color_manual(name = "",
                     values= c("Raa April" = "#fe2323", 
                               "Raa March" = "#fe2323",
                               "Soere Bjoerkum April" = "#369bbe", 
                               "Soere Bjoerkum March" = "#369bbe",
                               "Haugen April" = "#D0EDBF",
                               "Haugen March" = "#D0EDBF",
                               "Sprakehaug April" = "#4e4e4e",
                               "Sprakehaug March" = "#4e4e4e")) +
  guides(shape = FALSE,
         colour = guide_legend(override.aes = list(shape = c(16, 17, 16, 17, 16, 17, 16, 17),
                                                   size = 0.8,
                                                   linetype = c("blank", "blank", "blank", "blank", "blank", "blank", "blank", "blank")))) +
  theme(legend.title = element_blank()) # Remove legend title

all

png(paste(modelcov, "all.png", sep = ""), width=600, height=500, res=120)
print(all)
dev.off()