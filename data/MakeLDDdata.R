# Making 'Lærdal Deer Drone data' (structured for nimble input)

# Loading site-data
load(file = "data/UseData.rda")

# Loading counts
load(file = "data/Counts.rda")

# KOMMENTAR: UseData inneholder også counts (men i et annet format)

# Number of sites for each survey
N_sites <- table(UseData$Survey)

# Making y, Y, area and mean_field_dist 
y <- tapply(Counts$Count, list(Counts$Survey, Counts$Site, Counts$Cat), sum) # 3D array rekke = survey, colonne = site (bilde), lag = Cat
Y <- apply(y, c(1,2), sum) # Matrix with survey x no. deer
area <- tapply(UseData$Area, list(UseData$Survey, UseData$Site), sum)/10000 # Area in hectare
mean_field_dist <- tapply(UseData$mean_field_dist, list(UseData$Survey, UseData$Site), sum)

# Making index for 'study-area and month' (sam)
Counts$StudyArea_Month <- paste(Counts$Study_area, Counts$Month, sep="_")
tab <- table(Counts$Survey,Counts$StudyArea_Month)
sam <- apply(tab,1, function(i) which(i != 0))

####################################
# Simple estimates for prior means #
####################################

# Y per survey_area-and-month (sam)
Y_per_sam <- apply(Y, 1, sum, na.rm=TRUE) |> tapply(INDEX = sam, sum)
sam_names <- c("Haugen (April)", "Haugen (March)", "Raa (April)", "Raa (March)", "Søre Bjørkum (April)", "Søre Bjørkum (March)", "Sprakehaug (April)", "Sprakehaug (March)")
names(Y_per_sam) <- sam_names

# Simple estimate of detection probability
colsumy <- apply(y, 3, sum, na.rm=TRUE)
p1hat <- colsumy[3]/(colsumy[2]+colsumy[3])
p2hat <- colsumy[3]/(colsumy[1]+colsumy[3])
(phat_simple <- 1 - (1-p1hat)*(1-p2hat))

# Area per sam
area_per_sam <- apply(area, 1, sum, na.rm=TRUE) |> tapply(INDEX = sam, sum)
(names(area_per_sam) <- sam_names)

# Simple estimates of density
lambdahat <- (Y_per_sam/phat_simple)/area_per_sam

# For the sam's with zero detentions (see prior_zero_sites.R)
lambdahat[3] <- 0.0066 # Raa (April)
lambdahat[4] <- 0.0046 # Raa (March)
lambdahat[7] <- 0.0047 # Sprakehaug (April)

LDDdata <- list(
  data = list(
    y = y,
    Y = Y
  ),
  const = list(
    N_sites = N_sites,
    area = area,
    mean_field_dist = mean_field_dist,
    sam = sam, # study area and month
    lambdahat = lambdahat
  )
)

save(LDDdata, file = "data/LDDdata.rda")
