# Making 'Lærdal Deer Drone data' (structured for nimble input)

# Loading site-data
load(file = "data/derived/UseData.rda")

# Loading counts
load(file = "data/derived/Counts.rda")

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

LDDdata <- list(
  data = list(
    y = y,
    Y = Y
  ),
  const = list(
    N_sites = N_sites,
    area = area,
    mean_field_dist = mean_field_dist,
    sam = sam # study area and month
  )
)

save(LDDdata, file = "data/LDDdata.rda")
