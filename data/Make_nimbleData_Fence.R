# Loading the data
load("data/CountData_Fence.rda")

UseData_fence <- CountData_Fence[CountData_Fence$Area_Fence > 0.001,]
UseData_fence$Survey <- factor(UseData_fence$Survey, levels = unique(UseData_fence$Survey)) # sorted as in data

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

nimbleData_Fence <- list(
  data = list(y=y, Y=Y),
  constants = list(N_surv = length(N_sites), N_sites = N_sites, area = area)
)

save(nimbleData_Fence, file = 'data/nimbleData_Fence.rda')

