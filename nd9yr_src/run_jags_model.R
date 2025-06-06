# run jags model

rm(list=ls())
setwd("/mnt/sda2/matteoHDD/git_local_HDD/covid19-misinformation/jags/")

# load data
d_jags <- readRDS("../data/data_list.RDS")
str(d_jags)
# 
# metad_2conditions <- function (nR_S1_tot, nR_S2_tot) {
#   # nR_S1_tot and nR_S2_tot should be lists of size N subjects.
#   # Each list contain 2 responses vectors.
#   # Design matrix: Condition 1: (0) - Condition 2: (1)

nR_S1 <- list()
nR_S2 <- list()

for(n in 1:d_jags$J){
  # nR_S1_tot[[n]][i]
  # nR_S2_tot[[n]][i]
  nR_S1[[n]] <- list(unname(d_jags$F1[n,]), unname(d_jags$F2[n,]))
  nR_S2[[n]] <- list(unname(d_jags$T1[n,]), unname(d_jags$T2[n,]))
}

rm(d_jags)

# load custom functions
source("./wrapper_2conditions.R")



# Fit all data at once

output <- metad_2conditions_v4(nR_S1 = nR_S1, nR_S2 = nR_S2)
saveRDS(output, "hdmeta_jags_both_MLok.RDS")

