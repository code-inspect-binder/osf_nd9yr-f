# -------------------------------------------------------------
# ordinal regression analyses
# run models for all questions programmatically

rm(list=ls())
setwd("/mnt/sda2/matteoHDD/git_local_HDD/covid19-misinformation/")

library(tidyverse)
library(rstan)
options(mc.cores = parallel::detectCores())

# load data
d <- read.csv("./data/data_allfmt.csv",header=T ,sep="\t")
codes <- read.csv2("./data/code_book.csv")

# -------------------------------------------------------------
# run models

d_all <- d # bk all data before subsetting for individual analyses
rm(d)

# dependent variables
ordinal_y <- c("att_pubs","att_schools","att_social","behav_meetings","mask_stores","mask_transport","vaccine","vaccine2")

## run sampling
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

for(var_y in ordinal_y){
  
  # create data list for modelling in Stan based on chosen dependent variable
  y0 <- d_all[,var_y]
  d <- d_all[!is.na(y0),]
  N_cat <- max(d[,var_y])
  
  d_stan <- list(
    y = d[,var_y],
    d_prime = d$d_cov,
    X = log(d$M_cov),
    K = N_cat,
    N = nrow(d),
    age = d$age,
    gender = ifelse(d$profile_gender==2, -0.5, 0.5), # 2=female
    grade = d$profile_socialgrade_cie,
    marital_stat = d$profile_marital_stat,
    income = d$income,
    region = d$profile_GOR,
    covid_affected = ifelse(d$covid_affected==1,0.5,-0.5),
    edu = d$education,
    vote2019 = d$vote2019,
    voteEUref = d$voteEURef,
    n_income = max(d$income),
    n_marital_stat = max(d$profile_marital_stat),
    n_grade = max(d$profile_socialgrade_cie),
    n_region = max(d$profile_GOR),
    n_edu = max(d$education),
    n_vote2019 = max(d$vote2019),
    n_voteEUref = max(d$voteEURef)
  )
  # str(d_stan)
  
  m01 <- stan(file = "./stan/ordinal_model.stan", data = d_stan, 
              iter = 4000,  
              warmup = 2000, 
              chains = 4, 
              init_r=0.25, 
              control=list(adapt_delta = 0.99))
  
  saveRDS(m01, sprintf("./regression_results/ordinal_%s.RDS",var_y))
  
  print(m01, pars=c("b_d","b_X","b_age","b_age2","b_covid_affected","b_gender","sigma_grade","sigma_vote2019","sigma_voteEUref","sigma_edu","sigma_marital_stat","sigma_income","c"),probs=c(0.025, 0.975),digits=3)
  
}
