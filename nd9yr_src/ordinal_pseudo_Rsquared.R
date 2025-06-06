# Calculate pseudo R-squared for ordinal regression models.
# Each parameter is set to their posterior mean as point estimate.

rm(list=ls())
setwd("/mnt/sda2/matteoHDD/git_local_HDD/covid19-misinformation/control_analyses")

#
library(rstan)
library(tidyverse)
library(tidybayes)
source("../functions_bayes_ml.R")
source("../plot_wrapper_c01.R")

d <- read_delim("../data/data_allfmt.csv")
ordinal_y <- c("att_pubs","att_schools","att_social","behav_meetings","mask_stores","mask_transport","vaccine2")


## helper functions

# likelihood of null model
# r <- d$vaccine2
LLnull <- function(r){
  r <- r[which(!is.na(r))]
  p <- tapply(r,r,length)/length(r)
  return(sum(log(p[r])))
}

# ordered logistic
ordered_logistic <- function(eta, cutpoints){
  cutpoints <- c(cutpoints, Inf)
  k <- length(cutpoints)
  p <- rep(NA, k)
  p[1] <- plogis(cutpoints[1], location=eta, scale=1, lower.tail=TRUE)
  for(i in 2:k){
    p[i] <- plogis(cutpoints[i], location=eta, scale=1, lower.tail=TRUE) - 
      plogis(cutpoints[i-1], location=eta, scale=1, lower.tail=TRUE)
  }
  return(p)
}

# 
LLmodel <- function(r,m, sd2=T){
  
  index <- which(!is.na(r))
  d <- read_delim("../data/data_allfmt.csv")
  r <- r[index]
  d <- d[index,]
 
   m %>%
    spread_draws(c[cutpoint])%>%
    count(cutpoint) %>% nrow() -> n_cutpoints
  
  m %>%
    spread_draws(c[cutpoint]) %>%
    pivot_wider(names_from=cutpoint, values_from=c) %>%
    select(num_range("", 1:n_cutpoints))  %>%
    summarize_all(mean) -> cutpoints
    
  m %>%
    spread_draws(b_d, b_X, b_age, b_age2, b_covid_affected, b_gender) %>%
    select(starts_with("b_")) %>%
    summarise_all(mean)-> fixef
  
  m %>%
    spread_draws(sigma_grade,sigma_vote2019,sigma_voteEUref,sigma_edu,sigma_region, sigma_marital_stat, sigma_income) %>%
    select(starts_with("sigma_")) %>%
    summarise_all(mean)-> sigmas
  
  m %>%
    spread_draws(b_grade[grade]) %>%
    pivot_wider(names_from=grade, values_from=b_grade) %>%
    select(matches('\\d{1,}'))  %>%
    summarise_all(mean) -> b_grade
  
  m %>%
    spread_draws(b_vote2019[vote]) %>%
    pivot_wider(names_from=vote, values_from=b_vote2019) %>%
    select(matches('\\d{1,}'))  %>%
    summarise_all(mean)-> b_vote2019
  
  m %>%
    spread_draws(b_voteEUref[vote]) %>%
    pivot_wider(names_from=vote, values_from=b_voteEUref) %>%
    select(matches('\\d{1,}'))  %>%
    summarise_all(mean)-> b_voteEUref
  
  m %>%
    spread_draws(b_region[region]) %>%
    pivot_wider(names_from=region, values_from=b_region) %>%
    select(matches('\\d{1,}'))  %>%
    summarise_all(mean)-> b_region
  
  m %>%
    spread_draws(b_edu[edu]) %>%
    pivot_wider(names_from=edu, values_from=b_edu) %>%
    select(matches('\\d{1,}'))  %>%
    summarise_all(mean)-> b_edu
  
  m %>%
    spread_draws(b_marital_stat[marital_stat]) %>%
    pivot_wider(names_from=marital_stat, values_from=b_marital_stat) %>%
    select(matches('\\d{1,}'))  %>%
    summarise_all(mean)-> b_marital_stat
  
  m %>%
    spread_draws(b_income[income]) %>%
    pivot_wider(names_from=income, values_from=b_income) %>%
    select(matches('\\d{1,}'))  %>%
    summarise_all(mean)-> b_income
  
  # rescaling for continuous covariates?
  if(sd2==T){
    age_c <- (d$age - mean(d$age))/(2*sd(d$age))
    X <- log(d$M_cov)
    X <- (X - mean(X))/(2*sd(X))
    d_prime <- (d$d_cov - mean(d$d_cov))/(2*sd(d$d_cov))
  }else{
    age_c <- (d$age - mean(d$age))/sd(d$age)
    X <- log(d$M_cov)
    d_prime <- (d$d_cov - mean(d$d_cov))/(sd(d$d_cov))
  }
  
  gender <- ifelse(d$profile_gender==2, -0.5, 0.5) # 2=female
  pred_p <- rep(NA, nrow(d))
  
  for(i in 1:nrow(d)){
    gamma <- unlist(unname(fixef['b_d'] * d_prime[i] +
                             fixef['b_X'] * X[i] +
                             fixef['b_gender'] * gender[i] +
                             fixef['b_covid_affected']*d$covid_affected[i] +
                             fixef['b_age'] * age_c[i] + fixef['b_age2'] * age_c[i]^2 +
                             b_grade[as.character(d$profile_socialgrade_cie[i])]*sigmas['sigma_grade'] + 
                             b_region[as.character(d$profile_GOR[i])]*sigmas['sigma_region'] + 
                             b_edu[,as.character(d$education[i])]*sigmas[,'sigma_edu'] + 
                             b_vote2019[as.character(d$vote2019[i])]*sigmas['sigma_vote2019'] + 
                             b_voteEUref[as.character(d$voteEURef[i])]*sigmas['sigma_voteEUref'] +
                             b_marital_stat[as.character(d$profile_marital_stat[i])]*sigmas['sigma_marital_stat']+ 
                             b_income[as.character(d$income[i])]*sigmas['sigma_income']))
    
    pred_p[i] <- ordered_logistic(gamma, unlist(cutpoints))[r[i]]
  }
  
  #
  return(sum(log(pred_p)))
}

# pseudo_r2
pseudo_R2 <- function(r,m, sd2=T){
  LL_null <- LLnull(r)
  LL_model <- LLmodel(r,m)
  return(1 - LL_model/LL_null)
}

## calcualte pseudo R2
R2 <- rep(NA, length(ordinal_y))
for (QN in 1:length(ordinal_y)){
  type <- 2
  # if(type==1){
  #   m <- readRDS(sprintf("../control_analyses/ordinal_extended_control01_norm_%s.RDS",ordinal_y[QN]))
  # }else{
    m <- readRDS(sprintf("../control_analyses/ordinal_extended_control0%i_norm_%s.RDS",type, ordinal_y[QN]))
  # }
  
  r <- unlist(d[,ordinal_y[QN]])
  
  
  R2[QN] <- pseudo_R2(r,m)
  
}
names(R2) <- ordinal_y
data.frame(R2)
# R2
# att_pubs       0.13112820
# att_schools    0.07906952
# att_social     0.11457656
# behav_meetings 0.04995401
# mask_stores    0.08410775
# mask_transport 0.06869228
# vaccine2       0.21256426

range(R2) # 0.04995401 0.21256426
mean(R2) # 0.1057275
