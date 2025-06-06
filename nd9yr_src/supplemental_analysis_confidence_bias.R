# here we use an ordinal model to compare whether overall respondents rated higher confidence 
# in COVID-19 or Science items, as a way to check whether they were influenced the 
# possibly larger epistemic uncertainty around COVID-19 facts

rm(list=ls())
setwd("/mnt/sda2/matteoHDD/git_local_HDD/covid19-misinformation/")

library(tidyverse)
library(ggplot2)

# -----------------------------------------
# get unsigned confidence 

d <- read.csv("./data/UoEssex_Results_210414_Client.csv")

d_sub <- d[,21:(21+27)] # subset columns with answers to science/covid statements
colnames(d)[21:(21+27)] # sanity check: all questions
colnames(d)[21:(21+13)] # science items
colnames(d)[35:(35+13)] # covid items

d_sub$id <- d$ID # participant identifier
d_long <- as.data.frame(pivot_longer(d_sub, cols=-id, names_to="question",values_to = "rating_raw"))

# categorize question labels
Q <- read_csv("./data/statement_label_truth.csv")

# recode responses in long database
d_long$respT <- ifelse(d_long$rating_raw<=3,1,0)
d_long$type <- ifelse(as.numeric(substr(d_long$question,7,8))<=14,"science","covid")
d_long$isTrue <- Q$true[as.numeric(substr(d_long$question,7,8))]
d_long$correct <- ifelse(d_long$respT==d_long$isTrue,1,0)

# reverse code confidence rating (so that goes confident False >>> confident True)
d_long$rating_rev <- 7 - d_long$rating_raw

# transform into pure confidence levels, regardless of response true/false
d_long$rating_abs <- ifelse(d_long$rating_rev>3, 
                            d_long$rating_rev-3,
                            4 - d_long$rating_rev)

tapply(d_long$rating_abs, d_long$type, mean)

# sanity
tapply(d_long$rating_abs, list(d_long$rating_abs, d_long$type,d_long$correct), length)
tapply(d_long$rating_abs, list(d_long$type,d_long$isTrue), mean)
# tapply(d_long$rating_abs, list(d_long$type,d_long$correct,d_long$isTrue), mean)

# sanity 2
d_long %>%
  group_by(id, type) %>%
  summarise(confidence = mean(rating_abs)) %>%
  wilcox.test(data=., confidence ~ type, paired=T)

# -----------------------------------------
# compare confidence in science vs Covid-19

library(rstan)
options(mc.cores = parallel::detectCores())

# prep. data
d_stan <- list(
  N = nrow(d_long), # total number of trials
  J = max(d_long$id), # number of participants
  K = 3, # number of response options
  y = d_long$rating_abs, # vector of responses
  x = ifelse(d_long$type=='covid',1,0),
  correct = d_long$correct,
  id=d_long$id # subject identifier
)

m01 <- stan(file='./stan/ordinal_2conditions.stan', data = d_stan, iter = 2000,  chains = 4, init_r=0.1)
saveRDS(m01,"ordinal_2conditions.Rds")
m01 <- readRDS("ordinal_2conditions.Rds")
print(m01, pars=c("beta","sigma_u"),probs=c(0.025, 0.975),digits=3)


library(tidybayes)

m01 %>%
  spread_draws(beta[par]) %>%
  filter(par==1) %>%
  mutate(oddsratio = exp(beta)) %>%
  select(oddsratio) %>%
  mean_hdci() -> oddsr_ci


