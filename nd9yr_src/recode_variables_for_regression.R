# -------------------------------------------------------------
# ordinal regression analyses
# prepare dataset and recode variables

rm(list=ls())
setwd("/mnt/sda2/matteoHDD/git_local_HDD/covid19-misinformation/")

library(tidyverse)
library(rstan)
options(mc.cores = parallel::detectCores())

#
d <- read.csv("./data/data_wSDTpar.csv")
#codes <- read.csv2("./data/code_book.csv")
codebook <- read.csv("./data/codebook_stan_model.csv")

# -------------------------------------------------------------
## prepare data for modelling
## see codebook for numerical codes

##### -------------------------------------------------------------
##### ordinal dependent variables

## ATTITUDES TOWARD RESTRICTIONS
d$att_pubs <- ifelse(d$UoE_2_1<6, d$UoE_2_1, NA)
d$att_schools <- ifelse(d$UoE_2_2<6, d$UoE_2_2, NA)
d$att_social <- ifelse(d$UoE_2_3<6, d$UoE_2_3, NA)

## MEETING SOCIALLY WITH OTHERS (exclude don't know responses)
# for "UoE_3", 5 == 'don't know'
d$behav_meetings <- ifelse(d$UoE_3<5, d$UoE_3, NA)

## MASK WEARING
# masks - remove 'not applicable' (e.g. health reasons) responses
d$mask_stores <- ifelse(d$UoE_4_1<6, d$UoE_4_1, NA)
d$mask_transport <- ifelse(d$UoE_4_2<6, d$UoE_4_2, NA)

## COVID-19 VACCINES
# all response options (including people who were already vaccinated)
d$vaccine <- ifelse(d$UoE_6==6, NA, d$UoE_6) # remove answers 'medically exempt' (only 4)
d$vaccine <- ifelse(d$vaccine==7, 0, d$vaccine) # consider a certain 'yes' if they had already 1 dose
d$vaccine <- ifelse(is.element(d$vaccine,c(3,4,5)), d$vaccine + 1, d$vaccine) # increment by 1
d$vaccine <- ifelse(d$vaccine==6, 3, d$vaccine)
d$vaccine <- d$vaccine  + 1

# vaccines - now excluding people who had it already 
# (This is the analysis that is presented into the paper, as I think it makes more sense
# however even including people who already had vaccines - as done above - lead to similar results)
d$vaccine2 <- ifelse(d$UoE_6==6, NA, d$UoE_6) # remove answers 'medically exempt' (only 4)
d$vaccine2 <- ifelse(d$vaccine2==7, NA, d$vaccine2) # consider a certain 'yes' if they had already 1 dose
d$vaccine2 <- ifelse(is.element(d$vaccine2,c(3,4,5)), d$vaccine2 + 1, d$vaccine2) # increment by 1
d$vaccine2 <- ifelse(d$vaccine2==6, 3, d$vaccine2)
tapply(d$vaccine2,d$UoE_6,mean)
tapply(d$vaccine2,d$vaccine2,length)


##### -------------------------------------------------------------
##### predictors

## VOTE & POLITICAL ALIGNMENTS
d %>%
  group_by(pastvote_ge_2019,voted_ge_2019) %>%
  summarise(N=n()) -> v2019summary

v2019summary$N[v2019summary$voted_ge_2019==99]/sum(v2019summary$N) * 100


d$voteEURef <- d$pastvote_EURef #legacy
with(d, tapply(voteEURef, list(voteEURef), length))
31/1689 * 100 # respondents who can't remember who they voted in EU ref.

d %>%
  mutate(vote2019 = case_when(
    voted_ge_2019==1 & pastvote_ge_2019 < 90 ~ pastvote_ge_2019,
    voted_ge_2019==1 & pastvote_ge_2019 == 98 ~ 8L,
    voted_ge_2019==1 & pastvote_ge_2019 == 99 ~ 9L,
    voted_ge_2019==2  ~ 10L,
    voted_ge_2019==99 ~ 11L
  )) -> d

# vote in 2019 
with(d, tapply(voted_ge_2019 ,list(voted_ge_2019), length)) # sanity checks
with(d, tapply(vote2019, list(vote2019), length))
sum(with(d, tapply(vote2019, list(vote2019), length))) == nrow(d)

## AFFECTED BY COVID
d$covid_affected <- ifelse(d$UoE_5==1,1,0)


## INCOME
# http://www.nrs.co.uk/nrs-print/lifestyle-and-classification-data/social-grade/
# https://www.mrs.org.uk/resources/social-grade


# income (personal)
d$profile_gross_personal[which(is.na(d$profile_gross_personal))] <- max(d$profile_gross_personal,na.rm=T) +1
tapply(d$profile_gross_personal, d$profile_gross_personal, length)
d$income <- NA
for(i in 1:nrow(d)){
  d$income[i] <- codebook$code_alternative[codebook$question=="profile_gross_personal" & 
                                                codebook$code_n==d$profile_gross_personal[i] & 
                                                !is.na(codebook$code_n)]
}
round(tapply(d$income, d$income, length)/sum(tapply(d$income, d$income, length)) * 100, digits=2)


# marital status
d$profile_marital_stat[which(is.na(d$profile_marital_stat))] <- 8 # unknown / skipped / not asked
tapply(d$profile_marital_stat, d$profile_marital_stat, length)
79/1689 * 100 

# education
tapply(d$profile_education_level, d$profile_education_level, length)
d$education <- NA
for(i in 1:nrow(d)){
  d$education[i] <- codebook$code_alternative[codebook$question=="profile_education_level" & 
                                                codebook$code_n==d$profile_education_level[i] & 
                                                !is.na(codebook$code_n)]
}


### SAVE
write.table(d,file="./data/data_allfmt.csv", row.names=F,sep="\t",quote=F)


