# -------------------------------------------------------------
# load dataset and put it in format suitable for modelling
# Matteo Lisi 2020

# -------------------------------------------------------------
rm(list=ls())
setwd("~/git_local/covid19-misinformation")

# -------------------------------------------------------------
# load
d <- read.csv("./data/UoEssex_Results_210414_Client.csv")
str(d)

# -------------------------------------------------------------
library(tidyr)
library(dplyr)
library(ggplot2)

# -------------------------------------------------------------
d_sub <- d[,21:(21+27)] # subset columns with answers to science/covid statements
colnames(d)[21:(21+27)] # sanity check: all questions
colnames(d)[21:(21+13)] # science items
colnames(d)[35:(35+13)] # covid items

d_sub$id <- d$ID # participant identifier
d_long <- as.data.frame(pivot_longer(d_sub, cols=-id, names_to="question",values_to = "rating_raw"))

# categorize question labels
Q <- data.frame(labels = colnames(d)[21:(21+27)],
                type = c(rep("science",14),rep("covid",14)),
                true = c(c(0,0,0,0,0,0, 1,1,1,1,1,1,1,0), c(0,0,0,0,0,0,1,1,1,1,0,1,1,1)))

# sanity check with codebook
d_codes <- read.table("./data/code_book.csv",header=T,sep=";")

# recode responses in long database
d_long$respT <- ifelse(d_long$rating_raw<=3,1,0)
d_long$type <- ifelse(as.numeric(substr(d_long$question,7,8))<=14,"science","covid")
d_long$isTrue <- Q$true[as.numeric(substr(d_long$question,7,8))]
d_long$correct <- ifelse(d_long$respT==d_long$isTrue,1,0)

# reverse code confidence rating (so that goes confident False >>> confident True)
d_long$rating_rev <- 7 - d_long$rating_raw

# sanity check plot 1
ggplot(d_long,aes(x=rating_rev,color=respT,fill=respT,group=respT))+geom_bar()+facet_grid(isTrue~type)

# include also sampling weight to long dataset
d_long$WEIGHT <- NA
for(i in unique(d_long$id)){
  d_long$WEIGHT[d_long$id==i] <- d$WEIGHT[d$ID==i]
}

write.table(d_long,"./data/data_YouGov_April_long.csv",col.names=T,row.names=F,sep=";")

# -------------------------------------------------------------
# prepare data for modelling meta-d
T1 <- with(d_long[d_long$type=="science" & d_long$isTrue==1,], tapply(rating_rev, list(id, rating_rev), length))
F1 <- with(d_long[d_long$type=="science" & d_long$isTrue==0,], tapply(rating_rev, list(id, rating_rev), length))
T2 <- with(d_long[d_long$type=="covid" & d_long$isTrue==1,], tapply(rating_rev, list(id, rating_rev), length))
F2 <- with(d_long[d_long$type=="covid" & d_long$isTrue==0,], tapply(rating_rev, list(id, rating_rev), length))

T1[is.na(T1)] <- 0
F1[is.na(F1)] <- 0
T2[is.na(T2)] <- 0
F2[is.na(F2)] <- 0

# n sj & rating levels
J <- dim(T1)[1]
k <- dim(T1)[2]/2

d_prep <- list(
  J, k, d$WEIGHT,
  T1, F1, T2, F2
)
names(d_prep) <- c("J","k","weight","T1","F1","T2","F2")
str(d_prep)

saveRDS(d_prep, "./data/data_list.RDS")



