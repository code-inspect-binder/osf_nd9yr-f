# -------------------------------------------------------------
# Examine results of meta-d fit and create fig 1
# -------------------------------------------------------------
rm(list=ls())
setwd("/mnt/sda2/matteoHDD/git_local_HDD/covid19-misinformation/jags/")

library(rjags)
library(tidyverse)
library(tidybayes)

# load
output <- readRDS("hdmeta_jags_both_MLok.RDS")

# load helpers
source("wrapper_2conditions.R")

# -------------------------------------------------------------
# load and store individual parameters
d_prep <- readRDS("../data/data_list.RDS")
# d_i <- save_par_i(d_prep, output) # (uncomment this to run for first time)
# saveRDS(d_i, "metaD_i_par.RDS")
d_i <- readRDS("metaD_i_par.RDS")

# merge with dataset for further analyses
d <- read.csv("../data/UoEssex_Results_210414_Client.csv")
d <- merge(d, d_i, by.x="ID", by.y="ID",all=T)
str(d)
write.csv(d, file="../data/data_wSDTpar.csv", row.names=F)


# # recover also type-2 criteria for calculating predictions 
# # (uncomment this to run for first time)
# cS1 <- output[,grep('cS1',varnames(output))]
# cS2 <- output[,grep('cS2',varnames(output))]
# 
# cS1 <- summary(cS1, quantiles=NA) #$statistics[,1]
# cS1 <- cS1$statistics[,1]
# 
# cS2 <- summary(cS2, quantiles=NA) #$statistics[,1]
# cS2 <- cS2$statistics[,1]
# 
# cS12 <- data.frame(cS =unname(c(cS1,cS2)),
#                    resp = c(rep(1,length(cS1)),rep(2,length(cS2))),
#                    label = c(names(cS1),names(cS2)))
# rm(cS1, cS2)
# str(cS12)
# saveRDS(cS12, "type2_cS12.RDS")
cS12 <- readRDS("type2_cS12.RDS")

# tidy criteria parameters
cS12 %>% 
  mutate(ID = str_extract(label,'(?<=\\[)[[:digit:]]+(?=,)'),
         k = str_extract(label,'(?<=,)[[:digit:]](?=,)'),
         type = str_extract(label,'(?<=,)[[:digit:]](?=\\])'),
         type=ifelse(type==1,"science","covid")) %>%
  select(-label) %>%
  pivot_wider(names_from = c(k,resp,type),
              values_from = cS,
              names_prefix = "C",
              id_cols=ID) -> dC

# -------------------------------
# This also would do in theory but in practice run out of memory and crash R
# # reformat using tidybayes
# cS1 %>%
#   spread_draws(cS1[ID,k,type]) %>%
#   group_by(ID) %>%
#   summarise(cS1 = mean(cS1)) 
# -------------------------------

# make sure we have everything
str(d_i)
str(d_prep)
# note that rating in dataset are ordered such that (confident False >>> confident True)
# also S1==c(F1,F2), S2==c(T1,T2)

# # sanity check
# apply(d_prep$T1,2,sum)/(1689 * 7)

# -------------------------------------------------------------
# write function that gives predicted multinomial probabilities given set of parameters
predict_confidence <- function(d_prime, crit, m, b){
  
  nK <- length(b)+1
  
  # Means of SDT distributions
  mu <- m*d_prime
  Tmu <- mu/2
  Fmu <- -mu/2
  
  # p type 1
  pHIT <- 1-pnorm(crit - d_prime/2)
  pFA  <- 1-pnorm(crit + d_prime/2)
  pMIS <- pnorm(crit - d_prime/2)
  pCR  <- pnorm(crit + d_prime/2)
  
  # prob vectors
  theta_F <- rep(NA, nK)
  theta_T <- rep(NA, nK)  
  
  # p | F
  theta_F[1] <- pnorm(b[1] - Fmu)
  for(k in 1:(nK-1)){
    theta_F[k+1] <- pnorm(b[k+1] - Fmu) - pnorm(b[k] - Fmu)
  }
  theta_F[nK] <- 1-pnorm(b[nK-1] - Fmu)
  
  # p | T
  theta_T[1] <- pnorm(b[1] - Tmu)
  for(k in 1:(nK-1)){
    theta_T[k+1] <- pnorm(b[k+1] - Tmu) - pnorm(b[k] - Tmu)
  }
  theta_T[nK] <- 1-pnorm(b[nK-1] - Tmu)
  
  # weighted conditional on type1 probabilities
  k <- nK/2
  theta_T[(k+1):(2*k)] <- theta_T[(k+1):(2*k)]/sum(theta_T[(k+1):(2*k)]) * pHIT
  theta_T[1:k] <- theta_T[1:k]/sum(theta_T[1:k]) * pMIS
  theta_F[(k+1):(2*k)] <- theta_F[(k+1):(2*k)]/sum(theta_F[(k+1):(2*k)]) * pFA
  theta_F[1:k] <- theta_F[1:k]/sum(theta_F[1:k]) * pCR
  
  return(list(pT=theta_T,pF=theta_F))
}


# -------------------------------------------------------------
# now loop and calculate observed counts and predicted probabilities for each ID
d_prep$predT1 <- matrix(NA, nrow=d_prep$J, ncol=d_prep$k*2)
d_prep$predF1 <- matrix(NA, nrow=d_prep$J, ncol=d_prep$k*2)
d_prep$predT2 <- matrix(NA, nrow=d_prep$J, ncol=d_prep$k*2)
d_prep$predF2 <- matrix(NA, nrow=d_prep$J, ncol=d_prep$k*2)

for(id in d_i$ID){
  d_sci <- d_i$d_sci[id]
  c_sci <- d_i$c_sci[id]
  m_sci <- d_i$M_sci[id]
  b_sci <- c(dC$C1_1_science[id], dC$C2_1_science[id], c_sci, dC$C1_2_science[id], dC$C2_2_science[id])
  
  d_cov <- d_i$d_cov[id]
  c_cov <- d_i$c_cov[id]
  m_cov <- d_i$M_cov[id]
  b_cov <- c(dC$C1_1_covid[id], dC$C2_1_covid[id], c_cov, dC$C1_2_covid[id], dC$C2_2_covid[id])
  
  pred_science <- predict_confidence(d_sci, c_sci, m_sci, b_sci)
  pred_covid <- predict_confidence(d_cov, c_cov, m_cov, b_cov)
  
  # # NB: these are the counts
  # d_prep$T1[id,]
  # d_prep$F1[id,]
  # d_prep$T2[id,]
  # d_prep$F2[id,]
  
  d_prep$predT1[id,] <- pred_science$pT
  d_prep$predF1[id,] <- pred_science$pF
  d_prep$predT2[id,] <- pred_covid$pT
  d_prep$predF2[id,] <- pred_covid$pF
  
}


# -------------------------------------------------------------
# FINALLY make all the plots

library(ggplot2)

# nicer white ggplot2 theme
nice_theme <- theme_bw()+theme(text=element_text(family="Helvetica",size=9),panel.border=element_blank(),strip.background = element_rect(fill="white",color="white",size=0),strip.text=element_text(size=rel(0.8)),panel.grid.major.x=element_blank(),panel.grid.major.y=element_blank(),panel.grid.minor=element_blank(),axis.line.x=element_line(size=.4),axis.line.y=element_line(size=.4),axis.text.x=element_text(size=7,color="black"),axis.text.y=element_text(size=7,color="black"),axis.line=element_line(size=.4), axis.ticks=element_line(color="black"))

# pick colors
col_1 <- col2rgb("dark green")
col_2 <- "#2f2f75"
col_1_shade <- rgb(col_1[1],col_1[2],col_1[3],alpha=150, maxColorValue=255)
col_2_shade <- rgb(col2rgb(col_2)[1],col2rgb(col_2)[2],col2rgb(col_2)[3],alpha=150, maxColorValue=255)

# prepare data for plot
pl_d <- data.frame(N = c(apply(d_prep$T1,2,sum),
                             apply(d_prep$T2,2,sum),
                             apply(d_prep$F1,2,sum),
                             apply(d_prep$F2,2,sum)),
                   obs_p = c(apply(d_prep$T1,2,sum)/(1689 * 7),
                             apply(d_prep$T2,2,sum)/(1689 * 7),
                             apply(d_prep$F1,2,sum)/(1689 * 7),
                             apply(d_prep$F2,2,sum)/(1689 * 7)),
                   pred_p = c(apply(d_prep$predT1,2,mean),
                               apply(d_prep$predT2,2,mean),
                               apply(d_prep$predF1,2,mean),
                               apply(d_prep$predF2,2,mean)),
                   pred_p_SD = c(apply(d_prep$predT1,2,sd),
                              apply(d_prep$predT2,2,sd),
                              apply(d_prep$predF1,2,sd),
                              apply(d_prep$predF2,2,sd)),
                   type = c(rep("science",6),
                            rep("covid-19",6),
                            rep("science",6),
                            rep("covid-19",6)),
                   item = c(rep("true",12),rep("false",12)),
                   k = rep(1:6,4),
                   xlabel = c("extremely\nconfident\nFALSE",rep("",4),
                              "extremely\nconfident\nTRUE"))

# multinomial CI (these are simulataneous CI although assume full pooling)
obs_p_CI <- DescTools::MultinomCI(pl_d$N, conf.level=1-0.05/2, method="goodman") *4
pl_d$obs_p_lb <- obs_p_CI[,2]
pl_d$obs_p_ub <- obs_p_CI[,3]

pl_d %>%
  mutate(resp_true = ifelse(k>=4,1,0),
         fill_group = str_c(item, resp_true, type),
         line_group = str_c(item, type),
         item = factor(item, levels=c("true","false"))) %>%
  ggplot(aes(x=k,y=obs_p))+
  facet_grid(item ~ type) +
  nice_theme + 
  geom_col(aes(fill=fill_group))+
  geom_errorbar(aes(y=obs_p,ymin=obs_p_lb, ymax=obs_p_ub),color="dark grey",width=0,lwd=0.8)+ 
  #geom_ribbon(aes(y=pred_p,ymin=pred_p_lb, ymax=pred_p_ub,group=line_group),alpha=0.3,size=NA) +
  geom_line(aes(y=pred_p,group=line_group),size=1.5)+
  #scale_color_manual(values=c("blue","dark green"),guide="none") +
  scale_color_manual(values=c(col_2,"dark green",
                              col_2_shade,col_1_shade,
                              col_2_shade,col_1_shade,
                              col_2,"dark green"),
                     guide="none")+
  scale_fill_manual(values=c(col_2,"dark green",
                             col_2_shade,col_1_shade,
                             col_2_shade,col_1_shade,
                             col_2,"dark green"),
                    guide="none")+
  labs(y="",x="")+
  theme(panel.spacing = unit(1.5, "lines")) +
  scale_x_continuous(breaks=1:6,labels=c("extremely\nconfident\nFALSE",rep("",4),
                                         "extremely\nconfident\nTRUE"))+
  ggtitle(label="Confidence") -> pl01

# -------------------------------------------------------------
# M-ratio posterior density plot
log_mu_sci <- unlist(c(output[,grep('mu_logMratio_1',varnames(output))]))
log_mu_cov <- unlist(c(output[,grep('mu_logMratio_2',varnames(output))]))
d_groupM <- data.frame(logM =c(log_mu_sci,log_mu_cov),
                       type=c(rep('science',length(log_mu_sci)),rep('Covid-19',length(log_mu_sci))))

d_groupM %>%
  mutate(M=exp(logM)) %>%
  ggplot(aes(x=M,y=type,fill=type))+
  nice_theme +
  geom_vline(xintercept = 1, lty=2,size=0.4)+
  stat_halfeye(.width = c(.95, .95),aes(fill=type),alpha=0.8)+
  #scale_fill_manual(values=c(col_2_shade,col_1_shade),name="",guide='none') +
  scale_fill_manual(values=c(col_2,"dark green"),name="",guide="none") +
  #theme(legend.position=c(.9,.75))+
  scale_x_continuous("M-ratio\n \n ",limits=c(0.5,1.14))+
  labs(y="")+
  theme(axis.text.y = element_text(size=10))+
  ggtitle(label="Metacognitive efficiency")-> pl02


# -------------------------------------------------------------
# accuracy plot
d <- read.csv("../data/UoEssex_Results_210414_Client.csv")

d_sub <- d[,21:(21+27)]
colnames(d)[21:(21+27)]
colnames(d)[21:(21+13)]
colnames(d)[35:(35+13)]
col_science <- 21:(21+13)
col_covid <- 21:(21+13)

d_sub$id <- 1:nrow(d_sub)
d_long <- as.data.frame(pivot_longer(d_sub, cols=-id, names_to="question",values_to = "rating_raw"))

# get function to check accuracy
Q <- list(labels = colnames(d)[21:(21+27)],
          type = c(rep("science",14),rep("covid",14)),
          true = c(c(0,0,0,0,0,0, 1,1,1,1,1,1,1,0), c(0,0,0,0,0,0,1,1,1,1,0,1,1,1)))

# 
d_long$respT <- ifelse(d_long$rating_raw<=3,1,0)
d_long$type <- ifelse(as.numeric(substr(d_long$question,7,8))<=14,"science","covid")
d_long$isTrue <- Q$true[as.numeric(substr(d_long$question,7,8))]
d_long$correct <- ifelse(d_long$respT==d_long$isTrue,1,0)
d_long$confidence <- ifelse(d_long$rating_raw<=3,d_long$rating_raw,d_long$rating_raw-3)
d_long$confidence[d_long$rating_raw<=3] <- ifelse(d_long$confidence[d_long$rating_raw<=3]==3,1,ifelse(d_long$confidence[d_long$rating_raw<=3]==1,3,2))

d_long %>%
  group_by(id,type) %>%
  summarise(mean_correct = mean(correct),
            ci_lb=NA,ci_ub=NA) %>%
  mutate(type=ifelse(type=="covid","covid-19",type)) -> dag

d_long %>%
  group_by(id,type) %>%
  summarise(correct = mean(correct)) %>%
  group_by(type) %>%
  summarise(mean_correct = mean(correct),
            correct_se = mlisi::bootMeanSE(correct),
            ci_lb = mlisi::bootMeanCI(correct,nsim = 1e3,)[1],
            ci_ub = mlisi::bootMeanCI(correct,nsim = 1e3,)[2]) %>%
  mutate(type=c("covid-19","science")) -> dag2


dag %>%
  ggplot(aes(x=type,y=mean_correct)) +
  nice_theme + 
  scale_fill_manual(values=c(col_2,"dark green"),guide="none") +
  scale_color_manual(values=c(col_2,"dark green"),guide="none") +
  #geom_violin(aes(fill=type),color="white",bw=0.04,alpha=0.5) +
  geom_boxplot(aes(fill=type,color=type),notch=F,alpha=0.5,outlier.shape="*",size=0.8) +
  #coord_cartesian(ylim=c(0,1))+
  #geom_point(data=dag2,aes(color=type,fill=type),size=4,pch=21,color="white")+
  geom_point(data=dag2,aes(color=type,fill=type),size=3,pch=21,color="black")+
  #geom_errorbar(data=dag2,aes(ymin=ci_lb, ymax=ci_ub),width=0,size=1)+
  labs(x="",y="fraction correct")+
  theme(axis.text.x=element_text(angle = 55, hjust = 1))+
  scale_y_continuous(breaks=seq(0,1,0.1))+
  ggtitle(label="Accuracy",subtitle="") -> pl03
  

# -------------------------------------------------------------
# ROC curves
library(yardstick)

d_long %>%
  group_by(type) %>%
  mutate(correct = factor(ifelse(correct==1,"correct","error"))) %>%
  roc_curve(estiamted=confidence, truth = correct) %>%
  ggplot(aes(x = 1 - specificity, 
             y = sensitivity, 
             color = type)) + 
  geom_abline(slope = 1, intercept = 0, size = 0.6,lty=2)+
  geom_line(size = 2) +
  #geom_point(size=3.5,pch=21)  +
  scale_color_manual(values = c(col_2,"dark green"),name="",guide='none') +
  coord_fixed()+
  nice_theme +
  labs(x="type-2  False Alarms", y="type-2  Hits") +
  theme(legend.position=c(.9,.25))+
  ggtitle(label="type-2 ROC")  -> pl_ROC


d_long %>%
  group_by(type) %>%
  mutate(isTrue = factor(ifelse(isTrue==1,"true","false")),
         type = ifelse(type=="covid", "covid-19",type)) %>%
  roc_curve(estiamted=rating_raw, truth = isTrue) %>%
  ggplot(aes(x = 1 - specificity, 
             y = sensitivity, 
             color = type)) + 
  geom_abline(slope = 1, intercept = 0, size = 0.6,lty=2)+
  geom_line(size = 2) +
  #geom_point(size=3.5,pch=21)  +
  scale_color_manual(values = c(col_2,"dark green"),name="") +
  coord_fixed()+
  nice_theme +
  labs(x="type-1  False Alarms", y="type-1  Hits") +
  theme(legend.position=c(.8,.25))+
  ggtitle(label="type-1 ROC")  -> pl_ROC_type1


# -------------------------------------------------------------
# arrange in one single figure

library(ggpubr)

figure <- ggarrange(ggarrange(pl_ROC_type1, pl_ROC, nrow = 2, heights = c(1,1), align="h", labels = c("A", "B")),
                    NULL,
                    ggarrange(
                      ggarrange(pl03, pl01, nrow = 1, ncol=2, widths =  c(0.6,1.4), labels = c("C", "D")),
                      pl02, heights = c(1.3,0.7), nrow=2, ncol=1, labels=c("","E")),
                    nrow=1, ncol=3, align="v", widths =  c(0.6,0.08,1.4))

ggsave("../Fig1.pdf",plot=figure,width=8,height=6)


