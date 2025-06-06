# ------------------------------------------------------------
rm(list=ls())
setwd("/mnt/sda2/matteoHDD/git_local_HDD/covid19-misinformation")


# ------------------------------------------------------------
library(tidyverse)
library(ggplot2)

# ggplot theme
source("/mnt/D/Dropbox/sync/miscR/miscFunctions.R")


# ------------------------------------------------------------
# load long data
dLO <- read.table("./data/data_long.csv", header=T,sep=";")

# questions text
qtext <- read_delim("./data/code_book.csv")

# add condition & response labels
dLO$rating_label_full <- factor(dLO$rating_raw, levels=6:1, labels=c("extremely confident\nTRUE","fairly confident\nTRUE","not at all confident\nTRUE","not at all confident\nFALSE","fairly confident\nFALSE","extremely confident\nFALSE")[6:1],ordered=T)

dLO$rating_label <- factor(dLO$rating_raw, levels=6:1, labels=c("extremely\nconfident\nTRUE"," ","  ","   ","    ","extremely\nconfident\nFALSE")[6:1],ordered=T)

dLO$type_label <- ifelse(dLO$type=="covid","COVID-19","General science")
dLO$respT_lab <- ifelse(dLO$respT==1,"true","false")


# ------------------------------------------------------------
# calculate raw and weighted percentages

# sanity check
dLO %>%
  group_by(type_label, question) %>%
  summarise(N=length(rating_raw),
            Nw =sum(WEIGHT)) %>%
  pull(N) %>% unique() -> N_tot

# prepare data
dLO %>%
  group_by(type_label, respT_lab, question, rating_raw, rating_label, isTrue) %>%
  summarise(N=length(rating_raw),
            N_w =sum(WEIGHT)) %>%
  mutate(percentage = 100*N/N_tot, 
         percentage_w = 100*N_w/N_tot) -> d_plot

d_plot$Q_text <- NA
for(q_i in unique(d_plot$question)){
  d_plot$Q_text[d_plot$question==q_i] <- qtext$Label[qtext$Variable==q_i]
}

# make plots
d_plot %>%
  filter(type_label!="COVID-19", isTrue==1) %>%
  ggplot(aes(x=rating_label,fill=respT_lab))+
  geom_col(aes(y=percentage),width=0.95) +
  geom_col(aes(y=percentage_w),fill=NA,lty=2,size=0.4,color="black",width=0.95) +
  scale_fill_viridis_d(option="A",begin=0.65,end=0.2,name=NULL,guide='none')+
  facet_grid(~Q_text,  labeller = labeller(Q_text = label_wrap_gen(25)))+
  labs(y="%",x="")+
  nice_theme +
  theme(panel.spacing.x = unit(2, "lines")) +
  ggtitle("General science, True statements") -> sciT

d_plot %>%
  filter(type_label!="COVID-19", isTrue==0) %>%
  ggplot(aes(x=rating_label,fill=respT_lab))+
  geom_col(aes(y=percentage),width=0.95) +
  geom_col(aes(y=percentage_w),fill=NA,lty=2,size=0.4,color="black",width=0.95) +
  scale_fill_viridis_d(option="A",begin=0.65,end=0.2,name=NULL,guide='none')+
  facet_grid(~Q_text,  labeller = labeller(Q_text = label_wrap_gen(25)))+
  labs(y="%",x="")+
  nice_theme +
  theme(panel.spacing.x = unit(2, "lines")) +
  ggtitle("General science, False statements") -> sciF

d_plot %>%
  filter(type_label=="COVID-19", isTrue==1) %>%
  ggplot(aes(x=rating_label,fill=respT_lab))+
  geom_col(aes(y=percentage),width=0.95) +
  geom_col(aes(y=percentage_w),fill=NA,lty=2,size=0.4,color="black",width=0.95) +
  scale_fill_viridis_d(option="A",begin=0.65,end=0.2,name=NULL,guide='none')+
  facet_grid(~Q_text,  labeller = labeller(Q_text = label_wrap_gen(25)))+
  labs(y="%",x="")+
  nice_theme +
  theme(panel.spacing.x = unit(2, "lines")) +
  ggtitle("COVID-19, True statements") -> covT

d_plot %>%
  filter(type_label=="COVID-19", isTrue==0) %>%
  ggplot(aes(x=rating_label,fill=respT_lab))+
  geom_col(aes(y=percentage),width=0.95) +
  geom_col(aes(y=percentage_w),fill=NA,lty=2,size=0.4,color="black",width=0.95) +
  scale_fill_viridis_d(option="A",begin=0.65,end=0.2,name=NULL,guide='none')+
  facet_grid(~Q_text,  labeller = labeller(Q_text = label_wrap_gen(25)))+
  labs(y="%",x="")+
  nice_theme +
  theme(panel.spacing.x = unit(2, "lines")) +
  ggtitle("COVID-19, False statements") -> covF
  

# combine plots
library(ggpubr)
figure <- ggarrange(sciT,sciF,covT,covF, nrow = 4, labels = c("A", "B","C","D"), align="h")
figure

ggsave("responses_by_Q.pdf",plot=figure,width=11,height=11)



