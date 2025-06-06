# load in the results of all ordinal regression models and make plots for each of them
# (This corresponds to figure 2 and 3 in the main text, and all figures in the SI)

rm(list=ls())
setwd("/mnt/sda2/matteoHDD/git_local_HDD/covid19-misinformation/regression_results")

#
library(rstan)
library(tidyverse)
library(tidybayes)

#
source("../ordinal_helper_functions.R")
source("../ordinal_plot_functions.R")

# ggplot theme
#source("/home/matteo/sync/miscR/miscFunctions.R")
nice_theme <- theme_bw()+theme(text=element_text(family="Helvetica",size=9),panel.border=element_blank(),strip.background = element_rect(fill="white",color="white",size=0),strip.text=element_text(size=rel(0.8)),panel.grid.major.x=element_blank(),panel.grid.major.y=element_blank(),panel.grid.minor=element_blank(),axis.line.x=element_line(size=.4),axis.line.y=element_line(size=.4),axis.text.x=element_text(size=7,color="black"),axis.text.y=element_text(size=7,color="black"),axis.line=element_line(size=.4), axis.ticks=element_line(color="black"))


# loop over each model
for(i in 1:8){

  # plot of fixed-effects parameters
  make_big_plot_c01(i,sd2=T,error_type="data_ci")
  
  # plot of random intercepts
  plot_random_wrapper_c01(i)
}

# ------------------------------------------------------------------------------------ #
# Generate figure 3 i nthe main text

# make marginal effect plots
p_i <- c(1:6,8)
list_index <- seq(1, 2*length(p_i),2)
pl_list <- list()
for(i in 1:length(p_i)){
  pl_list[[list_index[i]]] <- make_single_stackedarea_plot_c01(p_i[i], type=2)
  pl_list[[list_index[i]+1]] <- NULL
}


# -------------------------------------
# distribution of M-ratio values

d <- read_delim("../data/data_allfmt.csv")

xQ <- quantile(d$M_cov, probs=c(0.05, 0.95))
d %>%
  ggplot(aes(x=M_cov)) +
  geom_segment(aes(x=0,xend=xQ[1],y=0.05, yend=0.05), lty=2, col="dark grey",size=0.4)+
  geom_segment(aes(x=0,xend=xQ[2],y=0.95, yend=0.95), lty=2, col="dark grey",size=0.4)+
  geom_segment(aes(x=xQ[1],xend=xQ[1],y=0.05, yend=-1), lty=2, col="dark grey",size=0.4)+
  geom_segment(aes(x=xQ[2],xend=xQ[2],y=0.95, yend=-1), lty=2, col="dark grey",size=0.4)+
  geom_point(data=NULL,aes(x=xQ[1],y=0.05), col="dark grey", pch=21,size=4, stroke=0.4)+
  geom_point(data=NULL,aes(x=xQ[2],y=0.95), col="dark grey", pch=21,size=4, stroke=0.4)+
  stat_ecdf(geom = "step", color="#2f2f75",size=1)+
  #stat_ecdf(aes(x=M_sci),geom = "step", color="dark green")+
  ggtitle("\n")+
  nice_theme +
  coord_cartesian(xlim=c(0.5,1.3),ylim=c(0,1)) +
  labs(x="M-ratio",y="cumulative density") -> Mdist_cumulative


# --------------------------------------------------------------------------- #
## plot effects of meta-efficiency (M-ratio) for all models in a single figure
## expressed as (ordinal) odds-ratio and alongside 95% credible intervals

ordinal_y <- c("att_pubs","att_schools","att_social","behav_meetings","mask_stores","mask_transport","vaccine2")
# "vaccine"

d_all <- {}
for(QN in 1:length(ordinal_y)){

  Q_text <- c(
    "Closing of bars, pubs\nand restaurants",
    "Closing of schools",
    "Restrictions on\nsocial meetings",
    "N. times met others outside\nhousedhold/support bubble",
    "Worn mask at\nsupermarket/grocery store",
    "Worn mask in public transport",
    "How likely to accept vaccine",
    "How likely to accept vaccine"
  )
  sel_Q_text <- Q_text[QN]
  
  ## load fitted model
   m <- readRDS(sprintf("../regression_results/ordinal_%s.RDS",ordinal_y[QN]))
  
  # calculate CI
  # note that I divide by 2, such that the resulting coefficinet represent the expected
  # change in odds when the predictor move by 1 Std. (since it was transformed by dividing by 2 Std.)
  m %>%
    spread_draws(b_X) %>%
    mutate(b_X = -b_X) %>% # invert such that positive == promote health protective behaviours/attitudes
    mutate(b_X = b_X/2) %>%
    mean_hdi() -> qn_i_ci
  
  
  d_i <- data.frame(
    Q = sel_Q_text,
    beta = exp(qn_i_ci$b_X),
    LB = exp(qn_i_ci$.lower),
    UB = exp(qn_i_ci$.upper)
  )
  d_all <- rbind(d_all, d_i)
}

d_all$category <- c(
  rep("Attitudes\ntoward\nrestrictions",3),
  rep("Self-reported\nbehaviours",3),
  "Vaccine\nintent"
)


d_all %>%
  ggplot(aes(y=Q,x=beta, xmin=LB, xmax=UB))+
  geom_point(size=3, color="#2f2f75")+
  geom_vline(xintercept=1,lty=2,color="dark grey")+
  geom_errorbarh(height=.25, color="#2f2f75",size=0.6)+
  labs(x="effect of metacognitive efficiency [odds-ratio]",y="") +
  facet_grid(category~.,scales="free_y",space="free_y")+
  coord_cartesian(xlim=c(1,2.1))+
  scale_x_continuous(breaks=seq(1,2.5, 0.2))+
  nice_theme +
  theme(panel.spacing = unit(0.6, "lines")) -> pl_slopes


# --------------------------------------------------------------------------- #

# unelegant hack to adjust spacing between panels in the figure
pl_list2 <- pl_list
pl_list2[[14]] <- NULL
pl_list2[[15]] <- Mdist_cumulative
pl_list2[[16]] <- NULL
spacing_x <- -0.03

library(ggpubr)

figure01 <- ggarrange(plotlist = pl_list2[1:6], align="hv",ncol = 6,nrow=1, widths = rep(c(1,spacing_x),3), labels = c("A",rep("",6-1)))
figure02 <- ggarrange(plotlist = pl_list2[(1:6)+6], align="hv",ncol = 6,nrow=1, widths = rep(c(1,spacing_x),2), labels = c("C","","D","",""))
figure03 <- ggarrange(plotlist = pl_list2[(1:4)+6*2], align="v",ncol = 6,nrow=1, widths = c(1,spacing_x,1.5,spacing_x,0.5,spacing_x), labels = c("B","","E",""))
figure <- ggarrange(figure01,figure03,figure02,pl_slopes, ncol=2, nrow=2, labels = c("","","","F"),heights=c(1,1.03))

ggsave("ordinal_marginal_effects.pdf",plot=figure,width=8,height=8*2/3)















