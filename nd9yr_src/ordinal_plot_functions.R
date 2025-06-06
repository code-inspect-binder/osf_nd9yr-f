# ---------------------------------------------------------------------------- 
# plotting functions for model with Type-1 d' as predictor

make_big_plot_c01 <- function(QN,sd2=T,  error_type="hdi89"){
  
  ## load data & codes
  d <- read_delim("./data/data_allfmt.csv",show_col_types = FALSE)
  codes <- read_delim("./data/codebook_stan_model.csv",show_col_types = FALSE)
  codelabels <- read_delim("./data/code_book.csv",show_col_types = FALSE)
  
  ## select dependent variable and analysis
  ordinal_y <- c("att_pubs","att_schools","att_social","behav_meetings","mask_stores","mask_transport","vaccine","vaccine2")
  var_y <- ordinal_y[QN]
  
  ## get question label
  Q_id <- unique(codes$question[codes$code_question==ordinal_y[QN]])
  Q_id <- Q_id[which(!is.na(Q_id))]
  
  
  Q_text <- c(
    "Closing of bars, pubs\nand restaurants",
    "Closing of schools",
    "Restrictions on social meetings",
    "N. times met others outside\nhousedhold/support bubble",
    "Worn mask at supermarket/grocery store",
    "Worn mask in public transport",
    "How likely to accept vaccine",
    "How likely to accept vaccine"
  )
  sel_Q_text <- Q_text[QN]
  
  ## losd fitted model
  m <-  readRDS(sprintf("../regression_results/ordinal_%s.RDS",ordinal_y[QN]))
  
  ## load prediction table (if present), otherwise create it
  if(any(dir()==str_c(ordinal_y[QN],"_Ptab.RDS"))){
      P <- readRDS(str_c(ordinal_y[QN],"_Ptab.RDS"))
    }else{
      P <- predicted_p_full_extended_c(m, d, sd2=sd2)
      saveRDS(P, file=str_c(ordinal_y[QN],"_Ptab.RDS"))
    }

  
  ## answer labels
  codes %>%
    filter(code_question==!!var_y) %>%
    arrange(code_n) %>%
    pull(Label) %>%
    str_replace_all("\x92","'") -> x_lab
  
  # fix vaccine answer label for plotting
  # (so that it fits in plot)
  if(QN==7 & var_y=="vaccine"){
    x_lab[1] <- "Already received\nat least 1 dose"
    x_lab[4] <- "Don't know"
  }
  
  ## assign a group variable to responses, to facilitate plotting
  d %>%
    mutate(group = .data[[!!var_y]]) -> d
  
  ## remove NA response (e.g. if people didn't take public transport at all)
  obs_index <- !is.na(d$group)
  P <- P[,obs_index,]
  d %>%
    filter(!is.na(group)) -> d
  
  ## calculate mean predictions & CI
  P_av <- apply(P, 3, mean)
  if(error_type=="hdi"){
    P_ci <- apply(P, 3, hdi)
  }else if(error_type=="se"){
    P_ci <- rbind(P_av - apply(P, 3, sd),
                  P_av + apply(P, 3, sd))
    P_ci[P_ci<0] <- 0
    P_ci[P_ci>1] <- 1
  }else if(error_type=="hdi89"){
    P_ci <- apply(P, 3, hdi, prob=0.89)
  }
  
  
  ## N response options
  k <- length(P_av)
  
  ## plot overall predictions
  library(ggtext)
  
  if(error_type!="data_ci"){
  d %>%
    #mutate(group = .data[[!!var_y]]) %>%
    group_by(group) %>%
    select(!!var_y) %>%
    summarise(N = length(.data[[!!var_y]]),
              P_obs = N/nrow(d)) %>%
    mutate(P_pred = P_av,
           P_lb = P_ci[1,],
           P_ub = P_ci[2,]) %>%
    ggplot(aes(x=group))+
    geom_col(aes(y=P_obs))+
    geom_ribbon(aes(ymin=P_lb, ymax=P_ub),alpha=0.3,color="black",lty=3,size=NA) +
    geom_line(aes(y=P_pred),size=1.5)+
    scale_x_continuous(breaks = 1:k, labels=x_lab)+
    nice_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title=element_text(face="bold"),
          plot.subtitle=element_text(face="italic")) + 
    labs(x='',y='proportion of responses') +
    ggtitle(label=sel_Q_text,subtitle="\n\n") -> pl01
  
  }else{
    
    d %>%
      #mutate(group = .data[[!!var_y]]) %>%
      group_by(group) %>%
      select(!!var_y) %>%
      summarise(N = length(.data[[!!var_y]]),
                P_obs = N/nrow(d)) %>%
      mutate(P_pred = P_av,
             LB = DescTools::MultinomCI(N, conf.level=1-0.05/2, method="goodman")[,2],
             UB = DescTools::MultinomCI(N, conf.level=1-0.05/2, method="goodman")[,3]) %>%
      ggplot(aes(x=group))+
      geom_col(aes(y=P_obs))+
      geom_errorbar(aes(ymin=LB, ymax=UB),color="dark grey",width=0,lwd=1) +
      geom_line(aes(y=P_pred),size=1.5)+
      scale_x_continuous(breaks = 1:k, labels=x_lab)+
      nice_theme +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title=element_text(face="bold"),
            plot.subtitle=element_text(face="italic")) + 
      labs(x='',y='proportion of responses') +
      ggtitle(label=sel_Q_text,subtitle="\n\n") -> pl01
  }
  
  # plot coefficients
  m %>%
    spread_draws(b_X,b_d, b_age, b_age2, b_gender, b_covid_affected) %>%
    pivot_longer(cols=starts_with("b_"),names_to='par',values_to='beta') %>%
    mutate(pcol=ifelse(par=="b_X","a","b"),
           par_label=case_when(
             par=="b_age" ~ "age",
             par=="b_d" ~ "d-prime",
             par=="b_age2" ~ "age^2",
             par=="b_gender" ~ "gender",
             par=="b_covid_affected" ~ "affected\nby Covid-19",
             par=="b_X" ~ "metacognitive\nefficiency"
           )) %>%
    mutate(par_label = factor(par_label,levels=unique(par_label))) %>%
    ggplot(aes(y = par_label, x = beta)) +
    geom_vline(xintercept = 0, lty=2) +
    stat_halfeye(.width = c(.95, .5),aes(fill=pcol)) +
    nice_theme+
    #scale_fill_manual(values=c(rgb(0,0,1,0.8),rgb(0,0,0,0.5)),guide="none")+
    scale_fill_manual(values=c("#2f2f75",rgb(0,0,0,0.5)),guide="none")+
    labs(x=expression(beta),y="") +
    ggtitle(label=" \n ",subtitle="Parameters\n(slopes)") -> pl02
  
  # plot random effect variances
  m %>%
    spread_draws(sigma_grade,sigma_vote2019,sigma_voteEUref,sigma_edu,sigma_region,sigma_marital_stat,sigma_income) %>%
    select(starts_with("sigma_")) %>%
    pivot_longer(cols=starts_with("sigma_"),names_to='par',values_to='sigma') %>%
    mutate(par_label=case_when(
      par=="sigma_grade" ~ "CIE social grade",
      par=="sigma_vote2019" ~ "Vote 2019 election",
      par=="sigma_voteEUref" ~ "Vote EU referendum",
      par=="sigma_income" ~ "Household income",
      par=="sigma_marital_stat" ~ "Marital status",
      par=="sigma_edu" ~ "Education",
      par=="sigma_region" ~ "Geographical\nregion"
    )) %>%
    ggplot(aes(y = par_label, x = sigma)) +
    stat_halfeye(.width = c(.95, .5)) +
    coord_cartesian(xlim=c(0,1.5))+
    nice_theme+
    scale_fill_manual(values=c(rgb(1,0,0,0.8),rgb(0,0,0,0.5)),guide="none")+
    labs(x=expression(sigma),y="") +
    ggtitle(label=" \n ",subtitle="Parameters\n(betweem groups SD)") -> pl02.2
  
  library(ggpubr)
  
  figure <- ggarrange(pl01, pl02, pl02.2, ncol = 3, labels = c("A", "B","C"))
  #figure
  ggsave(str_c(ordinal_y[QN],"_c0",type,".pdf"),plot=figure,width=9.8/1.2,height=9/2/1.2)
  return(figure)
}

# ------------------------------------------------------------------------------------ #
# wrapper to plot random intercepts - extended model with all predictors
# (include also income and marital status)

plot_random_wrapper_c01 <- function(QN, type=1){
  
  x_lim_range <- c(-4,4)
  prop_scal_limit <- c(0,0.5)
  
  d <- read_delim("../data/data_allfmt.csv",show_col_types = FALSE)
  
  ## load codes
  codes <- read_delim("../data/codebook_stan_model.csv",show_col_types = FALSE)
  codelabels <- read_delim("../data/code_book.csv",show_col_types = FALSE)
  
  ## select dependent variable and analysis
  ordinal_y <- c("att_pubs","att_schools","att_social","behav_meetings","mask_stores","mask_transport","vaccine","vaccine2")
  var_y <- ordinal_y[QN]
  
  ## get question label
  Q_id <- unique(codes$question[codes$code_question==ordinal_y[QN]])
  Q_id <- Q_id[which(!is.na(Q_id))]
  
  Q_text <- c(
    "Closing of bars, pubs and restaurants",
    "Closing of schools",
    "Restrictions on social meetings",
    "N. times met others outside housedhold/support bubble",
    "Worn mask at supermarket/grocery store",
    "Worn mask in public transport",
    "How likely to accept vaccine",
    "How likely to accept vaccine"
  )
  sel_Q_text <- Q_text[QN]
  
  ## losd fitted model
  m <-  readRDS(sprintf("../regression_results/ordinal_%s.RDS",ordinal_y[QN]))
  
  nice_theme <- theme_bw()+theme(text=element_text(family="Helvetica",size=9),panel.border=element_blank(),strip.background = element_rect(fill="white",color="white",size=0),strip.text=element_text(size=rel(0.8)),panel.grid.major.x=element_blank(),panel.grid.major.y=element_blank(),panel.grid.minor=element_blank(),axis.line.x=element_line(size=.4),axis.line.y=element_line(size=.4),axis.text.x=element_text(size=7,color="black"),axis.text.y=element_text(size=7,color="black"),axis.line=element_line(size=.4), axis.ticks=element_line(color="black"))
  
  # ------------------
  # CIE
  # 6 lines
  
  codes %>%
    filter(question=="profile_socialgrade_cie",
           !is.na(code_n)) %>%
    select(code_n, Label)  %>%
    mutate(Label = str_replace_all(Label,"\x92","'"))->  code_tab
  
  d %>%
    group_by(profile_socialgrade_cie) %>%
    summarise(N = n()) %>%
    ungroup() %>%
    mutate(prop = N/sum(N)) %>%
    rename(code_n = profile_socialgrade_cie) %>%
    right_join(.,code_tab, by="code_n") -> code_tab
  
  m %>%
    spread_draws(b_grade[code_n]) %>%
    left_join(., code_tab, by="code_n") %>%
    mutate(Label = factor(Label)) %>%
    ggplot(aes(y = Label, x = b_grade, fill=prop)) +
    stat_halfeye(.width = c(.95, .5),alpha=1)+
    scale_fill_viridis_c(limits=prop_scal_limit, name="fraction of\nresponses", guide='none') +
    geom_vline(xintercept=0,lty=2) +
    nice_theme +
    coord_cartesian(xlim=x_lim_range)+
    labs(x=expression(mu),y="") +
    ggtitle("CIE social grade") -> pl_grade
  
  
  
  # -------------------------
  # 2019GE
  # 11 lines
  
  codes %>%
    filter(question=="vote2019",
           !is.na(code_n)) %>%
    select(code_n, Label)  %>%
    mutate(Label = str_replace_all(Label,"\x92","'"))->  code_tab
  
  d %>%
    group_by(vote2019) %>%
    summarise(N = n()) %>%
    ungroup() %>%
    mutate(prop = N/sum(N,na.rm=T)) %>%
    rename(code_n = vote2019) %>%
    right_join(.,code_tab, by="code_n") -> code_tab
  
  m %>%
    spread_draws(b_vote2019[code_n]) %>%
    left_join(., code_tab, by="code_n") %>%
    ggplot(aes(y = fct_reorder(Label, b_vote2019, .fun=mean), x = b_vote2019, fill=prop)) +
    stat_halfeye(.width = c(.95, .5),alpha=1)+
    scale_fill_viridis_c(limits=prop_scal_limit, name="fraction of\nresponses", guide='none') +
    geom_vline(xintercept=0,lty=2) +
    nice_theme +
    coord_cartesian(xlim=x_lim_range)+
    labs(x=expression(mu),y="") +
    ggtitle("Vote in 2019\nGeneral Election") -> pl_v2019
  
  
  
  # ------------------
  # edu
  # 20 lines
  
  if(type==1){
  codes %>%
    filter(question=="profile_education_level",
           !is.na(code_n)) %>%
    select(code_n, Label)  %>%
    mutate(Label = str_replace_all(Label,"\x92","'"))->  code_tab
  
  d %>%
    group_by(profile_education_level) %>%
    summarise(N = n()) %>%
    ungroup() %>%
    mutate(prop = N/sum(N,na.rm=T)) %>%
    rename(code_n = profile_education_level) %>%
    right_join(.,code_tab, by="code_n") -> code_tab
  
  m %>%
    spread_draws(b_edu[code_n]) %>%
    left_join(., code_tab, by="code_n") %>%
    mutate(Label = factor(Label)) %>%
    ggplot(aes(y = Label, x = b_edu, fill=prop)) +
    stat_halfeye(.width = c(.95, .5),alpha=1)+
    scale_fill_viridis_c(limits=prop_scal_limit, name="fraction of\nresponses") +
    geom_vline(xintercept=0,lty=2) +
    nice_theme +
    coord_cartesian(xlim=x_lim_range)+
    labs(x=expression(mu),y="") +
    ggtitle("Education") -> pl_edu
  
  }else{
    codes %>%
      filter(question=="profile_education_level",
             !is.na(code_n)) %>%
      select(code_alternative, alternative_label)  %>%
      group_by(alternative_label) %>%
      summarise(code_alternative=mean(code_alternative)) %>%
      mutate(Label = str_replace_all(alternative_label,"\x92","'"),
             code_n = code_alternative)->  code_tab
    
    code_tab$Label <- fct_reorder(code_tab$Label, 1:7)
    
    d %>%
      group_by(education) %>%
      summarise(N = n()) %>%
      ungroup() %>%
      mutate(prop = N/sum(N,na.rm=T)) %>%
      rename(code_n = education) %>%
      right_join(.,code_tab, by="code_n") -> code_tab
    
    m %>%
      spread_draws(b_edu[code_n]) %>%
      left_join(., code_tab, by="code_n") %>%
      ggplot(aes(y = Label, x = b_edu, fill=prop)) +
      stat_halfeye(.width = c(.95, .5),alpha=1)+
      scale_fill_viridis_c(limits=prop_scal_limit, name="fraction of\nresponses") +
      geom_vline(xintercept=0,lty=2) +
      nice_theme +
      coord_cartesian(xlim=x_lim_range)+
      labs(x=expression(mu),y="") +
      ggtitle("Education") -> pl_edu
  }
  
  # ------------------
  # EU ref
  # 4 lines
  
  codes %>%
    filter(question=="pastvote_EURef",
           !is.na(code_n)) %>%
    select(code_n, Label) %>%
    mutate(Label = str_replace_all(Label,"\x92","'"))-> code_tab
  
  d %>%
    group_by(pastvote_EURef) %>%
    summarise(N = n()) %>%
    ungroup() %>%
    mutate(prop = N/sum(N,na.rm=T)) %>%
    rename(code_n = pastvote_EURef) %>%
    right_join(.,code_tab, by="code_n") -> code_tab
  
  m %>%
    spread_draws(b_voteEUref[code_n]) %>%
    left_join(., code_tab, by="code_n") %>%
    ggplot(aes(y = fct_reorder(Label, b_voteEUref, .fun=mean), x = b_voteEUref, fill=prop)) +
    stat_halfeye(.width = c(.95, .5),alpha=1)+
    scale_fill_viridis_c(limits=prop_scal_limit, name="fraction of\nresponses", guide='none') +
    geom_vline(xintercept=0,lty=2) +
    nice_theme +
    coord_cartesian(xlim=x_lim_range)+
    labs(x=expression(mu),y="") +
    ggtitle("Vote in EU referendum") -> pl_vEU
  
  # ------------------
  # Region
  # 11 lines
  
  codes %>%
    filter(question=="profile_GOR",
           !is.na(code_n)) %>%
    select(code_n, Label)  %>%
    mutate(Label = str_replace_all(Label,"\x92","'"))->  code_tab
  
  d %>%
    group_by(profile_GOR) %>%
    summarise(N = n()) %>%
    ungroup() %>%
    mutate(prop = N/sum(N,na.rm=T)) %>%
    rename(code_n = profile_GOR) %>%
    right_join(.,code_tab, by="code_n") -> code_tab
  
  m %>%
    spread_draws(b_region[code_n]) %>%
    left_join(., code_tab, by="code_n") %>%
    ggplot(aes(y = fct_reorder(Label, b_region, .fun=mean), x = b_region, fill=prop)) +
    stat_halfeye(.width = c(.95, .5),alpha=1)+
    scale_fill_viridis_c(limits=prop_scal_limit, name="fraction of\nresponses", guide='none') +
    geom_vline(xintercept=0,lty=2) +
    nice_theme +
    coord_cartesian(xlim=x_lim_range)+
    labs(x=expression(mu),y="") +
    ggtitle('Geographical region') -> pl_gor
  
  
  # ------------------
  # marital status
  # 8 lines
  
  codes %>%
    filter(question=="profile_marital_stat",
           !is.na(code_n)) %>%
    mutate(Label=ifelse(code_n==8,"Skipped / not asked", Label)) %>%
    group_by(Label) %>%
    summarise(code_n = mean(code_n))  %>%
    mutate(Label = str_replace_all(Label,"\x92","'"))->  code_tab
  
  d %>%
    group_by(profile_marital_stat) %>%
    summarise(N = n()) %>%
    ungroup() %>%
    mutate(prop = N/sum(N,na.rm=T)) %>%
    rename(code_n = profile_marital_stat) %>%
    right_join(.,code_tab, by="code_n") -> code_tab
  
  m %>%
    spread_draws(b_marital_stat[code_n]) %>%
    left_join(., code_tab, by="code_n") %>%
    ggplot(aes(y = fct_reorder(Label, b_marital_stat, .fun=mean), x = b_marital_stat, fill=prop)) +
    stat_halfeye(.width = c(.95, .5),alpha=1)+
    scale_fill_viridis_c(limits=prop_scal_limit, name="fraction of\nresponses", guide='none') +
    geom_vline(xintercept=0,lty=2) +
    nice_theme +
    coord_cartesian(xlim=x_lim_range)+
    labs(x=expression(mu),y="") +
    ggtitle('Marital status') -> pl_marit
  
  
  # ------------------
  # household income
  # 18 lines
  
  if(type==1){
  codes %>%
    filter(question=="profile_gross_household",
           !is.na(code_n)) %>%
    mutate(Label=ifelse(code_n==18,"Skipped / not asked", Label)) %>%
    group_by(Label) %>%
    summarise(code_n = mean(code_n))  %>%
    mutate(Label = str_replace_all(Label,"\x92","'"))->  code_tab
  
  d %>%
    group_by(profile_gross_household) %>%
    summarise(N = n()) %>%
    ungroup() %>%
    mutate(prop = N/sum(N,na.rm=T)) %>%
    rename(code_n = profile_gross_household) %>%
    right_join(.,code_tab, by="code_n") -> code_tab
  
  m %>%
    spread_draws(b_income[code_n]) %>%
    left_join(., code_tab, by="code_n") %>%
    mutate(Label=factor(Label)) %>%
    ggplot(aes(y = Label, x = b_income, fill=prop)) +
    stat_halfeye(.width = c(.95, .5),alpha=1)+
    scale_fill_viridis_c(limits=prop_scal_limit, name="fraction of\nresponses", guide='none') +
    geom_vline(xintercept=0,lty=2) +
    nice_theme +
    coord_cartesian(xlim=x_lim_range)+
    labs(x=expression(mu),y="") +
    ggtitle('Household income') -> pl_income
  
  }else{
    
    codes %>%
      filter(question=="profile_gross_household",
             !is.na(code_alternative)) %>%
      group_by(alternative_label) %>%
      summarise(code_n = mean(code_alternative))  %>%
      mutate(Label = str_replace_all(alternative_label,"\x92","'"))->  code_tab
    
    code_tab$Label[5] <- "Prefer not to answer/\nDon't know"
    
    d %>%
      group_by(income) %>%
      summarise(N = n()) %>%
      ungroup() %>%
      mutate(prop = N/sum(N,na.rm=T)) %>%
      rename(code_n = income) %>%
      right_join(.,code_tab, by="code_n") -> code_tab
    
    m %>%
      spread_draws(b_income[code_n]) %>%
      left_join(., code_tab, by="code_n") %>%
      mutate(Label=factor(Label)) %>%
      filter(!is.na(Label)) %>%
      ggplot(aes(y = Label, x = b_income, fill=prop)) +
      stat_halfeye(.width = c(.95, .5),alpha=1)+
      scale_fill_viridis_c(limits=prop_scal_limit, name="fraction of\nresponses", guide='none') +
      geom_vline(xintercept=0,lty=2) +
      nice_theme +
      coord_cartesian(xlim=x_lim_range)+
      labs(x=expression(mu),y="") +
      ggtitle('Household income') -> pl_income
  }
  
  
  # ------------------
  library(ggpubr)
  
  if(type==1){
    
  figure <- ggarrange(
    ggarrange(pl_edu, pl_income, ncol=2, nrow=1, labels=c("A","B"), align="h", widths=c(1.2,1),
              common.legend = TRUE, legend='right'),
    ggarrange(NULL,
              ggarrange(pl_grade, pl_gor, ncol = 1, nrow=2, labels = c("C", "D"), align="v", heights=c(1,1.7)), 
              ggarrange(pl_vEU, pl_v2019, ncol = 1, nrow=2, labels = c("E", "F"), align="v", heights=c(1,1.7)),
              ggarrange(NULL, pl_marit,NULL, ncol = 1, nrow=3, labels = c("", "G", ""), align="v", heights=c(0.5,1.5,0.5)),
              NULL,
              ncol=5,nrow=1,widths=c(0.1,1,1.1,1,0.2)
    ),
    align='v',ncol = 1, nrow=2)
  
  }else{
    
    figure <- ggarrange(
      ggarrange(pl_edu, pl_income, ncol=2, nrow=1, labels=c("A","B"), align="h", widths=c(1.2,1),
                common.legend = TRUE, legend='right'),
      ggarrange(NULL,
                ggarrange(pl_grade, pl_gor, ncol = 1, nrow=2, labels = c("C", "D"), align="v", heights=c(1,1.7)), 
                ggarrange(pl_vEU, pl_v2019, ncol = 1, nrow=2, labels = c("E", "F"), align="v", heights=c(1,1.7)),
                ggarrange(NULL, pl_marit,NULL, ncol = 1, nrow=3, labels = c("", "G", ""), align="v", heights=c(0.5,1.5,0.5)),
                NULL,
                ncol=5,nrow=1,widths=c(0.1,1,1.1,1,0.2)
      ),
      align='v',ncol = 1, nrow=2, heights=c(1,1.7))
    
  }
  
  figure <- annotate_figure(figure, top = text_grob(str_c(sel_Q_text,"\n"), face = "bold", size = 12, color="black"))
  #figure <- annotate_figure(figure, top = text_grob(str_c(Q_text[1]), size = 8, color="black"))
  
  #figure
  ggsave(str_c(ordinal_y[QN],"_randomMu_c0",type,".pdf"),plot=figure,width=12,height=14)
  return(figure)
}


# ------------------------------------------------------------------------------------ #
# custom function fior marginal effects plots

make_single_stackedarea_plot_c01 <- function(QN, type=1) {
  
  library(scales)
  
  # ggplot theme
  themeXstack <-
    theme_bw() + theme(
      text = element_text(family = "Helvetica", size = 9),
      #panel.border=element_blank(),
      panel.border = element_rect(size = .4),
      strip.background = element_rect(
        fill = "white",
        color = "white",
        size = 0
      ),
      strip.text = element_text(size = rel(0.8)),
      panel.grid.major.x = element_blank(),
      #panel.grid.major.y=element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(
        size = .3,
        linetype = 1,
        color = "dark grey"
      ),
      #panel.grid.minor.y=element_line(size=.2, linetype=1,color="dark grey"),
      axis.line.x = element_line(size = .4),
      axis.line.y = element_line(size = .4),
      axis.text.x = element_text(size = 7, color =
                                   "black"),
      axis.text.y = element_text(size = 7, color =
                                   "black"),
      axis.line = element_line(size = .4),
      axis.ticks = element_line(color = "black"),
      #legend.position = "bottom",
      legend.position=c(.47,.2),
      legend.background=element_blank(),
      legend.direction = "vertical",
      legend.text = element_text(size=5),
      legend.key.size = unit(0.25, 'cm'),
      plot.title = element_text(size=6),
      plot.subtitle = element_text(size=6),
      legend.margin = margin(r=0.1,l=0.1,t=0,b=0)
    )
  
  
  
  # load model
  ordinal_y <-
    c(
      "att_pubs",
      "att_schools",
      "att_social",
      "behav_meetings",
      "mask_stores",
      "mask_transport",
      "vaccine",
      "vaccine2"
    )
  var_y <- ordinal_y[QN]
  
  m <-  readRDS(sprintf("../regression_results/ordinal_%s.RDS",ordinal_y[QN]))
  
  # print(m, pars=c("b_X","b_age","b_age2","b_covid_affected","b_gender","sigma_vote2019","sigma_voteEUref","sigma_edu","sigma_marital_stat","sigma_income"),probs=c(0.025, 0.975),digits=3)
  # NS: "mask_transport"
  
  # load data
  d <- read_delim("../data/data_allfmt.csv")
  codes <- read_delim("../data/codebook_stan_model.csv")
  codelabels <- read_delim("../data/code_book.csv")
  
  # get point estimates
  m %>%
    spread_draws(c[cutpoint]) %>%
    group_by(cutpoint) %>%
    summarise(c_value = mean(c)) %>%
    .$c_value -> cutpoints
  
  m %>%
    spread_draws(b_X, b_age, b_age2, b_covid_affected, b_d) %>%
    select(starts_with("b_")) %>%
    with(., apply(., 2, mean)) -> fixef
  
  Q_text <- c(
    "Closing of bars, pubs and restaurants",
    "Closing of schools",
    "Restrictions on social meetings",
    "N. times met others outside housedhold/support bubble",
    "Worn mask at supermarket/grocery store",
    "Worn mask in public transport",
    "How likely to accept vaccine",
    "How likely to accept vaccine"
  )
  sel_Q_text <- Q_text[QN]
  
  ## answer labels
  codes %>%
    filter(code_question==!!var_y) %>%
    arrange(code_n) %>%
    pull(Label) %>%
    str_replace_all("\x92","'") -> x_lab
  
  
  # ordered_logistic stacked plot
  N_cat <- max(d[,var_y],na.rm=T)
  grn <- 250
  zLogM <- seq(-5,4,length.out=grn)
  P <- matrix(NA, nrow=grn, ncol=N_cat)
  for(i in 1:grn){
    P[i,] <- ordered_logistic(fixef['b_X']*zLogM[i],cutpoints)
  }
  
  if(QN==1 | QN==4 | QN==5 | QN==8){
    axis_labels <- c('fraction of responses','M-ratio')
  }else{
    axis_labels <- c(' ',' ')
  }
  
  data.frame(P) -> P
  colnames(P) <- x_lab #str_c("p_",1:N_cat)
  P$zlogM <- zLogM
  P$logM <- P$zlogM * 2*sd(log(d$M_cov)) + mean(log(d$M_cov))
  P$M <- exp(P$logM)
  
  library(viridis)
  P %>%
    pivot_longer(cols=x_lab, #starts_with("p_"),
                 values_to="p",
                 names_to="k") %>%
    mutate(k = as_factor(k),
           k = fct_rev(k)) %>%
    ggplot(aes(x=M, y=p, fill=k)) +
    scale_fill_viridis_d(name="",labels = wrap_format(20), option="rocket") +
    geom_area(alpha=0.8 , size=0.4, colour="black")+
    themeXstack +
    labs(y=axis_labels[1],x=axis_labels[2])+
    scale_y_continuous(breaks=seq(0,1,0.1))+
    coord_cartesian(xlim=c(0.5,1.3),ylim=c(0+0.045,1-0.045))+
    ggtitle(label=str_wrap(sel_Q_text,22)) -> pl_stack
  
  #
  return(pl_stack)
}
