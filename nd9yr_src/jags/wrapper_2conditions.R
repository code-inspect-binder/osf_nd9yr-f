#####################################

# Estimate metacognitive efficiency (Mratio) at the group level
#
# Adaptation in R of matlab function 'fit_meta_d_mcmc_groupCorr.m'
# by Steve Fleming 
# for more details see Fleming (2017). HMeta-d: hierarchical Bayesian 
# estimation of metacognitive efficiency from confidence ratings. 
#
# you need to install the following packing before using the function:
# tidyverse
# magrittr
# reshape2
# rjags
# coda
# lattice
# broom
# ggpubr
# ggmcmc
#
# nR_S1 and nR_S2 should be two lists of each nR_S1 or nR_S2 per task
# model output is a large mcmc list and two vectors for d1 and c1
#
# Author: Nicolas Legrand nicolas.legrand@cfin.au.dk

# modified 2021 Matteo Lisi

#####################################

## Packages
library(tidyverse)
#library(magrittr)
#library(reshape2)
library(rjags)
library(coda)
#library(lattice)
#library(broom)
# library(ggpubr)
# library(ggmcmc)

# Model -------------------------------------------------------------------

metad_2conditions <- function (nR_S1_tot, nR_S2_tot) {
  # nR_S1_tot and nR_S2_tot should be lists of size N subjects.
  # Each list contain 2 responses vectors.
  # Design matrix: Condition 1: (0) - Condition 2: (1)
  
  # Type 1 parameters
  nratings <- length(nR_S1_tot[[1]][[1]])/2
  nsubj <- length((nR_S1_tot))
  
  d1 <- matrix(ncol = 2, nrow = nsubj)
  c1 <- matrix(ncol = 2, nrow = nsubj)
  counts_total = array(dim = c(nsubj, nratings*4, 2))
  
  for (n in 1:(nsubj)) {
    for (i in 1:2) {
      
      nR_S1 = nR_S1_tot[[n]][i]
      nR_S2 = nR_S2_tot[[n]][i]
      
      # Adjust to ensure non-zero counts for type 1 d' point estimate
      adj_f <- 1/((nratings)*2)
      nR_S1_adj = unlist(nR_S1) + adj_f
      nR_S2_adj = unlist(nR_S2) + adj_f
      
      ratingHR <- matrix()
      ratingFAR <- matrix()
      
      for (c in 2:(nratings*2)) {
        ratingHR[c-1] <- sum(nR_S2_adj[c:length(nR_S2_adj)]) / sum(nR_S2_adj)
        ratingFAR[c-1] <- sum(nR_S1_adj[c:length(nR_S1_adj)]) / sum(nR_S1_adj)
      }
      
      d1[n, i] = qnorm(ratingHR[(nratings)]) - qnorm(ratingFAR[(nratings)])
      c1[n, i] = -0.5 * (qnorm(ratingHR[(nratings)]) + qnorm(ratingFAR[(nratings)]))
      
      counts_total[n, , i] <- t(nR_S1[[1]]) %>% 
        cbind(t(nR_S2[[1]]))
    }
  }
  
  Tol <- 1e-05
  
  data_f <- list(
    d1 = d1,  # [Subjects * Condition]
    c1 = c1,  # [Subjects * Condition]
    nsubj = nsubj,
    counts = counts_total,  # [Subjects * Counts * Condition]
    nratings = nratings,
    Tol = Tol
  )
  
  # Model using runJAGS
  # hdmeta_2cond_type2_ML.txt
  meta2_model <- jags.model("hdmeta_2cond_type2_ML.txt", data = data_f,
                            n.chains = 3, n.adapt= 1000, quiet = FALSE)
  update(meta2_model, n.iter=4000)

  # Sampling
  output <- coda.samples(
    model          = meta2_model,
    variable.names = c("mu_logMratio_1","mu_logMratio_2",
                       "Mratio",
                       "cS1","cS2",
                       "sigma_logMratio_1","sigma_logMratio_2"),
    n.iter         = 30000,
    thin           = 9 )
  
  return(output)
}



# -------------------------------------------------------------
# additional helper functions

# extract individual parameters and store for further analyses
save_par_i <- function (d_prep, output) {
  
  nR_S1_tot <- list()
  nR_S2_tot <- list()
  for(n in 1:d_prep$J){
    nR_S1_tot[[n]] <- list(unname(d_prep$F1[n,]), unname(d_prep$F2[n,]))
    nR_S2_tot[[n]] <- list(unname(d_prep$T1[n,]), unname(d_prep$T2[n,]))
  }
  
  # nR_S1_tot and nR_S2_tot should be lists of size N subjects.
  # Each list contain 2 responses vectors.
  # Design matrix: Condition 1: (0) - Condition 2: (1)
  
  # Type 1 parameters - hmeta-d style
  nratings <- length(nR_S1_tot[[1]][[1]])/2
  nsubj <- length((nR_S1_tot))
  
  d1 <- matrix(ncol = 2, nrow = nsubj)
  c1 <- matrix(ncol = 2, nrow = nsubj)
  counts_total = array(dim = c(nsubj, nratings*4, 2))
  
  for (n in 1:(nsubj)) {
    for (i in 1:2) {
      
      nR_S1 = nR_S1_tot[[n]][i]
      nR_S2 = nR_S2_tot[[n]][i]
      
      # Adjust to ensure non-zero counts for type 1 d' point estimate
      adj_f <- 1/((nratings)*2)
      nR_S1_adj = unlist(nR_S1) + adj_f
      nR_S2_adj = unlist(nR_S2) + adj_f
      
      ratingHR <- matrix()
      ratingFAR <- matrix()
      
      for (c in 2:(nratings*2)) {
        ratingHR[c-1] <- sum(nR_S2_adj[c:length(nR_S2_adj)]) / sum(nR_S2_adj)
        ratingFAR[c-1] <- sum(nR_S1_adj[c:length(nR_S1_adj)]) / sum(nR_S1_adj)
      }
      
      d1[n, i] = qnorm(ratingHR[(nratings)]) - qnorm(ratingFAR[(nratings)])
      c1[n, i] = -0.5 * (qnorm(ratingHR[(nratings)]) + qnorm(ratingFAR[(nratings)]))
      
      counts_total[n, , i] <- t(nR_S1[[1]]) %>% 
        cbind(t(nR_S2[[1]]))
    }
  }
  
  # Now extract Type 2 parameters from output
  par_i_sm <- summary(output[,grep('Mratio',varnames(output))]) # updated 
  par_i_sm <- par_i_sm$statistics[grep("Mratio",rownames(par_i_sm$statistics)),]
  M_i <- matrix(ncol = 2, nrow = nsubj)
  for (n in 1:(nsubj)) {
    for (i in 1:2) {
      M_i[n, i] <- par_i_sm[paste("Mratio[",n,",",i,"]",sep=""),1]
    }
  }
  
  # compile
  d_out <- data.frame(cbind(d1,c1,M_i))
  colnames(d_out) <- c("d_sci","d_cov","c_sci","c_cov","M_sci","M_cov")
  d_out$ID <- 1:nrow(d_out) # row/ID order was preserved throrough 
  return(d_out)
}


