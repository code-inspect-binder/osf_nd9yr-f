# -------------------------------------------------------------
# Estiamte from meta-d analysis for reporting
# -------------------------------------------------------------
rm(list=ls())
setwd("/mnt/sda2/matteoHDD/git_local_HDD/covid19-misinformation/jags/")

library(mlisi)
library(rjags)
library(tidyverse)
library(tidybayes)

## load jags results 
output <- readRDS("hdmeta_jags_both_MLok.RDS")
str(output)

# subset as in supplemental analysis 
# (excluding respondents who believe in covid bioweapon theory) 
# "hdmeta_jags_both_MLok_subset01.RDS"
#output <- readRDS("hdmeta_jags_both_MLok_subset01.RDS")
#   par       value .lower .upper .width .point .interval
#   <chr>     <dbl>  <dbl>  <dbl>  <dbl> <chr>  <chr>    
# 1 diffLogM -0.237 -0.319 -0.152   0.95 mean   hdci     
# 2 Mratio_1  1.03   0.961  1.10    0.95 mean   hdci     
# 3 Mratio_2  0.813  0.768  0.860   0.95 mean   hdci 

# visualize posteriors
plot(output[,"mu_logMratio_1"])
plot(output[,"mu_logMratio_2"])
plot(output[,"sigma_logMratio_1"])
plot(output[,"sigma_logMratio_2"])

# plot chains on separate facets, as in supplemental material
output %>%
  gather_draws(mu_logMratio_1,mu_logMratio_2,sigma_logMratio_1,sigma_logMratio_2) %>%
  group_by(.variable, .chain) %>%
  ggplot(aes(x=.iteration, y=.value, color=as.factor(.chain))) +
  geom_line(alpha=0.5,lwd=0.4) +
  facet_grid(.variable~.chain, scale='free_y') +
  geom_smooth(method='loess',lwd=0.6) + 
  scale_color_viridis_d(option = "C", begin=0.2, end=0.8)+
  labs(color='chain')


# group-level estimates and credibility intervals
output %>%
  spread_draws(mu_logMratio_1, 
               mu_logMratio_2) %>%
  mutate(diffLogM = mu_logMratio_2 - mu_logMratio_1,
         Mratio_1 = exp(mu_logMratio_1),
         Mratio_2 = exp(mu_logMratio_2)) %>%
  pivot_longer(c(Mratio_1, Mratio_2, diffLogM),
               values_to="value",
               names_to="par") %>%
  select(par,value) %>%
  group_by(par) %>%
  mean_hdci()

# calculate Bayes factor
# note we used informative prior for both logM:
#   dnorm(-0.3, 0.3)
# corresponds to 95% Pr. in [0.41, 1,33] 
# thus SD of prior difference is 
# sqrt(0.3^2 + 0.3^2) = 0.4242641
output %>%
  spread_draws(mu_logMratio_1, 
               mu_logMratio_2) %>%
  mutate(diffLogM = mu_logMratio_2 - mu_logMratio_1) %>%
  pull(diffLogM) -> diffLogM

# function to calculate Savage-Dickey BF (modified from mlisi library)
savage.dickey.bf <- function (x, x_0 = 0, prior.mean = 0, prior.sd = 1, plot = F, xlabel="parameter value") 
{
  require(polspline)
  fit.posterior <- logspline(x)
  posterior_w <- dlogspline(x_0, fit.posterior)
  if (plot) {

    plot(fit.posterior, n=500, xlab = xlabel, ylab = "density", 
         lwd = 2, xlim = c(prior.mean - 2*prior.sd, prior.mean + 2*prior.sd))
    x <- seq(prior.mean - 2*prior.sd, prior.mean + 2*prior.sd, length.out = 500)
    lines(x, dnorm(x, mean = prior.mean, sd = prior.sd), 
          col = "red", lwd = 2)
    abline(v = x_0, lty = 2)
    points(x_0, posterior_w, pch = 19, col = "black")
    points(x_0, dnorm(x_0, prior.mean, prior.sd), pch = 19, 
           col = "red")
    legend("topright", c("posterior", "prior"), lwd = 2, 
           col = c("black", "red"), pch = 19, bty = "n", inset = 0.02)
  }
  cat(paste0("Approximate BF (Savage-Dickey) in favor of null x=", 
             x_0, " : ", round(posterior_w/dnorm(x_0, prior.mean, 
                                                 prior.sd), digits = 2), "\n"))
  invisible(posterior_w/dnorm(x_0, prior.mean, prior.sd))
}

BF01 <- savage.dickey.bf(diffLogM, prior.sd=sqrt(0.3^2 + 0.3^2), plot=T, xlabel=expression(paste(Delta["log M-ratio"]," (general science minus COVID-19)")))
1/BF01



