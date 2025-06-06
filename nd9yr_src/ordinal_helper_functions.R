# Miscellaneous helper functions

# wrapper function to extract posterior samples (not used since I discovered tidybayes package)
extract_rstan <- function(fit, name, dim1=NULL, dim2=NULL, dim3=NULL){
  par_name <- case_when(
    is.null(dim1)  ~ name,
    !is.null(dim1)  ~ str_c(name,"[",as.character(dim1),"]")
  )
  X <- extract(fit, pars = par_name, inc_warmup = FALSE)[[1]]
  return(c(X))
}

# invlogit <- function (x) {1/(1+exp(-x))} # same as plogis()

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


# calculate prediction for regression model
predicted_p_full_extended_c <- function(m, d, sd2=T){
  
  # c("b_0","b_X","b_age","b_age2","b_covid_affected","b_grade","b_vote2019","b_voteEUref","b_edu")
  
  m %>%
    spread_draws(c[cutpoint])%>%
    count(cutpoint) %>% nrow() -> n_cutpoints
  
  m %>%
    spread_draws(c[cutpoint]) %>%
    pivot_wider(names_from=cutpoint, values_from=c) -> cutpoints
  
  m %>%
    spread_draws(b_d, b_X, b_age, b_age2, b_covid_affected, b_gender) -> fixef
  
  m %>%
    spread_draws(sigma_grade,sigma_vote2019,sigma_voteEUref,sigma_edu,sigma_region, sigma_marital_stat, sigma_income) -> sigmas
  
  m %>%
    spread_draws(b_grade[grade]) %>%
    pivot_wider(names_from=grade, values_from=b_grade) -> b_grade
  
  m %>%
    spread_draws(b_vote2019[vote]) %>%
    pivot_wider(names_from=vote, values_from=b_vote2019) -> b_vote2019
  
  m %>%
    spread_draws(b_voteEUref[vote]) %>%
    pivot_wider(names_from=vote, values_from=b_voteEUref) -> b_voteEUref
  
  m %>%
    spread_draws(b_region[region]) %>%
    pivot_wider(names_from=region, values_from=b_region) -> b_region
  
  m %>%
    spread_draws(b_edu[edu]) %>%
    pivot_wider(names_from=edu, values_from=b_edu) -> b_edu
  
  m %>%
    spread_draws(b_marital_stat[marital_stat]) %>%
    pivot_wider(names_from=marital_stat, values_from=b_marital_stat) -> b_marital_stat
  
  m %>%
    spread_draws(b_income[income]) %>%
    pivot_wider(names_from=income, values_from=b_income) -> b_income
  
  # codes[codes$question=='profile_education_level',]
  
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
  pred_p <- array(NA, dim=c(nrow(cutpoints), nrow(d), n_cutpoints+1))
  
  for(i in 1:nrow(d)){
    gamma <- unlist(unname(fixef[,'b_d'] * d_prime[i] +
                            fixef[,'b_X'] * X[i] +
                             fixef[,'b_gender'] * gender[i] +
                             fixef[,'b_covid_affected']*d$covid_affected[i] +
                             fixef[,'b_age'] * age_c[i] + fixef[,'b_age2'] * age_c[i]^2 +
                             b_grade[,as.character(d$profile_socialgrade_cie[i])]*sigmas[,'sigma_grade'] + 
                             b_region[,as.character(d$profile_GOR[i])]*sigmas[,'sigma_region'] + 
                             b_edu[,as.character(d$education[i])]*sigmas[,'sigma_edu'] + 
                             b_vote2019[,as.character(d$vote2019[i])]*sigmas[,'sigma_vote2019'] + 
                             b_voteEUref[,as.character(d$voteEURef[i])]*sigmas[,'sigma_voteEUref'] +
                             b_marital_stat[,as.character(d$profile_marital_stat[i])]*sigmas[,'sigma_marital_stat']+ 
                             b_income[,as.character(d$income[i])]*sigmas[,'sigma_income']))
    
    for(it in 1:nrow(cutpoints)){
      pred_p[it,i,] <- ordered_logistic(gamma[it], c(unlist(cutpoints[it,c(as.character(1:(n_cutpoints)))])))
    }
    
  }
  
  return(pred_p)
}


#' Highest density interval
#'
#' This is a function that will calculate the highest density interval from a
#' posterior sample.
#'
#' The default is to calcualte the highest 95 percent interval. It can be used
#' with any numeric vector instead of having to use one of the specific MCMC
#' classes. This function has been adapted from John K. Kruschke (2011). Doing
#' Bayesian Data Analaysis: A Tutorial with R and BUGS.
#'
#' @param x Numeric vector of a distribution of data, typically a posterior
#' sample
#' @param prob Width of the interval from some distribution. Defaults to `0.95`.
#' @param warn Option to turn off multiple sample warning message Must be in the
#' range of `[0,1]`.
#' @return Numeric range
#' @export
#' @examples
#' x <- qnorm(seq(1e-04, .9999, length.out=1001))
#' hdi_95 <- hdi(x, .95)
#' hdi_50 <- hdi(x, .50)
#'
#' hist(x, br=50)
#' abline(v=hdi_95, col="red")
#' abline(v=hdi_50, col="green")
#'
#' x <- exp(seq(pi * (1 - (1/16)), pi, len = 1000))
#' x <- c(x, rev(x)[-1])
#' x <- c(-x, x)
#' plot(sort(x), type="l")
#' plot(density(x, adjust=0.25))
#' abline(v=hdi(x, p=.49), col=2)
#' abline(v=hdi(x, p=.50), col=3)
hdi <- function(x, prob=0.95, warn=TRUE) {
  if (anyNA(x)) {
    stop("HDI: ", "x must not contain any NA values.", call.=FALSE)
  }
  
  N <- length(x)
  
  if (N < 3) {
    if (warn) {
      warning("HDI: ", "length of `x` < 3.", " Returning NAs", call.=FALSE)
    }
    return(c(NA_integer_, NA_integer_))
  }
  
  x_sort <- sort(x)
  window_size <- as.integer(floor(prob * length(x_sort)))
  
  if (window_size < 2) {
    if (warn) {
      warning("HDI: ", "window_size < 2.", " `prob` is too small or x does not ",
              "contain enough data points.", " Returning NAs.",
              call.=FALSE)
    }
    return(c(NA_integer_, NA_integer_))
  }
  
  lower <- seq_len(N - window_size)
  upper <- window_size + lower
  
  # vectorized difference between edges of cumulative distribution based on
  # scan_length. Values are arranged from left to right scanning.
  window_width_diff <- x_sort[upper] - x_sort[lower]
  
  # find minimum of width differences, check for multiple minima
  min_i <- which(window_width_diff == min(window_width_diff))
  n_candies <- length(min_i)
  
  if (n_candies > 1) {
    if (any(diff(sort(min_i)) != 1)) {
      if (warn) {
        warning("HDI: ", "Identical densities found along ",
                "different segments of the distribution.", " Choosing rightmost.",
                call.=FALSE)
      }
      min_i <- max(min_i)
    } else {
      min_i <- floor(mean(min_i))
    }
  }
  
  # get values based on minimum
  c(x_sort[min_i], x_sort[upper[min_i]])
}
