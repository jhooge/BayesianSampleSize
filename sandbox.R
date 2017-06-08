library(reshape2)
library(ggplot2)

n <- 100 ## sample size
success_prob <- .1 ## success probability
alpha <- 1 ## prior param
beta  <- 1 ## prior param

theta_u <- .5
theta_l <- .5
p_u <- .5
p_l <- .5

posterior <- function(theta, n, success_prob, alpha, beta) {
  # x <- rbinom(n, size = 1, success_prob)
  # k <- sum(x)
  k <- round(n*success_prob, 0)

  posterior <- dbeta(theta, alpha+k, beta+n-k)
  return(posterior)
}

posterior2 <- function(theta, n, k, alpha, beta) {
  n <- as.integer(n)
  k <- as.integer(k)
  # stopifnot(n>=k, "n>=k")
  
  posterior <- dbeta(theta, alpha+k, beta+n-k)
  return(posterior)
}

posteriorCDF <- function(k, theta, alpha, beta) {
  n <- max(k) ## sample size
  ## P(theta > theta_u|n) = 1 - probbeta(theta, alpha+k, beta+n-k)
  posterior_cdf <- sapply(1:n, function(k) posterior2(theta, n, k, alpha, beta))
  posterior_cdf <- posterior_cdf/sum(posterior_cdf)
  posterior_cdf <- cumsum(posterior_cdf)
  
  return(posterior_cdf)
}

plotPosteriorCDFs <- function(k, 
                              upper_cdf, lower_cdf, 
                              upper_crit, lower_crit) {
  
  data <- data.frame(k=1:n, 
                     Upper=upper_cdf,
                     Lower=lower_cdf)
  data.molten <- melt(data, id.vars = "k")
  colnames(data.molten) <- c("k", "Function", "Probability")
  
  fig <- ggplot(data.molten, aes(x=k, y=Probability, color=Function)) + 
    geom_line() +
    geom_vline(xintercept = c_u, 
               color="#F78E87", linetype="dashed") +
    geom_hline(yintercept = upper_posterior_pwr[upper_crit],
               color="#F78E87", linetype="dashed") +
    geom_vline(xintercept = c_l,
               color="#2FC9CD", linetype="dashed") +
    geom_hline(yintercept = lower_posterior_pwr[lower_crit],
               color="#2FC9CD", linetype="dashed")
  return(fig)
}

plotPowerCurves <- function(thetas, 
                            upper_power, lower_power, overall_power,
                            upper_crit, lower_crit) {
  
  data <- data.frame(theta=thetas, 
                     UpperPower=upper_power,
                     LowerPower=lower_power,
                     OverallPower=overall_power)
  
  data.molten <- melt(data, id.vars = "theta")
  colnames(data.molten) <- c("theta", "Function", "Power")
  
  fig <- ggplot(data.molten, aes(x=theta, y=Power, color=Function)) + 
    geom_line()
  
  return(fig)
}

## Compute critical values
## vars for critical value calculation
upper_cdf <- 1-posteriorCDF(1:n, theta_u, alpha, beta) # P(theta > theta_u|n) = 1 - probbeta(theta, alpha+k, beta+n-k)
lower_cdf <- posteriorCDF(1:n, theta_l, alpha, beta) # P(theta < theta_l|n) = probbeta(theta_l, alpha+k, beta+n-k)
upper_crit <- max(which(upper_cdf >= p_u))
lower_crit <- min(which(lower_cdf >= p_l))
cdf_fig <- plotPosteriorCDFs(1:n, upper_cdf, lower_cdf, upper_crit, lower_crit)

## Compute Power Curves at the determined critical value over theta
## vars for power calculation
thetas <- seq(0, 1, length.out=100)
upper_power <- posterior2(thetas, n, c_u, alpha, beta)
upper_power <- upper_power/length(thetas)
upper_power <- 1-cumsum(upper_power)

lower_power <- posterior2(thetas, n, c_l, alpha, beta)
lower_power <- lower_power/length(thetas)
lower_power <- cumsum(lower_power)

overall_power <- lower_power- (1-upper_power)
overall_power <- overall_power/length(thetas)
overall_power <- overall_power/sum(overall_power)
power_fig <- plotPowerCurves(thetas, 
                             upper_power, lower_power, overall_power,
                             upper_crit, lower_crit) 

cdf_fig
power_fig
