# k <- 25 ## number of successes
# success_prob <- .1 ## success probability
# alpha <- .5 ## prior param
# beta  <- .5 ## prior param
# 
# theta_u <- .8
# theta_l <- .8
# p_u <- .8
# p_l <- .7

library(ggplot2)
library(reshape2)
## p(theta|k, alpha, beta)
densbeta <- function(theta, n, k, alpha=1, beta=1) {
  stopifnot(n >= k)
  stopifnot((theta < 1 || theta > 0) )
  
  x <- dbeta(theta, alpha+k, beta+n-k)
  return(x)
}

## P(theta <= theta_0|k, alpha, beta)
probbeta <- function(theta, n, k, alpha=1, beta=1) {
  stopifnot(n >= k)
  stopifnot(theta >= 0)
  stopifnot(theta <= 1)
  
  # x <- pbeta(theta, alpha+k, beta+n-k)
  x <- cumsum(densbeta(theta, alpha, beta, n, k))
  return(x)
}

crit_k <- function(theta, prob, n, alpha=1, beta=1, type=c("go", "nogo")) {
  stopifnot(theta >= 0)
  stopifnot(theta <= 1)
  stopifnot(prob >= 0)
  stopifnot(prob <= 1)
  stopifnot(type %in% c("go", "nogo"))
  
  args <- list(theta=theta, alpha=alpha, beta=beta, n=n, k=1:n)
  x <- do.call(probbeta, args)
  x <- x/length(x)
  
  k <- ifelse(type=="go", 
              which.min(x < 1-prob),
              which.max(x >= prob))
  return(k)
}


#####
# thetas <- seq(0,1,by=0.01)
# args <- list(theta=thetas, 
#              alpha=alpha, beta=beta, 
#              n=n, k=k)
# ## plot f1 over all thetas
# x <- do.call(f1, args)
# x <- x/length(thetas)
# plot(thetas, x)
# ## plot F1 = P(theta <= theta_0|k, alpha, beta) over all thetas
# x <- do.call(F1, args)
# plot(thetas, x)
# ## plot 1-F1 = 1-P(theta <= theta_0|k, alpha, beta) over all thetas
# x <- 1-do.call(F1, args)
# plot(thetas, x)
# 
# theta <- .5
# ks <- 1:n
# args <- list(theta=theta, 
#              alpha=alpha, beta=beta, 
#              n=n, k=ks)
# ## plot f1 over all k's, with 1<=k<=n
# x <- do.call(f1, args)
# x <- x/length(x)
# plot(ks, x)
# ## plot F1 = P(theta <= theta_0|k, alpha, beta) over all over all k's, with 1<=k<=n
# x <- do.call(F1, args)
# x <- x/length(x)
# plot(ks, x)
# ## plot 1-F1 = 1-P(theta <= theta_0|k, alpha, beta) over all over all k's, with 1<=k<=n
# # x <- 1-do.call(F1, args)
# # plot(ks, x)
#####


n <- 100 ## sample size
alpha <- 1
beta  <- 1
k_go   <- crit_k(theta=.5, prob=.5, n=n, type="go")
k_nogo <- crit_k(theta=.8, prob=.5, n=n, type="nogo")

theta <- seq(0, 1, by=.01)
args <- list(theta=theta, 
             alpha=alpha, beta=beta, 
             n=n, k=k_go)
prob_go <- do.call(probbeta, args)
prob_go <- prob_go/length(prob_go)

args <- list(theta=theta, 
             alpha=alpha, beta=beta, 
             n=n, k=k_nogo)
prob_nogo <- do.call(probbeta, args)
prob_nogo <- prob_nogo/length(prob_nogo)
prob_nogo <- prob_nogo


data <- data.frame(theta=theta, prob_go=prob_go, prob_nogo=prob_nogo)
data.molten <- melt(data, id.vars = "theta")
colnames(data.molten) <- c("theta", "Function", "Probability")
ggplot(data.molten, aes(x=theta, y=Probability, color=Function)) + geom_line()




