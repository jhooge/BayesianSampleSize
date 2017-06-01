library(reshape2)
library(ggplot2)

n <- 48 ## sample size
n_test <- n/2
n_ref  <- n/2
# postprobu1 <- .8 ## upper posterior probability for scenario1
# postprobl1 <- .7 ## lower posterior probability for scenario1

pi_test    <- .5
pi_ref     <- .5 ## reference pi-value used for power calculations
alpha_test <- 1 ## value of a1 of the Beta(a1,b1)-prior for pi1
beta_test  <- 1 ## value of b1 of the Beta(a1,b1)-prior for pi1
alpha_ref  <- 1 ## value of a2 of the Beta(a2,b2)-prior for pi2
beta_ref   <- 1 ## value of b2 of the Beta(a2,b2)-prior for pi2

set.seed(42)
## Binomial distributed data
x_test <- sum(rbinom(n_test, 1, 1-pi_test))
x_ref  <- sum(rbinom(n_ref, 1, pi_ref))

compute_lower_crit <- function(n, k, 
                               alpha_test, alpha_ref,
                               beta_test, beta_ref) {
  N <- sum(1:n, 
            alpha_test, alpha_ref,
            beta_test, beta_ref, -2)
  M <- sum(n/2, alpha_test, beta_test, -1)
  n <- sum(k, alpha_test, alpha_ref, -1)
  crit <- cumsum(dhyper(n, N_, N_-M_, k)) 
  return(crit)
}

compute_upperCrit <- function(n, k, 
                              alpha_test, alpha_ref,
                              beta_test, beta_ref) {
  upper_crit <- compute_Crit(n, k, 
                             alpha_test, alpha_ref, 
                             beta_test, beta_ref)
  crit <- 1-compute_lower_crit
  return(crit)
}

compute_Crit(15, 5, alpha_test, alpha_ref, beta_test, beta_ref)

plot_phyper <- function(n, k, alpha_test, alpha_ref, beta_test, beta_ref) {
  N_ <- sum(n, 
            alpha_test, alpha_ref,
            beta_test, beta_ref, -2)
  M_ <- sum(n/2, alpha_test, beta_test, -1)
  n_ <- sum(k, alpha_test, alpha_ref, -1)
  
  data <- data.frame("N"=1:n, 
                     "upper_crit"=1-cumsum(dhyper(1:n, N_, N_-M_, k)),
                     "lower_crit"=cumsum(dhyper(1:n, N_, N_-M_, k)))
  data.molten <- melt(data, id.vars="N")
  colnames(data.molten) <- c("N", "Function", "Value")
  ggplot(data.molten) + 
    geom_line(aes(x=N, y=Value, colour=Function)) +
    ggtitle("Critical Values")
}


plot(cumsum(dbeta(seq(0,1,length.out=100),3,3)))


## Helper functions for power calculations
OR <- function(p1, p2) {
  y <- (p1*(1-p2))/(p2*(1-p1))
  return(y)
}

right_factor <- function(p1, p2, n1, n2, k) {
  p1 <- rep()
  
  y <- prod((1-p1)**n1, (1-p2)**n2, (p2/(1-p2))**k)
  return(y)
}

power_fun <- function(n_test, n_ref, x_test, k, p1, p2) {
  f <- choose(n_test, x_test)*choose(n_ref, k - x_test)*OR(p1, p2)
  return(f)
}


## Power calculation
p1 <- seq(0,1,length.out = 100)
# p2 <- seq(0,1,lenght.out = 100)
p2 <- .2
y <- power_fun(n_test, n_ref, x_test, k, p1, p2)
plot(p1, y)







n <- 40
p1 <- .5
k = qbinom(0.05, n, p1)
pbinom(k, n, p1)
pbinom(k - 1, n, p1)
k = qbinom(0.05, n, p1) - 1 ## rejection region for given n
pbinom(k, n, p1) ## power when prob=p1




# prod(choose(n_test, x_test),
#      choose(n_ref, k - x_test),
#      OR(p1, p2)**x_test,
#      right_factor(p1, p2, n/2, n/2, k))





# set.seed(42)
# ## Prior distribution 'pi_test' for a Go, No-Go decision (responder probability of test drug)
# plot(density(rbeta(1:n/2, alpha_test, beta_test)))
# mean_beta_test <- 1/(1-(beta_test/alpha_test))
# 
# set.seed(100)
# n <- 9
# z <- rbetabinom(1000, 0.5, size=n, theta=4)
# par(las=1,bty="l")
# plot(table(z)/length(z),ylim=c(0,0.34),col="gray",lwd=4,
#      ylab="probability")
# points(0:n,dbinom(0:n,size=n,prob=0.5),col=2,pch=16,type="b")
# points(0:n,dbetabinom(0:n,size=n,theta=4,
#                       prob=0.5),col=4,pch=17,type="b")
# ## correspondence with SuppDists 
# if (require(SuppDists)) {
#   d1a <- dghyper(0:5,a=-5,N=-10,k=5)
#   d1b <- dbetabinom(0:5,shape1=5,shape2=5,size=5)
#   max(abs(d1a-d1b))
#   p1a <- pghyper(0:5,a=-5,N=-10,k=5,lower.tail=TRUE)
#   p1b <- cumsum(d1b)
#   max(abs(p1a-p1b))
# } 
# 
# 
# 
# 
# 
# ## Prior distribution 'pi_ref' for a Go, No-Go decision (responder probability for reference drug)
# plot(density(rbeta(1:n/2, alpha_ref, alpha_ref)))