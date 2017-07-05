library(ggplot2)
library(reshape2)
## p(theta|k, alpha, beta)
densbeta <- function(theta, n, k, alpha=1, beta=1) {
  stopifnot((theta < 1 || theta > 0) )
  
  x <- dbeta(theta, alpha+k, beta+n-k)
  return(x)
}

## P(theta <= theta_0|k, alpha, beta)
probbeta <- function(theta, n, k, alpha=1, beta=1) {
  stopifnot(theta >= 0)
  stopifnot(theta <= 1)
  
  # x <- pbeta(theta, alpha+k, beta+n-k)
  x <- cumsum(densbeta(theta, alpha, beta, n, k))
  return(x)
}

#' Reimplementation of PROBNML function in SAS
#' @description computes the probability that an observation from a 
#'              binomial distribution Bin(n,theta) will be less than or equal to k. 
#'
#' @param theta is the probability of success for the 
#'              binomial distribution, where 0<=theta<=1. 
#'              In terms of acceptance sampling, is the 
#'              probability of selecting a nonconforming item.
#' @param n is the number of independent Bernoulli trials in the 
#'          binomial distribution, where n>=1. In terms of acceptance 
#'          sampling, n is the number of items in the sample.
#' @param k is the number of successes, where 0<=k<=n. In terms of acceptance 
#'          sampling, is the number of nonconforming items. 
#'
#' @return computes the probability that an observation from a 
#'         binomial distribution Bin(n,theta) will be less than or equal to k. 
#' @export
probbnml <- function(theta, n, k) {
  stopifnot(theta >= 0)
  stopifnot(theta <= 1)
  stopifnot(n >= 1)
  
  x <- sum(sapply(0:k, function(j) choose(n, j)*(theta**j)*((1-theta)**(n-j))))
  return(x)
}

crit_k <- function(theta, prob, n, alpha=1, beta=1, type=c("go", "nogo")) {
  stopifnot(theta >= 0)
  stopifnot(theta <= 1)
  stopifnot(prob >= 0)
  stopifnot(prob <= 1)
  stopifnot(type %in% c("go", "nogo"))
  
  x <- sapply(0:n, function(k) densbeta(theta, n, k, alpha, beta))
  x <- x/length(x)
  x <- cumsum(x)
  
  k <- ifelse(type=="go", 
              which.min(x <= prob),
              which.max(x >= 1-prob))
  return(k)
}

plotPower <- function(theta, go_power, nogo_power) {
  
  data <- data.frame(theta=theta, NoGo=no_go_pwr, Go=go_pwr,
                     Indecisive=indecisive_pwr)
  data.molten <- melt(data, id.vars = "theta")
  colnames(data.molten) <- c("theta", "Function", "Power")
  ggplot(data.molten, aes(x=theta, y=Power, color=Function)) +
    geom_line() +
    scale_x_continuous(breaks=seq(0,1,by=.1)) +
    scale_y_continuous(breaks=seq(0,1,by=.1))
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

n      <- 12 ## sample size
alpha  <- .5
beta   <- .5
theta  <- .3
p_go   <- .9
p_nogo <- .9

## Compute probability P(k|theta) and estimate critical k for which
## the following conditions are fullfilled
## 1) "Go" decision iff P(theta <= theta_go|k_go) < 1 - p_go
## 2) "No Go" decision iff P(theta <= theta_nogo|k_nogo) < p_nogo
##
## Condition 1 can be solved by searching the greatest k with 
## cumsum(dbeta(theta_go, alpha+k, beta+n-k)) <= 1 - p_go
##
## Condition 2 can be solved by searching the greatest k with 
## cumsum(dbeta(theta_nogo, alpha+k, beta+n-k)) >= p_nogo
k_nogo <- crit_k(theta=theta, prob=p_nogo, n=n, alpha=alpha, beta=beta, type="nogo")
k_go <- crit_k(theta=theta, prob=p_go, n=n, alpha=alpha, beta=beta, type="go")


## Go Power
theta <- seq(0, 1, by=.001)
go_pwr <- sapply(theta, function(theta) probbnml(theta, n, k_go))
go_pwr <- 1-go_pwr
plot(go_pwr)

## No Go Power
no_go_pwr <- sapply(theta, function(theta) probbnml(theta, n, k_nogo)) 
plot(no_go_pwr)

## Indecisive Power
indecisive_pwr <- 1 - no_go_pwr - go_pwr
plot(indecisive_pwr)

plotPower(theta, go_pwr, no_go_pwr)

# ## Sanity check
# a <- no_go_pwr[which.max(indecisive_pwr)]
# b <- go_pwr[which.max(indecisive_pwr)]
# c <- max(indecisive_pwr)
# sum(a, b, c)


sprintf("n=%i", n)
sprintf("p_go=%.2f", p_go)
sprintf("p_nogo=%.2f", p_nogo)
sprintf("k_go=%i", k_go)
sprintf("k_nogo=%i", k_nogo)

sprintf("Go (k>=%i)", k_go)
sprintf("No Go (k<=%i)", k_nogo)
sprintf("Indecisive (%i<=k<=%i)", k_nogo, k_go)

prob_go <- cumsum(sapply(0:n, function(k) densbeta(theta = .3, n, k, alpha, beta)))
prob_go <- prob_go/length(prob_go)
plot(prob_go)
print(k_go)

prob_nogo <- cumsum(sapply(0:n, function(k) densbeta(theta = .3, n, k, alpha, beta)))
prob_nogo <- prob_nogo/length(prob_nogo)
plot(prob_nogo)
print(k_nogo)

data <- data.frame(k=0:n,
                   Go=prob_go,
                   NoGo=prob_nogo)
data.molten <- melt(data, id.vars = "k")
colnames(data.molten) <- c("k", "Function", "Probability")


font <- list(
  # family = "sans serif",
  size = 14,
  color = '#EBEBE9')

k_nogo_text <- paste0("NoGo Crit. (", k_nogo, "|",
                      round(data$NoGo[k_nogo], 4), ")")

k_nogo_val <- list(
  x = k_nogo,
  y = data$NoGo[k_nogo+1],
  text = paste0("NoGo Crit. (", k_nogo, "|",
                round(data$NoGo[k_nogo], 4), ")"),
  xref = "x",
  yref = "y",
  showarrow = TRUE,
  arrowhead = 7,
  arrowcolor = "#C5C9CB",
  ax = 50,
  ay = -40
)

k_go_text <- paste0("Go Crit. (", k_go, "|",
                    round(data$Go[k_go], 4), ")")

k_go_val <- list(
  x = k_go,
  y = data$Go[k_go],
  text = k_go_text,
  xref = "x",
  yref = "y",
  showarrow = TRUE,
  arrowhead = 7,
  arrowcolor = "#C5C9CB",
  ax = -50,
  ay = -40
)

annotations <- list(k_nogo_val,
                    k_go_val)

print(annotations)

plot_ly(data.molten, x = ~k, y= ~Probability,
               type = 'scatter', mode = 'lines',
               line = list(width = 5),
               color = ~Function) %>%
  layout(title = 'Cummulative Densities',
         hovermode="all",
         xaxis = list(title = 'Number of Successes',
                      range = c(0, max(data$k)),
                      showgrid = F),
         yaxis = list(title = 'Probability',
                      range = c(0, 1),
                      showgrid = F),
         # legend = list(orientation = 'h'),
         legend = list(x = 0.8, y = 0.5),
         annotations = annotations,
         font=font,
         plot_bgcolor="transparent",
         paper_bgcolor="transparent")



