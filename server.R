
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(reshape2)
library(ggplot2)
library(plotly)

betaMode <- function(alpha, beta) {
  return((alpha - 1)/(alpha+beta-2))
}

betaMean <- function(alpha, beta) {
  return((alpha)/(alpha+beta))
}

betaStd <- function(alpha, beta) {
  return(sqrt((alpha * beta)/(((alpha + beta)**2) * (alpha + beta + 1))))
}

posterior2 <- function(theta, n, k, alpha, beta) {
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
  
  data <- data.frame(k, 
                     Upper=upper_cdf,
                     Lower=lower_cdf)
  data.molten <- melt(data, id.vars = "k")
  colnames(data.molten) <- c("k", "Function", "Probability")
  
  fig <- ggplot(data.molten, aes(x=k, y=Probability, color=Function)) + 
    geom_line() +
    geom_vline(xintercept = upper_crit, 
               color="#F78E87", linetype="dashed") +
    geom_hline(yintercept = upper_cdf[upper_crit],
               color="#F78E87", linetype="dashed") +
    geom_vline(xintercept = lower_crit,
               color="#2FC9CD", linetype="dashed") +
    geom_hline(yintercept = lower_cdf[lower_crit],
               color="#2FC9CD", linetype="dashed")
  return(fig)
}

plotPosteriorCDFs2 <- function(k, 
                               upper_cdf, lower_cdf, 
                               upper_crit, lower_crit) {
  data <- data.frame(k, 
                     Upper=upper_cdf,
                     Lower=lower_cdf)
  data.molten <- melt(data, id.vars = "k")
  colnames(data.molten) <- c("k", "Function", "Probability")
  
  fig <- plot_ly(data.molten, x = ~k, y= ~Probability,
                 type = 'scatter', mode = 'lines', 
                 color = ~Function) %>%
    layout(title = 'Cummulative Densities',
           hovermode="all",
           xaxis = list(title = 'Number of Successes'),
           yaxis = list(title = 'Probability',
                        range = c(0, 1)),
           legend = list(orientation = 'h'))
  
  return(fig)
}

plotPowerCurves <- function(thetas, 
                            upper_power, lower_power, overall_power,
                            upper_crit, lower_crit) {
  
  data <- data.frame(theta=thetas, 
                     NoGo=upper_power,
                     Go=lower_power,
                     Indecisive=overall_power)
  
  data.molten <- melt(data, id.vars = "theta")
  colnames(data.molten) <- c("theta", "Function", "Power")
  
  fig <- ggplot(data.molten, aes(x=theta, y=Power, color=Function)) + 
    geom_line()
  
  return(fig)
}

plotPowerCurves2 <- function(thetas, 
                            upper_power, lower_power, overall_power,
                            upper_crit, lower_crit) {
  
  data <- data.frame(theta=thetas, 
                     NoGo=upper_power,
                     Go=lower_power,
                     Indecisive=overall_power)
  
  data.molten <- melt(data, id.vars = "theta")
  colnames(data.molten) <- c("theta", "Function", "Power")
  
  fig <- plot_ly(data.molten, x = ~theta, y= ~Power,
                 type = 'scatter', mode = 'lines', 
                 color = ~Function) %>%
    layout(title = 'Power Curves',
           hovermode="all",
           xaxis = list(title = 'Theta',
                        tick0 = 0,
                        dtick = 0.1),
           yaxis = list(title = 'Power',
                        range = c(0, 1)),
           legend = list(orientation = 'h'))
  
  return(fig)
}


shinyServer(function(input, output) {
  
  ## reactive vars
  x <- reactive({
    set.seed(42)
    x <- rbinom(input$n, size = 1, input$pi)
  })
  
  output$samplingInput <- renderUI({
    n <- input$n
    
    probInput <- numericInput("pi", withMathJax("$$\\textbf{Success}\\ \\textbf{Probability}\\ \\pi$$"), 
                              min = 0, max = 1, step=.1, value = .5)
    successInput <- numericInput("k", withMathJax("$$\\textbf{Number of Successes}\\ k$$"), 
                                 min = 0, max = n, step=1, value = floor(n/2))
    
    uiElement <- switch(input$radioSample,
                        "prob" = probInput,
                        "successes" = successInput,
                        probInput)
    return(uiElement)
    
  })
  
  ## plots
  output$triPlot <- renderPlotly({
    
    prob <- input$prob ## success probability
    x <- x()
    n <- length(x)
    
    alpha <- input$alpha
    beta  <- input$beta
    pi <- seq(0, 1, length.out=1000)
    
    ## Data
    k <- switch(input$radioSample,
                        "prob" = sum(x),
                        "successes" = input$k,
                        probInput)
    
    # k <- sum(x) ## number of successes
    
    # Likelihood p(x|pi_test) with x ~ Bin(pi_test, alpha_test, beta)
    likelihood <- dbinom(k, n, pi)
    likelihood <- likelihood/(sum(likelihood)/length(pi)) ## Normalize Density

    ## Prior p(pi) based on Beta(pi, alpha, beta)
    prior <- dbeta(pi, alpha, beta)
    # prior <- prior/(sum(prior[!is.infinite(prior)])/length(pi)) ## Normalize Density
    
    ## Posterior Distribution p(pi|x)
    posterior <- dbeta(pi, alpha+k, beta+n-k)
    # posterior <- posterior/(sum(posterior)/length(pi)) ## Normalize Density
    
    data <- data.frame(Pi=pi, 
                       Posterior=posterior, Likelihood=likelihood, Prior=prior)
    data.molten <- melt(data, id.vars = "Pi")
    colnames(data.molten) <- c("Pi", "Function", "Density")
    
    fig <- plot_ly(data.molten, x = ~Pi, y= ~Density,
                   type = 'scatter', mode = 'lines', 
                   color = ~Function) %>%
      layout(title = 'Density Functions',
             hovermode="all",
             xaxis = list(title = 'Pi',
                          tick0 = 0,
                          dtick = 0.1),
             yaxis = list(title = 'Density',
                          range = c(0, max(data$Posterior, data$Likelihood))),
             legend = list(orientation = 'h'))
    
    return(fig)
  })
  
  output$posteriorProbPlot <- renderPlotly({
    
    prob <- input$prob ## success probability
    x <- x()
    n <- length(x)
    
    alpha <- input$alpha
    beta  <- input$beta
    pi <- seq(0, 1, length.out=1000)
    
    ## Data
    k <- switch(input$radioSample,
                "prob" = sum(x),
                "successes" = input$k,
                probInput)
    
    ## Posterior Distribution p(pi|x)
    posterior <- dbeta(pi, alpha+k, beta+n-k)
    posterior <- posterior/sum(posterior)
    go_prob   <- cumsum(posterior)
    nogo_prob <- 1-cumsum(posterior)
    
    data <- data.frame(Pi=pi,
                       Go=go_prob, 
                       NoGo=nogo_prob)
    data.molten <- melt(data, id.vars = "Pi")
    colnames(data.molten) <- c("Pi", "Function", "Probability")
    
    fig <- plot_ly(data.molten, x = ~Pi, y= ~Probability,
                   type = 'scatter', mode = 'lines', 
                   color = ~Function) %>%
      # add_trace(x = c(.3, .6), y=c(.3, .6), mode = "lines") %>%
      layout(title = 'Probability Functions',
             shapes=list(type='line', 
                         x0=0.2, x1=0.2, #
                         y0=0.3, y1=0.3, 
                         line=list(dash='dot', width=1)),
             hovermode="all",
             xaxis = list(title = 'Pi',
                          tick0 = 0,
                          dtick = 0.1),
             yaxis = list(title = 'Probability'),
             legend = list(orientation = 'h'))
    
    return(fig)
  })
  
  output$posteriorCDFPlot <- renderPlot({
    
    n <- input$n
    alpha <- input$alpha
    beta <- input$beta
    theta_u <- input$pi_u
    theta_l <- input$pi_l
    p_u <- input$p_u
    p_l <- input$p_l
    
    ## Compute critical values
    ## vars for critical value calculation
    upper_cdf  <- 1-posteriorCDF(1:n, theta_u, alpha, beta) # P(theta > theta_u|n) = 1 - probbeta(theta, alpha+k, beta+n-k)
    lower_cdf  <- posteriorCDF(1:n, theta_l, alpha, beta) # P(theta < theta_l|n) = probbeta(theta_l, alpha+k, beta+n-k)
    upper_crit <- max(which(upper_cdf >= p_u))
    lower_crit <- min(which(lower_cdf >= p_l))
    fig <- plotPosteriorCDFs(1:n, upper_cdf, lower_cdf, upper_crit, lower_crit)
    return(fig)
  })
  
  output$posteriorCDFPlot2 <- renderPlotly({
    
    n <- input$n
    alpha <- input$alpha
    beta <- input$beta
    theta_u <- input$pi_u
    theta_l <- input$pi_l
    p_u <- input$p_u
    p_l <- input$p_l
    
    ## Compute critical values
    ## vars for critical value calculation
    upper_cdf  <- 1-posteriorCDF(1:n, theta_u, alpha, beta) # P(theta > theta_u|n) = 1 - probbeta(theta, alpha+k, beta+n-k)
    lower_cdf  <- posteriorCDF(1:n, theta_l, alpha, beta) # P(theta < theta_l|n) = probbeta(theta_l, alpha+k, beta+n-k)
    upper_crit <- max(which(upper_cdf >= p_u))
    lower_crit <- min(which(lower_cdf >= p_l))
    
    print(upper_crit)
    fig <- plotPosteriorCDFs2(1:n, upper_cdf, lower_cdf, upper_crit, lower_crit)
    return(fig)
  })
  
  output$powerCurvePlot <- renderPlot({
    
    n <- input$n
    alpha <- input$alpha
    beta  <- input$beta
    theta_u <- input$pi_u
    theta_l <- input$pi_l
    p_u <- input$p_u
    p_l <- input$p_l
    
    upper_cdf <- 1-posteriorCDF(1:n, theta_u, alpha, beta) # P(theta > theta_u|n) = 1 - probbeta(theta, alpha+k, beta+n-k)
    lower_cdf <- posteriorCDF(1:n, theta_l, alpha, beta) # P(theta < theta_l|n) = probbeta(theta_l, alpha+k, beta+n-k)
    upper_crit <- max(which(upper_cdf >= p_u))
    lower_crit <- min(which(lower_cdf >= p_l))
    
    thetas <- seq(0, 1, length.out=100)
    upper_power <- posterior2(thetas, n, upper_crit, alpha, beta)
    upper_power <- upper_power/length(thetas)
    upper_power <- 1-cumsum(upper_power)
    
    lower_power <- posterior2(thetas, n, lower_crit, alpha, beta)
    lower_power <- lower_power/length(thetas)
    lower_power <- cumsum(lower_power)
    
    overall_power <- lower_power- (1-upper_power)
    overall_power <- overall_power/length(thetas)
    overall_power <- overall_power/sum(overall_power)
    fig <- plotPowerCurves(thetas, 
                           upper_power, lower_power, overall_power,
                           upper_crit, lower_crit)
    return(fig)
  })
  
  output$powerCurvePlot2 <- renderPlotly({
    
    n <- input$n
    alpha <- input$alpha
    beta  <- input$beta
    theta_u <- input$pi_u
    theta_l <- input$pi_l
    p_u <- input$p_u
    p_l <- input$p_l
    
    upper_cdf <- 1-posteriorCDF(1:n, theta_u, alpha, beta) # P(theta > theta_u|n) = 1 - probbeta(theta, alpha+k, beta+n-k)
    lower_cdf <- posteriorCDF(1:n, theta_l, alpha, beta) # P(theta < theta_l|n) = probbeta(theta_l, alpha+k, beta+n-k)
    upper_crit <- max(which(upper_cdf >= p_u))
    lower_crit <- min(which(lower_cdf >= p_l))
    
    thetas <- seq(0, 1, length.out=100)
    upper_power <- posterior2(thetas, n, upper_crit, alpha, beta)
    upper_power <- upper_power/length(thetas)
    upper_power <- 1-cumsum(upper_power)
    
    lower_power <- posterior2(thetas, n, lower_crit, alpha, beta)
    lower_power <- lower_power/length(thetas)
    lower_power <- cumsum(lower_power)
    
    overall_power <- lower_power- (1-upper_power)
    overall_power <- overall_power/length(thetas)
    overall_power <- overall_power/sum(overall_power)
    fig <- plotPowerCurves2(thetas, 
                           upper_power, lower_power, overall_power,
                           upper_crit, lower_crit)
    return(fig)
  })
  
  ## latex elements
  output$priorDistFormula <- renderUI({
    alpha <- input$alpha
    beta <- input$beta
    uiElement <- withMathJax(helpText(sprintf('$$\\begin{align}
                                      p(\\theta)&=\\frac{\\theta^{\\alpha-1}(1-\\theta)^{\\beta-1}}{B(\\alpha, \\beta)}
                                            \\\\&=Beta(\\theta|\\alpha, \\beta)
                                            \\\\&=Beta(\\theta|\\textbf{%.2f}, \\textbf{%.2f})
                                      \\end{align}$$', alpha, beta, alpha, beta, alpha, beta)))
    return(uiElement)
  })
  
  output$likelihoodFormula <- renderUI({
    n <- input$n
    k <- sum(x())
    alpha <- input$alpha
    beta <- input$beta
    uiElement <- list(withMathJax(helpText(sprintf('$$X\\sim Bin(n, \\theta) = Bin(\\textbf{%i}, \\theta)$$', n))),
                      withMathJax(helpText(sprintf('$$\\begin{align}
                                      p(x|\\theta)&={n\\choose{x}}\\theta^{x}(1-\\theta)^{n-x}
                                              \\\\&={\\textbf{%i}\\choose{\\textbf{%i}}}\\theta^{\\textbf{%i}}(1-\\theta)^{\\textbf{%i}}
                                      \\end{align}$$', n, k, k, n-k))))
    return(uiElement)
  })
  
  output$posteriorFormula <- renderUI({
    n <- input$n
    k <- sum(x())
    alpha <- input$alpha
    beta <- input$beta
    uiElement <- withMathJax(helpText(sprintf('$$\\begin{align}
                                               p(\\theta|x)&=p(x|\\theta)p(\\theta)
                                                       \\\\&=\\theta^{x}(1-\\theta)^{n-x}\\theta^{\\alpha-1}(1-\\theta)^{\\beta-1}
                                                       \\\\&=\\theta^{(\\alpha+x)-1}(1-\\theta)^{(\\beta+n-x)-1}
                                                       \\\\&=Beta(\\theta|\\alpha+x, \\beta+n-x)
                                                       \\\\&=Beta(\\theta|\\textbf{%.2f}, \\textbf{%.2f})
                                                   \\end{align}$$', sum(alpha, k), sum(beta, n, -k))))
    return(uiElement)
  })
  
  
  ## point estimate calculations
  output$pointEst_Prior <- renderTable({
    validate(
      need(!(input$n==0), NULL)
    )
    
    alpha <- input$alpha
    beta  <- input$beta
    
    beta_mode <- betaMode(alpha, beta)
    beta_mean <- betaMean(alpha, beta)
    beta_std  <- betaStd(alpha, beta)
    
    pE <- data.frame(Type=c("Mode", "Mean", "Std"), 
                     PointEstimate=c(beta_mode, beta_mean, beta_std))
    return(pE)
  })
  
  output$pointEst_Likelihood <- renderTable({
    x <- x()
    n <- input$n
    validate(
      need(!(input$n==0), NULL)
    )
    
    k <- switch(input$radioSample,
                "prob" = sum(x),
                "successes" = input$k,
                probInput)
    
    alpha <- k + 1
    beta <- n - k + 1
    
    beta_mode <- betaMode(alpha, beta)
    beta_mean <- betaMean(alpha, beta)
    beta_std  <- betaStd(alpha, beta)
    
    pE <- data.frame(Type=c("Mode", "Mean", "Std"), 
                     PointEstimate=c(beta_mode, beta_mean, beta_std))
    return(pE)
  })
  
  output$pointEst_Posterior <- renderTable({
    x <- x()
    n <- input$n
    validate(
      need(n!=0, NULL)
    )
    
    k <- switch(input$radioSample,
                "prob" = sum(x),
                "successes" = input$k,
                probInput)
    
    alpha <- input$alpha + k
    beta <- input$beta + n - k
    
    beta_mode <- betaMode(alpha, beta)
    beta_mean <- betaMean(alpha, beta)
    beta_std  <- betaStd(alpha, beta)
    
    pE <- data.frame(Type=c("Mode", "Mean", "Std"), 
                     PointEstimate=c(beta_mode, beta_mean, beta_std))
    return(pE)
    
  })
  
})
