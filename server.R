
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

shinyServer(function(input, output) {
  
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
    go_prob <- cumsum(posterior)
    nogo_prob <- 1-cumsum(posterior)
    
    data <- data.frame(Pi=pi,
                       Go=go_prob, 
                       NoGo=nogo_prob)
    data.molten <- melt(data, id.vars = "Pi")
    colnames(data.molten) <- c("Pi", "Function", "Probability")
    
    fig <- plot_ly(data.molten, x = ~Pi, y= ~Probability,
                   type = 'scatter', mode = 'lines', 
                   color = ~Function) %>%
      layout(title = 'Probability Functions',
             hovermode="all",
             xaxis = list(title = 'Pi',
                          tick0 = 0,
                          dtick = 0.1),
             yaxis = list(title = 'Probability'),
             legend = list(orientation = 'h'))
    
    return(fig)
  })
  
  
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
