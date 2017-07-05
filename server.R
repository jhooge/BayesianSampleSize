
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

posterior <- function(theta, n, k, alpha=1, beta=1) {
  x <- dbeta(theta, alpha+k, beta+n-k)
  return(x)
}

## P(theta <= theta_0|k, alpha, beta)
posterior_cdf <- function(theta, n, k, alpha=1, beta=1) {
  stopifnot(theta >= 0)
  stopifnot(theta <= 1)
  
  # x <- pbeta(theta, alpha+k, beta+n-k)
  x <- cumsum(posterior(theta, alpha, beta, n, k))
  return(x)
}

densbeta <- function(theta, n, k, alpha=1, beta=1) {
  stopifnot((theta < 1 || theta > 0) )
  
  x <- dbeta(theta, alpha+k, beta+n-k)
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
  k <- k-1 ## k starts at 0
  
  # print(paste0("Type=", type))
  # print(sprintf("k=%i",k))
  # print(sprintf("p_type=%.2f", prob))
  # print(x)
  
  return(k)
}

plotPosteriorCDF <- function(k, go_cdf, nogo_cdf,
                             k_go, k_nogo) {

  data <- data.frame(k,
                     go=go_cdf,
                     nogo=nogo_cdf)
  data.molten <- melt(data, id.vars = "k")
  colnames(data.molten) <- c("k", "Function", "Probability")

  font <- list(
    # family = "sans serif",
    size = 14,
    color = '#EBEBE9')

  k_nogo_text <- paste0("NoGo Crit. (", k_nogo, "|",
                        round(nogo_cdf[k_nogo], 4), ")")

  k_nogo_val <- list(
    x = k_nogo,
    y = lower_cdf[k_nogo],
    text = paste0("Lower Crit. (", k_nogo, "|",
                  round(nogo_cdf[k_nogo], 4), ")"),
    xref = "x",
    yref = "y",
    showarrow = TRUE,
    arrowhead = 7,
    arrowcolor = "#C5C9CB",
    ax = 50,
    ay = -40
  )

  k_go_text <- paste0("Upper Crit. (", k_go, "|",
                            round(go_cdf[k_go], 4), ")")

  k_go_val <- list(
    x = k_go,
    y = go_cdf[k_go],
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

  fig <- plot_ly(data.molten, x = ~k, y= ~Probability,
                 type = 'scatter', mode = 'lines',
                 line = list(width = 5),
                 color = ~Function
                 ) %>%
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

  return(fig)
}


plotPowerCurves <- function(thetas,
                            go_power, nogo_power, indecisive_power) {

  data <- data.frame(theta=thetas,
                     Go=go_power,
                     NoGo=nogo_power,
                     Indecisive=indecisive_power)

  data.molten <- melt(data, id.vars = "theta")
  colnames(data.molten) <- c("theta", "Function", "Power")

  font <- list(
    # family = "sans serif",
    size = 14,
    color = '#EBEBE9')

  fig <- plot_ly(data.molten, x = ~theta, y= ~Power,
                 type = 'scatter', mode = 'lines',
                 line = list(width = 5),
                 color = ~Function) %>%
    layout(title = 'Power Curves',
           hovermode="all",
           xaxis = list(title = 'θ',
                        tick0 = 0,
                        dtick = 0.1,
                        showgrid = F),
           yaxis = list(title = 'Power',
                        dtick = 0.1,
                        # range = c(0, 1),
                        showgrid = F),
           legend = list(orientation = 'h'),
           font = font,
           plot_bgcolor="transparent",
           paper_bgcolor="transparent")

  return(fig)
}

plotDensities <- function(theta, posterior, likelihood, prior) {

  data <- data.frame(Theta=theta,
                     Posterior=posterior, Likelihood=likelihood, Prior=prior)
  data.molten <- melt(data, id.vars = "Theta")
  colnames(data.molten) <- c("Theta", "Function", "Density")

  font <- list(
    # family = "sans serif",
    size = 14,
    color = '#EBEBE9')

  fig <- plot_ly(data.molten, x = ~Theta, y= ~Density,
                 type = 'scatter', mode = 'lines',
                 line = list(width = 5),
                 color = ~Function) %>%
    layout(title = 'Density Functions',
           hovermode="all",
           xaxis = list(title = "θ",
                        tick0 = 0,
                        dtick = 0.1,
                        showgrid = F),
           yaxis = list(title = 'Density',
                        range = c(0, max(data$Posterior, data$Likelihood)),
                        showgrid = F),
           legend = list(orientation = 'h'),
           font=font,
           plot_bgcolor="transparent",
           paper_bgcolor="transparent")

  return(fig)
}

shinyServer(function(input, output) {

  ## reactive vars
  x <- reactive({
    set.seed(42)
    x <- rbinom(input$n, size = 1, input$pi)
  })
  
  dens <- reactive({
    n <- input$n 
    m <- 101
    
    alpha <- input$alpha
    beta <- input$beta
    thetas <- seq(0, 1, length.out=m)
    
    d <- c()
    for (theta in thetas) {
      for (k in 1:n) {
        d <- c(d, posterior(theta, n, k, alpha, beta))
      } 
    }
    
    d <- matrix(d, nrow = n, ncol=m)
    # d <- apply(d, 1, function(x) x/sum(x))
    
    return(d)
  })

  output$samplingInput <- renderUI({
    n <- input$n

    probInput <- numericInput("pi", withMathJax("$$\\textbf{Success}\\ \\textbf{Probability}\\ \\theta$$"),
                              min = 0, max = 1, step=.1, value = .5,
                              width = "40%")
    successInput <- numericInput("k", withMathJax("$$\\textbf{Number of Successes}\\ k$$"),
                                 min = 0, max = n, step=1, value = floor(n/2),
                                 width = "40%")

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

    fig <- plotDensities(pi, posterior, likelihood, prior)

    return(fig)
  })

  output$posteriorProbPlot <- renderPlotly({

    prob <- input$prob ## success probability
    x <- x()
    n <- length(x)

    alpha <- input$alpha
    beta  <- input$beta
    theta <- seq(0, 1, length.out=1000)

    ## Data
    k <- switch(input$radioSample,
                "prob" = sum(x),
                "successes" = input$k,
                probInput)

    ## Posterior Distribution p(theta|x)
    posterior <- dbeta(theta, alpha+k, beta+n-k)
    posterior <- posterior/sum(posterior)
    go_prob   <- cumsum(posterior)
    nogo_prob <- 1-cumsum(posterior)

    data <- data.frame(Theta=theta,
                       Go=go_prob,
                       NoGo=nogo_prob)
    data.molten <- melt(data, id.vars = "Theta")
    colnames(data.molten) <- c("Theta", "Function", "Probability")

    fig <- plot_ly(data.molten, x = ~Theta, y= ~Probability,
                   type = 'scatter', mode = 'lines',
                   line = list(width = 5),
                   color = ~Function) %>%
      # add_trace(x = c(.3, .6), y=c(.3, .6), mode = "lines") %>%
      layout(title = 'Probability Functions',
             shapes=list(type='line',
                         x0=0.2, x1=0.2, #
                         y0=0.3, y1=0.3,
                         line=list(dash='dot', width=1)),
             hovermode="all",
             xaxis = list(title = 'theta',
                          tick0 = 0,
                          dtick = 0.1),
             yaxis = list(title = 'Probability'),
             legend = list(orientation = 'h'))

    return(fig)
  })
  
  output$critValPlot <- renderPlotly({
    
    n <- input$n
    theta <- input$theta_0
    alpha <- input$alpha
    beta  <- input$beta
    p_go <- input$p_go
    p_nogo <- input$p_nogo
    
    validate(need(sum(p_go, p_nogo) >= 1, "The sum of Go and NoGo probabilities should be larger or equal to 1."))
    
    k_go   <- crit_k(theta=theta, prob=p_go, n=n, alpha=alpha, beta=beta, type="go")
    k_nogo <- crit_k(theta=theta, prob=p_nogo, n=n, alpha=alpha, beta=beta, type="nogo")
    
    prob_go <- cumsum(sapply(0:n, function(k) densbeta(theta, n, k, alpha, beta)))
    prob_go <- prob_go/length(prob_go)
    
    prob_nogo <- cumsum(sapply(0:n, function(k) densbeta(theta, n, k, alpha, beta)))
    prob_nogo <- 1-prob_nogo/length(prob_nogo)
    
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
                          round(data$NoGo[k_nogo+1], 4), ")")
    
    k_nogo_val <- list(
      x = k_nogo,
      y = data$NoGo[k_nogo+1],
      text = paste0("NoGo Crit. (", k_nogo, "|",
                    round(data$NoGo[k_nogo+1], 4), ")"),
      xref = "x",
      yref = "y",
      showarrow = TRUE,
      arrowhead = 7,
      arrowcolor = "#C5C9CB",
      ax = 50,
      ay = -40
    )
    
    k_go_text <- paste0("Go Crit. (", k_go, "|",
                        round(data$Go[k_go+1], 4), ")")
    
    k_go_val <- list(
      x = k_go,
      y = data$Go[k_go+1],
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
    
    fig <- plot_ly(data.molten, x = ~k, y= ~Probability,
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
             legend = list(orientation = 'h'),
             # legend = list(x = 0.8, y = 0.5),
             annotations = annotations,
             font=font,
             plot_bgcolor="transparent",
             paper_bgcolor="transparent")
    
    return(fig)
  })

  output$posterior3DPlot <- renderPlotly({
    
    dens <- dens()
    
    n <- nrow(dens)
    m <- ncol(dens)
    
    print(dim(dens))
    thetas <- seq(0, 1, length.out=m)
    
    font <- list(
      # family = "sans serif",
      size = 14,
      color = '#EBEBE9')
    
    fig <- plot_ly(x=thetas, y=1:nrow(dens), z=dens, showscale=FALSE, source="posterior3D") %>% 
      add_surface() %>%
      layout(title = 'Density',
             hovermode="all",
             scene = list(
               xaxis = list(
                            title = "θ",
                            tick0 = 0,
                            dtick = 0.2,
                            # range = c(0, 1),
                            showgrid = T),
               yaxis = list(title = 'k',
                            # range = c(0, n),
                            showgrid = T),
               zaxis = list(title = 'Density',
                            range = c(0, 20),
                            showgrid= F)),
             font=font,
             plot_bgcolor="transparent",
             paper_bgcolor="transparent",
             showlegend=FALSE)
    
    return(fig)
  })
  
  output$posterior2D_k <- renderPlotly({
    
    dens <- dens()
    n <- nrow(dens)
    m <- ncol(dens)
    
    hover_event <- event_data("plotly_hover", source="posterior3D")
    click_event <- event_data("plotly_click", source="posterior3D")
    
    validate(
      need(!is.null(click_event), "Please click in the density plot above to select a projection.")
    )
    
    theta <- click_event$x
    thetas <- seq(0, 1, length.out=m)
    
    i <- match(theta, thetas)
    dens <- dens[, i]
    dens <- dens/sum(dens)
    
    data <- data.frame(k=1:n,
                       Density=dens)
    
    font <- list(
      # family = "sans serif",
      size = 14,
      color = '#EBEBE9')
    
    fig <- plot_ly(data, x = ~k, y= ~Density,
                   type = 'scatter', mode = 'lines',
                   line = list(width = 5,
                               color="#510E61")) %>%
      layout(title = sprintf("For θ=%.2f", theta),
             hovermode="all",
             xaxis = list(title = 'k',
                          tick0 = 1,
                          showgrid = F),
             yaxis = list(title = 'Density',
                          showgrid = F),
             legend = list(orientation = 'h'),
             font=font,
             plot_bgcolor="transparent",
             paper_bgcolor="transparent",
             showlegend=FALSE)
    return(fig)
  })
  
  output$posterior2D_theta <- renderPlotly({
    
    dens <- dens()
    n <- nrow(dens)
    m <- ncol(dens)
    hover_event <- event_data("plotly_hover", source="posterior3D")
    click_event <- event_data("plotly_click", source="posterior3D")
    
    validate(
      need(!is.null(click_event), "Please click in the density plot above to select a projection.")
    )
    
    
    k <- click_event$y
    thetas <- seq(0, 1, length.out=m)
    
    dens <- dens[k, ]
    dens <- dens/sum(dens)
    
    # dens <- dens/m ## scale based on theta resolution
    
    data <- data.frame(theta=thetas,
                       Density=dens)
    
    font <- list(
      # family = "sans serif",
      size = 14,
      color = '#EBEBE9')
    
    fig <- plot_ly(data, x = ~theta, y= ~Density,
                   type = 'scatter', mode = 'lines',
                   line = list(width = 5,
                               color="#510E61")) %>%
      layout(title = sprintf("For k=%i", k),
             hovermode="all",
             xaxis = list(title = 'θ',
                          tick0 = 0,
                          dtick = 0.1,
                          showgrid = F),
             yaxis = list(title = 'Density',
                          showgrid = F),
             legend = list(orientation = 'h'),
             font=font,
             plot_bgcolor="transparent",
             paper_bgcolor="transparent",
             showlegend=FALSE)
    return(fig)
  })
  
  output$powerCurvePlot <- renderPlotly({

    n <- input$n
    alpha <- input$alpha
    beta  <- input$beta
    theta <- input$theta_0
    p_go <- input$p_go
    p_nogo <- input$p_nogo
    
    validate(need(sum(p_go, p_nogo) >= 1, "The sum of Go and NoGo probabilities should be larger or equal to 1."))
    
    ## Compute probability P(k|theta) and estimate critical k for which
    ## the following conditions are fullfilled
    ## 1) "Go" decision iff P(theta <= theta_0|k_go) < 1 - p_go
    ## 2) "No Go" decision iff P(theta <= theta_0|k_nogo) < p_nogo
    ##
    ## Condition 1 can be solved by searching the greatest k with 
    ## cumsum(dbeta(theta_0, alpha+k, beta+n-k)) <= 1 - p_go
    ##
    ## Condition 2 can be solved by searching the greatest k with 
    ## cumsum(dbeta(theta_0, alpha+k, beta+n-k)) >= p_nogo
    k_go   <- crit_k(theta=theta, prob=p_go, 
                     n=n, alpha=alpha, beta=beta, 
                     type="go")
    k_nogo <- crit_k(theta=theta, prob=p_nogo, n=n, 
                     alpha=alpha, beta=beta, 
                     type="nogo")
    
    print(k_go)
    print(k_nogo)
    
    ## Go Power
    theta <- seq(0, 1, by=.001)
    go_pwr <- sapply(theta, function(theta) probbnml(theta, n, k_go))
    # go_pwr <- 1 - go_pwr
    go_pwr <- 1 - go_pwr
    
    ## No Go Power
    # no_go_pwr <- sapply(theta, function(theta) probbnml(theta, n, k_nogo)) 
    no_go_pwr <- sapply(theta, function(theta) probbnml(theta, n, k_nogo)) 
    
    ## Indecisive Power
    # indecisive_pwr <- 1 - no_go_pwr - go_pwr
    indecisive_pwr <- 1 - no_go_pwr - go_pwr
    
    # ## Sanity check
    a <- no_go_pwr[which.max(indecisive_pwr)]
    b <- go_pwr[which.max(indecisive_pwr)]
    c <- max(indecisive_pwr)
    print("Sane?")
    print(sum(a, b, c) == 1)
    
    fig <- plotPowerCurves(theta, go_pwr, no_go_pwr, indecisive_pwr)
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
