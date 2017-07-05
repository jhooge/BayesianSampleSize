
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

# plotPosteriorCDFs <- function(k,
#                               upper_cdf, lower_cdf,
#                               upper_crit, lower_crit) {
# 
#   data <- data.frame(k,
#                      Upper=upper_cdf,
#                      Lower=lower_cdf)
#   data.molten <- melt(data, id.vars = "k")
#   colnames(data.molten) <- c("k", "Function", "Probability")
# 
#   fig <- ggplot(data.molten, aes(x=k, y=Probability, color=Function)) +
#     geom_line() +
#     geom_vline(xintercept = upper_crit,
#                color="#F78E87", linetype="dashed") +
#     geom_hline(yintercept = upper_cdf[upper_crit],
#                color="#F78E87", linetype="dashed") +
#     geom_vline(xintercept = lower_crit,
#                color="#2FC9CD", linetype="dashed") +
#     geom_hline(yintercept = lower_cdf[lower_crit],
#                color="#2FC9CD", linetype="dashed")
#   return(fig)
# }

plotPosteriorCDFs2 <- function(k,
                               upper_cdf, lower_cdf,
                               upper_crit, lower_crit) {

  data <- data.frame(k,
                     Upper=upper_cdf,
                     Lower=lower_cdf)
  data.molten <- melt(data, id.vars = "k")
  colnames(data.molten) <- c("k", "Function", "Probability")

  font <- list(
    # family = "sans serif",
    size = 14,
    color = '#EBEBE9')

  lower_crit_text <- paste0("Lower Crit. (", lower_crit, "|",
                            round(lower_cdf[lower_crit], 4), ")")

  lower_crit_val <- list(
    x = lower_crit,
    y = lower_cdf[lower_crit],
    text = paste0("Lower Crit. (", lower_crit, "|",
                  round(lower_cdf[lower_crit], 4), ")"),
    xref = "x",
    yref = "y",
    showarrow = TRUE,
    arrowhead = 7,
    arrowcolor = "#C5C9CB",
    ax = 50,
    ay = -40
  )

  upper_crit_text <- paste0("Upper Crit. (", upper_crit, "|",
                            round(upper_cdf[upper_crit], 4), ")")

  upper_crit_val <- list(
    x = upper_crit,
    y = upper_cdf[upper_crit],
    text = upper_crit_text,
    xref = "x",
    yref = "y",
    showarrow = TRUE,
    arrowhead = 7,
    arrowcolor = "#C5C9CB",
    ax = -50,
    ay = -40
  )

  annotations <- list(lower_crit_val,
                      upper_crit_val)

  fig <- plot_ly(source = "cdfPlot",
                 data.molten, x = ~k, y= ~Probability,
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

# plotPowerCurves <- function(thetas,
#                             upper_power, lower_power, overall_power,
#                             upper_crit, lower_crit) {
#   ##
#   data <- data.frame(theta=thetas,
#                      NoGo=upper_power,
#                      Go=lower_power,
#                      Indecisive=overall_power)
# 
#   data.molten <- melt(data, id.vars = "theta")
#   colnames(data.molten) <- c("theta", "Function", "Power")
# 
#   fig <- ggplot(data.molten, aes(x=theta, y=Power, color=Function)) +
#     geom_line()
# 
#   return(fig)
# }

plotPowerCurves2 <- function(thetas,
                            upper_power, lower_power, overall_power,
                            upper_crit, lower_crit) {

  data <- data.frame(theta=thetas,
                     Go=upper_power,
                     NoGo=lower_power,
                     Indecisive=overall_power)

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
           xaxis = list(title = 'theta_0',
                        tick0 = 0,
                        dtick = 0.1,
                        showgrid = F),
           yaxis = list(title = 'Power',
                        range = c(0, 1),
                        showgrid = F),
           # legend = list(orientation = 'h')
           legend = list(x = 0.8, y = 0.5),
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
           xaxis = list(title = "theta_0",
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

# plotTrials <- function(trials) {
#
#   fig <- plot_ly(
#     x = c("Successes", "Fails"),
#     y = c(sum(trials), sum(!trials)),
#     type = "bar"
#   )
#
#   return(fig)
# }

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
        d <- c(d, posterior2(theta, n, k, alpha, beta))
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

  ## plots

  # output$trialPlot <- renderPlotly({
  #   x <- x()
  #   fig <- plotTrials(x)
  # })

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

  # output$posteriorCDFPlot <- renderPlot({
  # 
  #   n <- input$n
  #   alpha <- input$alpha
  #   beta <- input$beta
  #   theta_u <- input$pi_u
  #   theta_l <- input$pi_l
  #   p_u <- input$p_u
  #   p_l <- input$p_l
  # 
  #   ## Compute critical values
  #   ## vars for critical value calculation
  #   upper_cdf  <- 1-posteriorCDF(1:n, theta_u, alpha, beta) # P(theta > theta_u|n) = 1 - probbeta(theta, alpha+k, beta+n-k)
  #   lower_cdf  <- posteriorCDF(1:n, theta_l, alpha, beta) # P(theta < theta_l|n) = probbeta(theta_l, alpha+k, beta+n-k)
  #   upper_crit <- max(which(upper_cdf >= p_u))
  #   lower_crit <- min(which(lower_cdf >= p_l))
  #   fig <- plotPosteriorCDFs(1:n, upper_cdf, lower_cdf, upper_crit, lower_crit)
  #   return(fig)
  # })

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
    upper_cdf  <- posteriorCDF(1:n, theta_u, alpha, beta) # P(theta > theta_u|n) = 1 - probbeta(theta, alpha+k, beta+n-k)
    lower_cdf  <- 1-posteriorCDF(1:n, theta_l, alpha, beta) # P(theta < theta_l|n) = probbeta(theta_l, alpha+k, beta+n-k)
    upper_crit <- min(which(upper_cdf >= p_u))
    lower_crit <- max(which(lower_cdf >= p_l))

    fig <- plotPosteriorCDFs2(1:n, upper_cdf, lower_cdf, upper_crit, lower_crit)
    return(fig)
  })
  
  output$posterior3DPlot <- renderPlotly({
    
    # alpha <- input$alpha
    # beta <- input$beta
    # theta_u <- input$pi_u
    # theta_l <- input$pi_l
    # p_u <- input$p_u
    # p_l <- input$p_l
    # 
    # ## Compute critical values
    # ## vars for critical value calculation
    # upper_cdf  <- posteriorCDF(1:n, theta_u, alpha, beta) # P(theta > theta_u|n) = 1 - probbeta(theta, alpha+k, beta+n-k)
    # lower_cdf  <- 1-posteriorCDF(1:n, theta_l, alpha, beta) # P(theta < theta_l|n) = probbeta(theta_l, alpha+k, beta+n-k)
    # upper_crit <- min(which(upper_cdf >= p_u))
    # lower_crit <- max(which(lower_cdf >= p_l))

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
    
    s <- event_data("plotly_click", source="posterior3D")
    theta <- s$x
    thetas <- seq(0, 1, length.out=m)
    
    i <- match(theta, thetas)
    dens <- dens[, i]
    dens <- dens/sum(dens)
    dens <- cumsum(dens)
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
    s <- event_data("plotly_click", source="posterior3D")
    k <- s$y
    thetas <- seq(0, 1, length.out=m)
    
    dens <- dens[k, ]
    dens <- dens/sum(dens)
    dens <- cumsum(dens)
    
    print(dens[10])
    print(pbinom(.1, 1,1))
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
             xaxis = list(title = 'theta_0',
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
  

  # output$powerCurvePlot <- renderPlot({
  # 
  #   n <- input$n
  #   alpha <- input$alpha
  #   beta  <- input$beta
  #   theta_u <- input$pi_u
  #   theta_l <- input$pi_l
  #   p_u <- input$p_u
  #   p_l <- input$p_l
  # 
  #   upper_cdf <- 1-posteriorCDF(1:n, theta_u, alpha, beta) # P(theta > theta_u|n) = 1 - probbeta(theta, alpha+k, beta+n-k)
  #   lower_cdf <- posteriorCDF(1:n, theta_l, alpha, beta) # P(theta < theta_l|n) = probbeta(theta_l, alpha+k, beta+n-k)
  #   upper_crit <- max(which(upper_cdf >= p_u))
  #   lower_crit <- min(which(lower_cdf >= p_l))
  # 
  #   thetas <- seq(0, 1, length.out=100)
  #   upper_power <- posterior2(thetas, n, upper_crit, alpha, beta)
  #   upper_power <- upper_power/length(thetas)
  #   upper_power <- 1-cumsum(upper_power)
  # 
  #   lower_power <- posterior2(thetas, n, lower_crit, alpha, beta)
  #   lower_power <- lower_power/length(thetas)
  #   lower_power <- cumsum(lower_power)
  # 
  #   overall_power <- lower_power- (1-upper_power)
  #   overall_power <- overall_power/length(thetas)
  #   overall_power <- overall_power/sum(overall_power)
  #   fig <- plotPowerCurves(thetas,
  #                          upper_power, lower_power, overall_power,
  #                          upper_crit, lower_crit)
  #   return(fig)
  # })

  output$powerCurvePlot2 <- renderPlotly({

    # eventdata <- event_data("plotly_hover", source = "cdfPlot")
    # validate(need(!is.null(eventdata), "Hover over the cdf chart to populate this power plot"))

    n <- input$n
    alpha <- input$alpha
    beta  <- input$beta
    theta_u <- input$pi_u
    theta_l <- input$pi_l
    p_u <- input$p_u
    # p_u <- subset(eventdata, curveNumber==0)$y
    p_l <- input$p_l
    # p_l <- subset(eventdata, curveNumber==1)$y

    upper_cdf <- 1-posteriorCDF(1:n, theta_u, alpha, beta) # P(theta > theta_u|n) = 1 - probbeta(theta, alpha+k, beta+n-k)
    lower_cdf <- posteriorCDF(1:n, theta_l, alpha, beta) # P(theta < theta_l|n) = probbeta(theta_l, alpha+k, beta+n-k)
    upper_crit <- max(which(upper_cdf >= p_u))
    lower_crit <- min(which(lower_cdf >= p_l))

    thetas <- seq(0, 1, length.out=100)
    upper_power <- posterior2(thetas, n, upper_crit, alpha, beta)
    upper_power <- upper_power/length(thetas)
    upper_power <- cumsum(upper_power)

    lower_power <- posterior2(thetas, n, lower_crit, alpha, beta)
    lower_power <- lower_power/length(thetas)
    lower_power <- 1-cumsum(lower_power)

    overall_power <- (1-lower_power) - upper_power
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
