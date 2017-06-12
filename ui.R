
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(plotly)

shinyUI(fluidPage(
  
  # includeCSS(file.path("www/bootstrap3_3_7.css")),
  # includeCSS(file.path("www/normalize.css")),

  # includeCSS(file.path("www/ion.rangeSlider.adepro.css")),
  # includeCSS(file.path("www/ion.rangeSlider.skinModern.css")),
  # includeCSS(file.path("www/ion.rangeSlider.skinHTML5.css")),
  # includeCSS(file.path("www/ion.rangeSlider.skinNice.css")),
  # includeCSS(file.path("www/ion.rangeSlider.skinSimple.css")),
  # includeCSS(file.path("www/ion.rangeSlider.skinFlat.css")),
  # includeCSS(file.path("www/ion.rangeSlider.css")),

  # includeCSS(file.path("www/styles.css")),
  
  withMathJax(),
  
  # Application title
  titlePanel("Sample Size Calculation"),
  
  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      h3("Likelihood Parameters"),
      h3(withMathJax("$$X\\sim Bin(n, \\pi), with$$")),
      # fluidRow(
      #   column(4,
      #          numericInput("seed", withMathJax("$$\\textbf{Seed}\\ \\textbf{for}\\ \\textbf{sampling}$$"), 
      #                       value = 42))),
      fluidRow(
        column(4,
               numericInput("n", withMathJax("$$\\textbf{Number}\\ \\textbf{of}\\ \\textbf{samples}\\ \\textbf{n}$$"), 
                            min = 1, value = 10))),
      fluidRow(
        column(4,
          radioButtons("radioSample", "Sample by:",
                       choices=c("Success Probability"="prob",
                                 "Number of Successes"="successes")))),
      fluidRow(
        column(4,
          uiOutput("samplingInput"))),
      h3("Prior Parameters"),
      h3(withMathJax("$$p(\\pi)=Beta(\\pi|\\alpha, \\beta), with$$")),
      fluidRow(
        column(3, 
               numericInput("alpha", withMathJax("$$\\alpha$$"), 
                            min = 0, step=.05, value = 1)),
        column(3, 
               numericInput("beta", withMathJax("$$\\beta$$"), 
                            min = 0, step=.05, value = 1))
      )
      # h3("Hypothesis Parameters"),
      # h3(withMathJax("$$P(\\pi\\geq\\pi_{u})>p_{u}:\\ Go$$")),
      # fluidRow(
      #   column(4,
      #          numericInput("pi_u", withMathJax("$$\\pi_{u}$$"),
      #                       min = 0, max = 1, step=.01, value = .5))
        # column(4,
        #        numericInput("p_u", withMathJax("$$p_{u}$$"),
        #        min = 0, max = 1, step=.01, value = .5))
      # ),
      # h3(withMathJax("$$P(\\pi<\\pi_{l})>p_{l}:\\ No\\ Go$$"))
      # fluidRow(
      #   column(4,
      #          numericInput("pi_l", withMathJax("$$\\pi_{l}$$"),
      #          min = 0, max = 1, step=.01, value = .5))
        # column(4,
        #        numericInput("p_l", withMathJax("$$p_{l}$$"),
        #        min = 0, max = 1, step=.01, value = .5))
      # )
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      fluidRow(column(12, plotlyOutput("triPlot"))),
      fluidRow(column(6, plotlyOutput("posteriorCDFPlot2")),
               h3("Hypothesis Parameters"),
               column(3, h4(withMathJax("$$P(\\pi\\geq\\pi_{u})>p_{u}:\\ Go$$")),
                         sliderInput("pi_u", withMathJax("$$\\pi_{u}$$"),
                                     min = 0, max = 1, step=.01, value = .5,
                                     animate=T),
                         sliderInput("p_u", withMathJax("$$p_{u}$$"),
                                     min = 0, max = 1, step=.01, value = .5,
                                     animate=T)),
               column(3, h4(withMathJax("$$P(\\pi<\\pi_{l})>p_{l}:\\ No\\ Go$$")),
                         sliderInput("pi_l", withMathJax("$$\\pi_{l}$$"),
                                     min = 0, max = 1, step=.01, value = .5,
                                     animate=T),
                         sliderInput("p_l", withMathJax("$$p_{l}$$"),
                                     min = 0, max = 1, step=.01, value = .5,
                                     animate=T))),
      fluidRow(column(6, plotlyOutput("powerCurvePlot2"))),
      # fluidRow(column(6, plotOutput("posteriorCDFPlot")),
      #          column(6, plotOutput("powerCurvePlot"))
      # ),
      h2("Parametrization"),
      wellPanel(
        fluidRow(
          column(4, 
                 h3("Prior"), 
                 withMathJax(uiOutput("priorDistFormula"))),
          column(4, 
                 h3("Likelihood"), 
                 withMathJax(uiOutput("likelihoodFormula"))),
          column(4, 
                 h3("Posterior"), 
                 withMathJax(uiOutput("posteriorFormula")))
        )
      ),
      h2("Point Estimates"),
      wellPanel(
        fluidRow(
          column(4, 
                 h3("Prior"),
                 tableOutput("pointEst_Prior")),
          column(4, 
                 h3("Likelihood"),
                 tableOutput("pointEst_Likelihood")),
          column(4, 
                 h3("Posterior"),
                 tableOutput("pointEst_Posterior")))
      )
    )
  )
))
