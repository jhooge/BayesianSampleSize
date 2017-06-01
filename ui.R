
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(plotly)

shinyUI(fluidPage(
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
      # fluidRow(
      #   column(4,
      #          numericInput("pi0", withMathJax("$$\\pi_{0}, with\\ H_{0}: \\pi\\leq\\pi_{0}, H_{1}: \\pi>\\pi_{0}$$"), 
      #                       min = 0, max = 1, step=.01, value = .5))),
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
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      fluidRow(column(11, 
                      plotlyOutput("triPlot"),
                      plotlyOutput("posteriorProbPlot"))
               ),
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
