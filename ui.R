
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(shinythemes)
library(plotly)


shinyUI(fluidPage(
  # themeSelector(), ## dynamic theme
  theme = shinytheme("superhero"), ## fixed theme
  withMathJax(),

  # list(tags$head(
  #   HTML(
  #     '<link rel="icon", href="Logo_Cross_Screen_RGB.png", type="image/png" />'
  #   )
  # )),
  #
  # logo <- img(src="Logo_Cross_Screen_RGB.png",
  #             height = 80,
  #             width = 80,
  #             style = "margin:10px 10px;
  #                      float:right;"),


  # ## Title without logo
  titlePanel(title = "Sample Size Calculation",
             windowTitle = "Sample Size Calculation"),

  # includeCSS(file.path("www/bootstrap3_3_7.css")),
  # includeCSS(file.path("www/normalize.css")),

  # includeCSS(file.path("www/ion.rangeSlider.adepro.css")),
  # includeCSS(file.path("www/ion.rangeSlider.skinModern.css")),
  # includeCSS(file.path("www/ion.rangeSlider.skinHTML5.css")),
  # includeCSS(file.path("www/ion.rangeSlider.skinNice.css")),
  # includeCSS(file.path("www/ion.rangeSlider.skinSimple.css")),
  # includeCSS(file.path("www/ion.rangeSlider.skinFlat.css")),
  # includeCSS(file.path("www/ion.rangeSlider.css")),

  includeCSS(file.path("www/styles.css")),

  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      h3("Likelihood Parameters"),
      h3(withMathJax("$$X\\sim Bin(n, \\Theta), with$$")),
      # fluidRow(
      #   column(4,
      #          numericInput("seed", withMathJax("$$\\textbf{Seed}\\ \\textbf{for}\\ \\textbf{sampling}$$"),
      #                       value = 42))),
      fluidRow(
        column(4,
               numericInput("n", withMathJax("$$\\textbf{Number}\\ \\textbf{of}\\ \\textbf{samples}\\ \\textbf{n}$$"),
                            min = 1, value = 10))),
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
      tabsetPanel(type = "tabs",
                  tabPanel("Posterior Estimation",
                           br(),
                           wellPanel(
                             fluidRow(column(8, plotlyOutput("triPlot")
                                                # plotlyOutput("trialPlot")
                                             ),
                                      fluidRow(
                                        column(4,
                                              radioButtons("radioSample", "Sample by:",
                                                           choices=c("Success Probability"="prob",
                                                                     "Number of Successes"="successes")),
                                              uiOutput("samplingInput")))
                             )
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
                           )),
                  # tabPanel("Power Calculation",
                  #          br(),
                  #          wellPanel(
                  #            fluidRow(column(6, plotlyOutput("posteriorCDFPlot2")),
                  #                     # h3("Hypothesis Parameters"),
                  #                     column(3, h4(withMathJax("$$P(\\pi\\geq\\pi_{u})>p_{u}:\\ Go$$")),
                  #                            sliderInput("pi_u", withMathJax("$$\\pi_{u}$$"),
                  #                                        min = 0, max = 1, step=.01, value = .5,
                  #                                        animate=T),
                  #                            sliderInput("p_u", withMathJax("$$p_{u}$$"),
                  #                                        min = 0, max = 1, step=.01, value = .5,
                  #                                        animate=T)),
                  #                     column(3, h4(withMathJax("$$P(\\pi<\\pi_{l})>p_{l}:\\ No\\ Go$$")),
                  #                            sliderInput("pi_l", withMathJax("$$\\pi_{l}$$"),
                  #                                        min = 0, max = 1, step=.01, value = .5,
                  #                                        animate=T),
                  #                            sliderInput("p_l", withMathJax("$$p_{l}$$"),
                  #                                        min = 0, max = 1, step=.01, value = .5,
                  #                                        animate=T)))
                  #          ),
                  #          wellPanel(
                  #           fluidRow(column(6, plotlyOutput("powerCurvePlot2")))
                  #          )),
                  tabPanel("Exploration",
                           br(),
                           wellPanel(
                             fluidRow(column(12, plotlyOutput("posterior3DPlot")))
                           ),
                           wellPanel(
                             fluidRow(column(6, plotlyOutput("posterior2D_k")),
                                      column(6, plotlyOutput("posterior2D_theta")))
                           )),
                  tabPanel("Power Calculation",
                           # h3("Hypothesis Parameters"),
                           br(),
                           wellPanel(
                             fluidRow(
                               column(8, plotlyOutput("powerCurvePlot")),
                               column(4, 
                                      sliderInput("theta_0", withMathJax("$$\\theta_{0}$$"),
                                                  min = 0, max = 1, step=.01, value = .5,
                                                  animate=T),
                                      h4(withMathJax("$$P(\\theta\\geq\\theta_{0})>p_{u}:\\ Go$$")),
                                      sliderInput("p_go", withMathJax("$$p_{u}$$"),
                                                  min = 0, max = 1, step=.01, value = .5,
                                                  animate=T),
                                      h4(withMathJax("$$P(\\theta<\\theta_{0})>p_{l}:\\ No\\ Go$$")),
                                      sliderInput("p_nogo", withMathJax("$$p_{l}$$"),
                                                  min = 0, max = 1, step=.01, value = .5,
                                                  animate=T)
                                      )
                             )),
                           wellPanel(
                             fluidRow(
                               column(8, plotlyOutput("critValPlot"))
                             )
                           )
                           ),
                  tabPanel("About",
                           includeHTML("www/about.html")
                           )
                  )
    )
  )
))
