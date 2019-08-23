##########################
# Simulate tab UI

tabPanel(
  title = "Power",
  id    = "powerTab",
  value = "powerTab",
  name  = "powerTab",
  class = "fade in",
  icon  = icon("line-chart"),
  h1("Calculate power", style = "color:#f8c954"),
  #shiny::wellPanel(
  #  downloadButton('p_downloadCode', "Download R code", class = "start-button"),
  #  h4("Download a reproducible R script for the code to generate the results of this tab")),
  wellPanel(
    fluidRow(
      column(6, sliderInput("p_R", "Number of simulations", min = 10, max = 1e4, value = 10, step = 50))),
    fluidRow(
      column(6, sliderInput("p_r0", "Allocation ratio", min = 0, max = 1, value = 0.5))
      ,column(6, numericInput("p_eventEnde", "Maximal number of events", min = 10, max = 1e6, value = 450))
     ), fluidRow(
      column(4, selectInput("p_lambdaRekr_unit", "Recruitment rate unit",
                             choices = c("day", "month", "year"), selected = "year"))
      ,column(4, numericInput("p_lambdaRekr", "Recruitment rate (Poisson assumption)", min = 0, max = 1e6, value = 300))
     ), fluidRow(
      column(4, selectInput("p_lambdaZens_unit", "Censoring rate unit", 
                              choices = c("day", "month", "year"), selected = "year"))
      ,column(4, numericInput("p_lambdaZens", "Censoring rate (Exponential assumption)", min = 0, max = 1e6, value = 0.013))
     ), fluidRow(
      column(4, selectInput("p_maxRekrKalenderZeit_unit", "Maximal duration of recruitment unit", 
                              choices = c("day", "month", "year"), selected = "year"))
      ,column(4, numericInput("p_maxRekrKalenderZeit", "Maximal duration of recruitment", min = 1, max = 1e6, value = 3))
     ), fluidRow(
      column(4, selectInput("p_maxKalender_unit", "Maximal study duration unit", 
                              choices = c("day", "month", "year"), selected = "year"))
      ,column(4, numericInput("p_maxKalender", "Maximal study duration", min = 1, max = 1e6, value = 4))
     ), fluidRow(
      column(6, numericInput("p_alpha",
                              "Significance level for hypothesis test of equal hazard functions in both groups", 
                              min = 0, max = 1, value = 0.05))
      ,column(6, numericInput("p_seed",
                              "Seed for simulations", 
                              min = 0, max = 1e5, value = floor(runif(1, min = 0, max = 1e5))))
    )
    )
  ,wellPanel(
    actionButton("p_calculate", "Calculate power", icon = icon("calculator")),
    conditionalPanel(
      condition = "input.p_calculate > 0",
      fluidRow(column(width = 12,
                      withSpinner(tableOutput("p_results")), 
                                  align="center")),
      fluidRow(column(width = 6,
                      h3("Create report"),
                      radioButtons('format', 'Document format', c('PDF', 'HTML', 'Word'),
                                   inline = TRUE),
                      downloadButton('downloadReport'))))
  )
)
