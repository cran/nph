##########################
# Simulate tab UI

tabPanel(
  title = "Conditional Power",
  id    = "condpowerTab",
  value = "condpowerTab",
  name  = "condpowerTab",
  class = "fade in",
  icon  = icon("line-chart"),
  h1("Calculate conditional power", style = "color:#f8c954"),
  #shiny::wellPanel(
   # downloadButton('cp_downloadCode', "Download R code", class = "start-button"),
  #  h4("Download a reproducible R script for the code to generate the results of this tab")),
  wellPanel(
    h2("1. Upload your data"),
    fluidRow(column(5,fileInput("outputfile",label="Select data file in csv")),
             # h5("Preview first 6 rows of the dataset"),
             fluidRow(column(6,tableOutput("head_table")))),
    fluidRow(column(12, "Note: The data should have an observation per row and four columns, one indicating the group (treatment or control), one indicating the time of the event, the event indicator, and the calendar time of the event.\n")),
    fluidRow(column(3, downloadButton('downloadData', 'Download an example dataset')))),
  wellPanel(
    h2("2. Choose the corresponding variables"),
    fluidRow(column(12,
                    uiOutput("choose_Vars"),
                    uiOutput("choose_groups"))
    )),
  wellPanel(
    h2("3. Select design parameters"),
    fluidRow(
      column(6, sliderInput("cp_R", "Number of simulations", min = 10, max = 1e4, value = 10, step = 50))),
    fluidRow(
      column(6, sliderInput("cp_r0", "Allocation ratio", min = 0, max = 1, value = 0.5))
      # ,column(6, numericInput("cp_eventinterim", "Maximal number of events in First stage data", min = 10, max = 1e6, value = 150))
      ,column(6, numericInput("cp_eventEnde", "Maximal number of events", min = 10, max = 1e6, value = 450))
     ), fluidRow(
      column(4, selectInput("cp_lambdaRekr_unit", "Recruitment rate unit",
                             choices = c("day", "month", "year"), selected = "year"))
      ,column(4, numericInput("cp_lambdaRekr", "Recruitment rate (Poisson assumption)", min = 0, max = 1e6, value = 300))
     ), fluidRow(
      column(4, selectInput("cp_lambdaZens_unit", "Censoring rate unit", 
                              choices = c("day", "month", "year"), selected = "year"))
      ,column(4, numericInput("cp_lambdaZens", "Censoring rate (Exponential assumption)", min = 0, max = 1e6, value = 0.013))
     ), fluidRow(
      column(4, selectInput("cp_maxRekrKalenderZeit_unit", "Maximal duration of recruitment unit", 
                              choices = c("day", "month", "year"), selected = "year"))
      ,column(4, numericInput("cp_maxRekrKalenderZeit", "Maximal duration of recruitment", min = 1, max = 1e6, value = 3))
     ), fluidRow(
      column(4, selectInput("cp_maxKalender_unit", "Maximal study duration unit", 
                              choices = c("day", "month", "year"), selected = "year"))
      ,column(4, numericInput("cp_maxKalender", "Maximal study duration", min = 1, max = 1e6, value = 4))
     ), fluidRow(
      column(6, numericInput("cp_alpha",
                              "Significance level for hypothesis test of equal hazard functions in both groups", 
                              min = 0, max = 1, value = 0.05))
      ,column(6, numericInput("cp_seed",
                              "Seed for simulations", 
                              min = 0, max = 1e5, value = floor(runif(1, min = 0, max = 1e5))))
    )
    )
  ,wellPanel(
    actionButton("cp_calculate", "Calculate power", icon = icon("calculator")),
    conditionalPanel(
      condition = "input.cp_calculate > 0",
      withSpinner(fluidRow(column(width = 12,
                      plotOutput(width = "50%", "cp_interimplot"), 
                      div(style = "min-height: 20px;"),
                      tableOutput("cp_results"), 
                      align="center"))),
      fluidRow(column(width = 6,
             h3("Create report"),
          radioButtons('cp_format', 'Document format', c('PDF', 'HTML', 'Word'),
                       inline = TRUE),
          downloadButton('cp_downloadReport')),
      column(width = 6, h3("Download interim data"),
          radioButtons('cp_dataformat', 'Data format', c('csv', 'excel', "txt"),
                       inline = TRUE),
          downloadButton('cp_downloadData'))),
      div(style = "min-height: 20px;"))
  )
)
