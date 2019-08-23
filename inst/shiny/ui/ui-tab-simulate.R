##########################
# Simulate tab UI

tabPanel(
  title = "Simulate",
  id    = "simulateTab",
  value = "simulateTab",
  name  = "simulateTab",
  class = "fade in",
  icon  = icon("table"),
  
  h1("Population parameters", style = "color:#f8c954"),
  #shiny::wellPanel(
   # downloadButton('s_downloadCode', "Download R code", class = "start-button"),
 #   h4("Download a reproducible R script for the code to #generate the results of this tab")),
  shiny::wellPanel(
    fluidRow(tags$h3("Time intervals"),
             column(1, "Number of time intervals"),
             column(4, sliderInput("s_T_timepoints", label = "", 
                                   value = 1, min = 1, max = 2, step = 1)),
             column(6, tags$b("Limits:"), tableOutput("s_T_time")),
             column(1, "")
    )),
  wellPanel(
    fluidRow(
      column(6, tags$h3("Treatment Group")),
      column(6, tags$h3("Control Group"))
    ),
    fluidRow(
      column(6, tags$h4("Specify parameters for Treatment Group")),
      column(6, tags$h4("Specify parameters for Control Group"))
    ),
    fluidRow(
      column(6,
             fluidRow(
               fluidRow(column(4,
                               sliderInput("s_T_subgr", label = "Number of subgroups", 
                                           value = 1, min = 1, max = 2, step = 1))),
               tags$b("Subgroup proportions:"), tableOutput("s_T_subgr_p"),
               tags$b("Subgroup labels:"), tableOutput("s_T_subgr_labels")),
             fluidRow(column(6,
                             fluidRow(
                               column(12,
                                      tags$b("lambdaMat1:"), 
                                      tableOutput("s_T_lambdaMat1")
                               )
                             ),
                             fluidRow(
                               column(12,
                                      tags$b("lambdaMat2:"), 
                                      tableOutput("s_T_lambdaMat2")
                               )
                             ),
                             fluidRow(
                               column(12,
                                      tags$b("lambdaProgMat:"),
                                      tableOutput("s_T_lambdaProgMat")
                               )
                             ))),
             fluidRow(column(12, checkboxInput("show_trt_diag", label = "Show diagram for treatment arm", value = FALSE))),
             conditionalPanel(condition = "input.show_trt_diag > 0",
                              fluidRow(column(12, plotOutput("trt_plot"))))),
      column(6,
             fluidRow(
               fluidRow(column(4,
                               sliderInput("s_C_subgr", label = "Number of subgroups", 
                                           step = 1, value = 1, min = 1, max = 2))),
               tags$b("Subgroup proportions:"), tableOutput("s_C_subgr_p"),
               tags$b("Subgroup labels:"), tableOutput("s_C_subgr_labels")),
             fluidRow(column(6,
                             fluidRow(
                               column(12,
                                      tags$b("lambdaMat1:"), 
                                      tableOutput("s_C_lambdaMat1")
                               )),
                             fluidRow(
                               column(12,
                                      tags$b("lambdaMat2:"), 
                                      tableOutput("s_C_lambdaMat2")
                               )),
                             fluidRow(
                               column(12,
                                      tags$b("lambdaProgMat:"),
                                      tableOutput("s_C_lambdaProgMat")
                               )
                               ))),
              fluidRow(column(12, checkboxInput("show_ctrl_diag", label = "Show diagram for control arm", value = FALSE))),
              conditionalPanel(condition = "input.show_ctrl_diag > 0",
                              fluidRow(column(12, plotOutput("ctrl_plot")))))
    )),
    wellPanel(
    fluidRow(
      fluidRow(column(5, ""),
               column(2, actionButton(inputId = "s_calculate", label = "Plot!", icon = icon("line-chart"))),
               column(5, "")),
      conditionalPanel(
        condition = "input.s_calculate > 0",
        withSpinner(fluidRow(column(width = 12,
                                    h3("Results"),
                                    plotOutput("panplot"),
                                    div(style = "min-height: 20px;"),
                                    fluidRow(column(width = 12, tableOutput("median_table"), align="center")))
        ))))
  )
)



