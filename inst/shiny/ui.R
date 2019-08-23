# Load packages ----------------------------------------------------------------
library(markdown)
library(shiny)
# library(shinyBS)
library(shinycssloaders)
library(xlsx)
# Load data --------------------------------------------------------------------

# Define UI for application ----------------------------------------------------
ui <- shinyUI(
  bootstrapPage(title = '', theme = "merck.css",
                navbarPage(
                  title = tags$b("nph"),
                  windowTitle = "nph",
                  id = "mainNav",
                  inverse = TRUE,
                  fluid = FALSE,
                  collapsible = FALSE,
                  # theme = "bootstrap.min.css",
                  source(file.path("ui", "ui-tab-home.R"),     local = TRUE)$value,
                  source(file.path("ui", "ui-tab-simulate.R"), local = TRUE)$value,
                  source(file.path("ui", "ui-tab-power.R"), local = TRUE)$value,
                  source(file.path("ui", "ui-tab-condpower.R"), local = TRUE)$value,
                  source(file.path("ui", "ui-tab-about.R"),    local = TRUE)$value
                  )
  )
)

