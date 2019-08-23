nph.env<-new.env()

library(ggplot2)
library(dplyr)
library(formatR)

#library(nph)
server <- shinyServer(function(input, output, session) {
  options(stringsAsFactors = F)
  # # include logic for each tab
  source(file.path("server", "srv-tab-home.R"),    local = TRUE)$value
  source(file.path("server", "srv-tab-simulate.R"),   local = TRUE)$value
  source(file.path("server", "srv-tab-power.R"),local = TRUE)$value
  source(file.path("server", "srv-tab-condpower.R"),local = TRUE)$value
  source(file.path("server", "helper_functions.R"),local = TRUE)$value
  source(file.path("server", "paket_nph_1-9.R"),local = TRUE)$value
  source(file.path("server", "additional_functions.R"),local = TRUE)$value
  source(file.path("server", "continuous_functions.R"),local = TRUE)$value
  


})


