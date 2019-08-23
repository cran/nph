##########################
# About tab UI

tabPanel(
  title = "About",
  id    = "aboutTab",
  value = "aboutTab",
  name  = "aboutTab",
  class = "home",
  icon  = icon("info-circle"),
  includeMarkdown(file.path("text", "about.md"))
  # includeMarkdown("inst/shiny/text/about.md")
)