# Home tab UI ##################################################################

tabPanel(
  title = "Home",
  id    = "homeTab",
  value = "homeTab",
  name  = "homeTab",
  class = "home",
  icon  = icon("home"),
  div("nph: Tools for planning and analysing surivival studies under non-proportional hazards", 
      class = "cover-heading"),
  div("In randomized controlled
  studies, survival functions are typically compared using the log-rank test, which
  is most powerful under the assumption of proportional hazards. However, when the
  proportionality assumption is violated, hypothesis tests and visual inspection of estimated
  survival functions may be misleading or hard to interpret.",
      class = "cover-text"),
  div(class="start-button-div",
      actionButton("toSimulate", "START", 
                   class="start-button")),
  div(
  # div("A collaboration by:"),
  div(#img(src="merck.png", align = "left", width="30%"),
      a(href="https://cemsiis.meduniwien.ac.at", target="_blank",
        # img(src="muw.png",   align = "left", width="30%")),
        img(src="cemsiis_en.png",   align = "left", width="60%")),
      class = "cover-footing-img"),
      class = "cover-footing"),
  div(class = "cover-footing-end")
)


