##########################
# Power tab server

# Create ---------------------------------------------------------------------
# Initial definitions --------------------------------------------------------
p_lambdaRekr <- reactive({
  switch(input$p_lambdaRekr_unit,
         "day"   = input$p_lambdaRekr,
         "month" = input$p_lambdaRekr / (365.25/12),
         "year"  = input$p_lambdaRekr / 365.25)
})
p_lambdaZens <- reactive({
  switch(input$p_lambdaZens_unit,
         "day"   = input$p_lambdaZens,
         "month" = input$p_lambdaZens / (365.25/12),
         "year"  = input$p_lambdaZens / 365.25)
})
p_maxRekrKalenderZeit <- reactive({
  switch(input$p_maxRekrKalenderZeit_unit,
         "day"   = input$p_maxRekrKalenderZeit,
         "month" = input$p_maxRekrKalenderZeit * (365.25/12),
         "year"  = input$p_maxRekrKalenderZeit * 365.25)
})
p_maxKalender <- reactive({
  switch(input$p_maxKalender_unit,
         "day"   = input$p_maxKalender,
         "month" = input$p_maxKalender * (365.25/12),
         "year"  = input$p_maxKalender * 365.25)
})


# Calculate power --------------------------------------------------------
p_pow <- eventReactive(input$p_calculate, {
  set.seed(input$p_seed)
    colMeans(do.call(rbind, lapply(1:input$p_R, function(i){
      if (i %% 10 == 0) {
        cat(i, "\n")
        flush.console()
      }
      
      # Draw random data
      dat <- sample_fun(
        A = K5(),
        B = B5(),
        r0 = input$p_r0,
        eventEnd = input$p_eventEnde,
        lambdaRecr = p_lambdaRekr(),
        lambdaCens = p_lambdaZens(),
        maxRecrCalendarTime = p_maxRekrKalenderZeit(),
        maxCalendar = p_maxKalender()
      )
      LRT <- logrank.maxtest(
        time = dat$y,
        event = dat$event,
        group = dat$group,
        rho = c(0, 0, 1, 1),
        gamma = c(0, 1, 0, 1)
      )
      c(LRT$tests$p <= input$p_alpha, LRT$pmult <= input$p_alpha)
  })))
})


# Output results --------------------------------------------------------
output$p_results <- renderTable({
  data.frame(
    Test = c("Maximum Logrank", rep("Weighted Logrank", 4)),
    rho = c("0011", 0, 0, 1, 1),
    gamma = c("0101", 0, 1, 0, 1),
    power = p_pow())
}, 
caption = "<h4 style='text-align:center;color:black;'><b>Empirical power</b></h4>", 
caption.placement = getOption("xtable.caption.placement", "top"), 
caption.width = getOption("xtable.caption.width", NULL),
align = 'lccr', 
width = "50%")




###########---------------------------------------------------------------------

output$downloadReport <- downloadHandler(
  filename = function() {
    paste('my-report', sep = '.', switch(
      input$format, PDF = 'pdf', HTML = 'html', Word = 'docx'
    ))
  },
  
  content = function(file) {
    src <- normalizePath('nph-report.Rmd')
    
    # temporarily switch to the temp dir, in case you do not have write
    # permission to the current working directory
    owd <- setwd(tempdir())
    on.exit(setwd(owd))
    file.copy(src, 'nph-report.Rmd', overwrite = TRUE)
    # Set up parameters to pass to Rmd document
    params <- list(input = input,
                   p_pow = p_pow())
    
    library(rmarkdown)
    out <- render('nph-report.Rmd', switch(
      input$format,
      PDF = pdf_document(), HTML = html_document(), Word = word_document()
    ),  params = params,  envir = new.env(parent = globalenv()))
    file.rename(out, file)
  }
)


# Download code ----------------------------------------------------------------
output$p_downloadCode <- downloadHandler(
  filename = paste('nph-power-rcode.r'),
  content = function(file) {
    
    input1 = list()
    input1$times = 1
    input1$lambda = 2
    a = sprintf(
      "################################################################################
      ## nph package - Reproducible R code from shiny app.
      
      library('nph')
      
      times = %s                          # Define timebreaks for intervals
      
      # Create piecewise constant hazard objects for Treatment and Control arms
      # The lambda matrices are given in terms of the median survival times,
      # Therefore we use the nph::m2r function to translate from median to rate
      TRT <- pop_pchaz(
      Tint = times,
      lambdaMat1 = %s,   
      lambdaMat2 = %s,   
      lambdaProgMat = %s, 
      p = %s,
      timezero = FALSE, 
      discrete_approximation = TRUE
      )
      
      CRL <- pop_pchaz(
      Tint = times,
      lambdaMat1 = %s, 
      lambdaMat2 = %s,
      lambdaProgMat = %s, 
      p = %s,
      timezero = TRUE, 
      discrete_approximation = TRUE
      )
      
      print(TRT)
      print(CRL)
      
      # Create Survival, Hazard and hazard ratio plots
      plot_shhr(A = TRT, B = CRL)
      
      set.seed(%s)

      ### HINT: To speed up the computation, you can also use mclapply for parallel computing instead of lapply
      results = lapply(1:%s, function(i){
      if (i %%%% 10 == 0) {  # Counter to display in the console every 10 replications
      cat(i, '\\n')
      flush.console()
      }
      
      # Create one dataset
      dat <- sample_fun(
      A = TRT,
      B = CRL,
      r0 = %s,
      eventEnd = %s,
      lambdaRecr = %s,                       # Rate per day for recruiting patients, assuming recruitung follows a Poisson process
      lambdaCens = %s,                       # Rate per day for random censoring, assuming censoring times are exponential
      maxRecrCalendarTime = %s,              # Maximal duration of recruitment in days
      maxCalendar = %s                       # Maximal total study duration in days, after which the study stops
      )

      # Perform weighted log.rank tests
      LRT <- logrank.maxtest(
      time  = dat$y,
      event = dat$event,
      group = dat$group,
      rho   = c(0, 0, 1, 1),
      gamma = c(0, 1, 0, 1)
      )
      # Output the decision for each test
      c(LRT$tests$p <= %s, LRT$pmult <= %s)
      })
      
      power_table = power_table = data.frame(
      Test = c('Maximum Logrank', rep('Weighted Logrank', 4)),
      rho = c('0011', 0, 0, 1, 1),
      gamma = c('0101', 0, 1, 0, 1),
      power = colMeans(do.call(rbind, results)))
      
      print(power_table)
      ", 
    sprintf("c(%s)", paste0(s_T_times(), collapse = ", ")),   # The time intervals
    sprintf("m2r(matrix(c(%s), nrow = %s, ncol = %s))", paste0(m2r(s_T_lambdaMat1()), collapse = ", "), input$s_T_subgr, input$s_T_timepoints),      
    sprintf("m2r(matrix(c(%s), nrow = %s, ncol = %s))", paste0(m2r(s_T_lambdaMat2()), collapse = ", "), input$s_T_subgr, input$s_T_timepoints),      
    sprintf("m2r(matrix(c(%s), nrow = %s, ncol = %s))", paste0(m2r(s_T_lambdaProgMat()), collapse = ", "), input$s_T_subgr, input$s_T_timepoints),      
    sprintf("c(%s)", paste0(s_T_p(), collapse = ", ")),       # The proportions in the subgroups for treatment arm
    sprintf("m2r(matrix(c(%s), nrow = %s, ncol = %s))", paste0(m2r(s_T_lambdaMat1()), collapse = ", "), input$s_T_subgr, input$s_T_timepoints),      
    sprintf("m2r(matrix(c(%s), nrow = %s, ncol = %s))", paste0(m2r(s_T_lambdaMat2()), collapse = ", "), input$s_T_subgr, input$s_T_timepoints),      
    sprintf("m2r(matrix(c(%s), nrow = %s, ncol = %s))", paste0(m2r(s_T_lambdaProgMat()), collapse = ", "), input$s_T_subgr, input$s_T_timepoints),      
    sprintf("c(%s)", paste0(s_C_p(), collapse = ", ")),
    input$p_seed,
    input$p_R,
    input$p_r0,
    input$p_eventEnde,
    p_lambdaRekr(),
    p_lambdaZens(),
    p_maxRekrKalenderZeit(),
    p_maxKalender(),
    input$p_alpha, input$p_alpha
    )
    
    owd <- setwd(tempdir())
    on.exit(setwd(owd))
    capture.output(cat(a), file = paste0("test.R"))
    # tidy_file(file = "test.R", width.cutoff = 80)   # this is to format nicely the R code https://yihui.name/formatr/
    styler::style_file(path = "test.R", scope = "indention")  # this is to format nicely the R code
    
    out <- paste0("test.R")
    file.rename(out, file)
  }
)
