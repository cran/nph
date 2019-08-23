##########################
# Power tab server


# downloadHandler() takes two arguments, both functions.
# The content function is passed a filename as an argument, and
#   it should write out data to that filename.
output$downloadData <- downloadHandler(
  
  # This function returns a string which tells the client
  # browser what name to use when saving the file.
  filename = function() {
    "nph_dataexample.csv"
  },
  
  # This function should write data to a file given to it by
  # the argument 'file'.
  content = function(file) {
    set.seed(input$cp_seed)
    data_example = sample_fun(
      A = K5(),
      B = B5(),
      r0 = input$cp_r0,
      eventEnd = 100,
      lambdaRecr = cp_lambdaRekr(),
      lambdaCens = cp_lambdaZens(),
      maxRecrCalendarTime = cp_maxRekrKalenderZeit(),
      maxCalendar = cp_maxKalender()
    )
    print(data_example)
    # data_example = data_example[, c(1,3,5)]
    # names(data_example) = c("group", "time", "event")
    # Write to a file specified by the 'file' argument
    write.csv(data_example, file,
              row.names = FALSE)
  }
)


nph.env$outputData<-NULL

Dataset <- shiny::reactive({
  if(is.null(input$outputfile)){
    return(data.frame())
  }
  dat <- data.frame(do.call("read.csv",
                            list(input$outputfile$datapath)))
  nph.env$outputData <- dat
  return(nph.env$outputData)
  
})

output$head_table <- renderTable(head(Dataset()),
                                 caption = "Preview first 6 rows of the dataset",
                                 caption.placement = getOption("xtable.caption.placement", "top"),
                                 caption.width = getOption("xtable.caption.width", NULL))

output$choose_Vars<-shiny::renderUI({
  if(is.null(input$outputfile))
    return("Please upload a dataset first")
  if(identical(Dataset(),'')||identical(Dataset(),data.frame()))
    return("Please upload a dataset first")
  nph.env$outputData <- Dataset()
  nph.env$NUM     <- dim(nph.env$outputData)[2] #Num of all variable  ##
  nph.env$Class   <- sapply(apply(nph.env$outputData,2,unique),length)  ##
  nph.env$Colnames <- colnames(nph.env$outputData) ##
  out <- column(6,
                selectInput("group", label="Grouping variable", c("",nph.env$Colnames), selected = ifelse("group" %in% nph.env$Colnames, "group", "")),
                selectInput("event", label="Event variable",    c("",nph.env$Colnames), selected = ifelse("event" %in% nph.env$Colnames, "event", "")),
                selectInput("y",     label="Time variable",     c("",nph.env$Colnames), selected = ifelse("y" %in% nph.env$Colnames, "y", "")),
                selectInput("inclusion", label="Calendar inclusion time variable", c("",nph.env$Colnames), selected = ifelse("inclusion" %in% nph.env$Colnames, "inclusion", "")),
                selectInput("yCalendar", label="Calendar event time variable",  c("",nph.env$Colnames), selected = ifelse("yCalendar" %in% nph.env$Colnames, "yCalendar", "")),
		    selectInput("adminCens", label="Administrative interim censoring indicator",  c("",nph.env$Colnames), selected = ifelse("adminCens" %in% nph.env$Colnames, "adminCens", ""))
		)
  out
})

group_opts <- shiny::eventReactive(input$group, {
  if(is.null(input$outputfile))
    return()
  if(identical(Dataset(),'')||identical(Dataset(),data.frame()))
    return(NULL)
  if(is.null(input$group)){
    return(NULL)
  }
  dat = Dataset()
  unique(dat[input$group])
})

output$choose_groups <- shiny::renderUI({
  if(is.null(input$outputfile))
    return()
  if(identical(Dataset(),'')||identical(Dataset(),data.frame()))
    return(NULL)
  dat = Dataset()
  nph.env$outputData <- dat
  nph.env$NUM     <- dim(nph.env$outputData)[2] #Num of all variable  ##
  nph.env$Class   <- sapply(apply(nph.env$outputData,2,unique),length)  ##
  nph.env$Colnames <- colnames(nph.env$outputData) ##
  # if(length(input$group) == 0){
  #   return(NULL)
  # }
  out <- column(6,
                selectInput("data_treatment", label="Treatment value in group variable", choices = c("",group_opts()), selected = ifelse(1 %in% group_opts(), "1", "")))
                # selectInput("data_control",   label="Control",   choices = c("",group_opts()), selected = ifelse(0 %in% group_opts(), "0", "")))
  out
})


# Create ---------------------------------------------------------------------
# Initial definitions --------------------------------------------------------
cp_lambdaRekr <- reactive({
  switch(input$cp_lambdaRekr_unit,
         "day"   = input$cp_lambdaRekr,
         "month" = input$cp_lambdaRekr / (365.25/12),
         "year"  = input$cp_lambdaRekr / 365.25)
})
cp_lambdaZens <- reactive({
  switch(input$cp_lambdaZens_unit,
         "day"   = input$cp_lambdaZens,
         "month" = input$cp_lambdaZens / (365.25/12),
         "year"  = input$cp_lambdaZens / 365.25)
})
cp_maxRekrKalenderZeit <- reactive({
  switch(input$cp_maxRekrKalenderZeit_unit,
         "day"   = input$cp_maxRekrKalenderZeit,
         "month" = input$cp_maxRekrKalenderZeit * (365.25/12),
         "year"  = input$cp_maxRekrKalenderZeit * 365.25)
})
cp_maxKalender <- reactive({
  switch(input$cp_maxKalender_unit,
         "day"   = input$cp_maxKalender,
         "month" = input$cp_maxKalender * (365.25/12),
         "year"  = input$cp_maxKalender * 365.25)
})


# Simulation of a data set until interim analysis at 150 events
# dat <- eventReactive(input$cp_calculate, {
#   set.seed(input$s_seed)
#   sample_fun(
#     A = K5(),
#     B = B5(),
#     r0 = input$cp_r0,
#     eventEnd = input$cp_eventinterim,
#     lambdaRecr = cp_lambdaRekr(),
#     lambdaCens = cp_lambdaZens(),
#     maxRecrCalendarTime = cp_maxRekrKalenderZeit(),
#     maxCalendar = cp_maxKalender()
#   )
# })

# Calculate power --------------------------------------------------------
cp_pow <- eventReactive(input$cp_calculate, {
  set.seed(input$cp_seed)
  
  ## We do not simulate one dataset but use the one provided by the user
  # dat = sample_fun(
  #   A = K5(),
  #   B = B5(),
  #   r0 = input$cp_r0,
  #   eventEnd = input$cp_eventinterim,
  #   lambdaRecr = cp_lambdaRekr(),
  #   lambdaCens = cp_lambdaZens(),
  #   maxRecrCalendarTime = cp_maxRekrKalenderZeit(),
  #   maxCalendar = cp_maxKalender()
  # )
  dat = Dataset()
  group = 1*(dat[input$group] == input$data_treatment)
  event = dat[input$event] 
  y = dat[input$y] 
  inclusion = dat[input$inclusion] 
  yCalendar = dat[input$yCalendar] 
  #neu: adminCens
  adminCens = dat[input$adminCens]
  dat = data.frame(group = group, inclusion = inclusion, y = y, yCalendar = yCalendar, event = event, adminCens=adminCens)
  dat$cumEvents <- cumsum(dat$event)
  
  colMeans(do.call(rbind, lapply(1:input$cp_R, function(i){
    if (i %% 10 == 0) {
      cat(i, "\n")
      flush.console()
    }
    
    # Draw random data
    dat_cond <- sample_conditional_fun(
      dat = dat,
      A = K5(),
      B = B5(),
      r0 = input$cp_r0,
      eventEnd = input$cp_eventEnde,
      lambdaRecr = cp_lambdaRekr(),
      lambdaCens = cp_lambdaZens(),
      maxRecrCalendarTime = cp_maxRekrKalenderZeit(),
      maxCalendar = cp_maxKalender()
    )
    LRT <- logrank.maxtest(
      time = dat_cond$y,
      event = dat_cond$event,
      group = dat_cond$group,
      rho = c(0, 0, 1, 1),
      gamma = c(0, 1, 0, 1)
    )
    c(LRT$tests$p <= input$cp_alpha, LRT$pmult <= input$cp_alpha)
  })))
})


# Output results --------------------------------------------------------
output$cp_results <- renderTable({
  data.frame(
    Test = c("Maximum Logrank", rep("Weighted Logrank", 4)),
    rho = c("0011", 0, 0, 1, 1),
    gamma = c("0101", 0, 1, 0, 1),
    power = cp_pow())
},
caption = "<h4 style='text-align:center;color:black;'><b>Empirical conditional power</b></h4>", 
caption.placement = getOption("xtable.caption.placement", "top"), 
caption.width = getOption("xtable.caption.width", NULL),
align = 'lccr', 
width = "50%")

output$cp_interimplot <- renderPlot({
  
  
  dat = Dataset()
  group = 1*(dat[input$group] == input$data_treatment)
  event = dat[input$event] 
  y = dat[input$y] 
  inclusion = dat[input$inclusion] 
  yCalendar = dat[input$yCalendar] 
  dat = data.frame(group = group, inclusion = inclusion, y = y, yCalendar = yCalendar, event = event)
  dat$cumEvents <- cumsum(dat$event)
  
  intplot(dat, K5(), B5(), main = "Interim data")
})


###########---------------------------------------------------------------------

output$cp_downloadReport <- downloadHandler(
  filename = function() {
    paste('my-report-condpower', sep = '.', switch(
      input$cp_format, PDF = 'pdf', HTML = 'html', Word = 'docx'
    ))
  },
  
  content = function(file) {
    src <- normalizePath('nph-report-condpower.Rmd')
    
    # temporarily switch to the temp dir, in case you do not have write
    # permission to the current working directory
    owd <- setwd(tempdir())
    on.exit(setwd(owd))
    file.copy(src, 'nph-report-condpower.Rmd', overwrite = TRUE)
    # Set up parameters to pass to Rmd document
    params <- list(input = input,
                   cp_pow = cp_pow(),
                   intplot = intplot,
                   dat = Dataset(),
                   K5 = K5(),
                   B5 = B5())
    
    library(rmarkdown)
    out <- render('nph-report-condpower.Rmd', switch(
      input$format,
      PDF = pdf_document(), HTML = html_document(), Word = word_document()
    ),  params = params,  envir = new.env(parent = globalenv()))
    file.rename(out, file)
  }
)




# Download code ----------------------------------------------------------------
output$cp_downloadCode <- downloadHandler(
  filename = paste('nph-condpower-rcode.r'),
  content = function(file) {
    
    input1 = list()
    input1$times = 1
    input1$lambda = 2
    a = sprintf(
      "################################################################################
      ## nph package - Reproducible R code from shiny app.
      
      library('nph')
      
      # We first define the dataset to be used. 
      # For reproducibility, we have copied it here, but you could also read it from the file
      group = %s
      event = %s
      y = %s
      inclusion = %s
      yCalendar = %s
      
      dat = data.frame(group = group, inclusion = inclusion, y = y, yCalendar = yCalendar, event = event)
      dat$cumEvents <- cumsum(dat$event)

      # Define model parameters
      times = %s                              # Define timebreaks for intervals
      
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
      
      # Now plot model curves together with Kaplan-Meier curves
      intplot(dat, TRT, CRL, main = 'Interim data')
      
      set.seed(%s)
      
      ### HINT: To speed up the computation, you can also use mclapply for parallel computing instead of lapply
      results = lapply(1:%s, function(i){
      if (i %%%% 10 == 0) {  # Counter to display in the console every 10 replications
      cat(i, '\\n')
      flush.console()
      }
      
      # Draw random data
      dat_cond <- sample_conditional_fun(
      dat = dat,
      A = TRT,
      B = CRL,
      r0 = %s,
      eventEnd = %s,
      lambdaRecr = %s,                       
      lambdaCens = %s,                       
      maxRecrCalendarTime = %s,              
      maxCalendar = %s                       
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
      sprintf("c(%s)", paste0(1*(nph.env$outputData[input$group] == input$data_treatment), collapse = ", ")),   
      sprintf("c(%s)", paste0(nph.env$outputData[input$event], collapse = ", ")),  
      sprintf("c(%s)", paste0(nph.env$outputData[input$y], collapse = ", ")),   
      sprintf("c(%s)", paste0(nph.env$outputData[input$inclusion], collapse = ", ")),   
      sprintf("c(%s)", paste0(nph.env$outputData[input$yCalendar], collapse = ", ")),  
      
      sprintf("c(%s)", paste0(s_T_times(), collapse = ", ")),   # The time intervals
      sprintf("m2r(matrix(c(%s), nrow = %s, ncol = %s))", paste0(m2r(s_T_lambdaMat1()), collapse = ", "), input$s_T_subgr, input$s_T_timepoints),      
      sprintf("m2r(matrix(c(%s), nrow = %s, ncol = %s))", paste0(m2r(s_T_lambdaMat2()), collapse = ", "), input$s_T_subgr, input$s_T_timepoints),      
      sprintf("m2r(matrix(c(%s), nrow = %s, ncol = %s))", paste0(m2r(s_T_lambdaProgMat()), collapse = ", "), input$s_T_subgr, input$s_T_timepoints),      
      sprintf("c(%s)", paste0(s_T_p(), collapse = ", ")),       # The proportions in the subgroups for treatment arm
      sprintf("m2r(matrix(c(%s), nrow = %s, ncol = %s))", paste0(m2r(s_T_lambdaMat1()), collapse = ", "), input$s_T_subgr, input$s_T_timepoints),      
      sprintf("m2r(matrix(c(%s), nrow = %s, ncol = %s))", paste0(m2r(s_T_lambdaMat2()), collapse = ", "), input$s_T_subgr, input$s_T_timepoints),      
      sprintf("m2r(matrix(c(%s), nrow = %s, ncol = %s))", paste0(m2r(s_T_lambdaProgMat()), collapse = ", "), input$s_T_subgr, input$s_T_timepoints),      
      sprintf("c(%s)", paste0(s_C_p(), collapse = ", ")),
      input$cp_seed,
      input$cp_R,
      input$cp_r0,
      input$cp_eventEnde,
      cp_lambdaRekr(),
      cp_lambdaZens(),
      cp_maxRekrKalenderZeit(),
      cp_maxKalender(),
      input$p_alpha, input$p_alpha
    )
    
    owd <- setwd(tempdir())
    on.exit(setwd(owd))
    
    cat(a)
    capture.output(cat(a), file = paste0("test.R"))
    tidy_file(file = "test.R", width.cutoff = 80)   # this is to format nicely the R code https://yihui.name/formatr/
    # styler::style_file(path = "test.R", scope = "indention")  # this is to format nicely the R code
    
    out <- paste0("test.R")
    file.rename(out, file)
  }
    )


# Download dataset at interim analysis -----------------------------------------

output$cp_downloadData <- downloadHandler(
  filename = function() {
    paste('interim-data', sep = '.', switch(
      input$cp_dataformat, csv = 'csv', excel = 'xls', txt = 'txt'))
  },
  
  content = function(file) {
    # Write to a file specified by the 'file' argument
    savefunc = switch(
      input$cp_dataformat, csv = write.csv, excel = write.xlsx, txt = write.table)
    savefunc(Dataset(), file,
              row.names = FALSE)
  }
)
