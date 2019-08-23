##########################
# Simulate tab server
# Create ---------------------------------------------------------------------
# Initial definitions --------------------------------------------------------

## Treatment group parameters --------------------------------------------------

# Number of subgroups ----------------------------------------------------------
observe({
  isolate({
    output$s_T_subgr_p <-renderTable({
      if (input$s_T_subgr<1) return("Number of subgroups should be > 0")
      if (input$s_T_subgr>=10) return("Number of subgroups should be < 10")
      num.inputs.col1 <- paste0("<input id='s_T_p", 1:input$s_T_subgr, 
                                "' class='shiny-bound-input' ",
                                "type='number' max = 1 min = 0 step = 0.1 value=", 
                                1/input$s_T_subgr, " />")
      out = as.data.frame(t(num.inputs.col1),
                    row.names = "Subgroup proportions")
      names(out) = paste0("S", 1:input$s_T_subgr)
      out
    },
    sanitize.text.function = function(x) x, 
    # rownames = TRUE, 
    width = '100%')
    
    output$s_T_subgr_labels <-renderTable({
      if (input$s_T_subgr<1) return("Number of subgroups should be > 0")
      if (input$s_T_subgr>=10) return("Number of subgroups should be < 10")
      num.inputs.col1 <- paste0("<div class='form-group shiny-input-container' style='width: 20%;'><input id='s_T_label", 
                                1:input$s_T_subgr, 
                                "' class='shiny-bound-input' ",
                                "type='text' value='Subgr", 
                                1:input$s_T_subgr, "' /></div>")
      out = as.data.frame(t(num.inputs.col1),
                          row.names = "Subgroup labels")
      names(out) = paste0("S", 1:input$s_T_subgr)
      out
    },
    sanitize.text.function = function(x) x, 
    # rownames = TRUE, 
    width = '100%')
  })
})


s_T_p <- reactive({
  unlist(lapply(1:input$s_T_subgr, function(i){
    input[[paste0("s_T_p", i)]]
  }))
})
s_T_labels <- reactive({
  unlist(lapply(1:input$s_T_subgr, function(i){
    input[[paste0("s_T_label", i)]]
  }))
})

# output$s_T_p <- renderText(s_T_p())

# Number of time intervals -----------------------------------------------------
observe({
  isolate({
    output$s_T_time <-renderTable({
      if (input$s_T_timepoints<1) return("Number of timepoints should be > 0")
      if (input$s_T_timepoints>=10) return("Number of timepoints should be < 10")
      num.inputs.col1 <- 
        c("<input id='s_T_time0' class='shiny-bound-input' type='number' max = 999999 min = 0 step = 0 value=0 disabled=true>",
          paste0("<input id='s_T_time", 1:input$s_T_timepoints, 
                 "' class='shiny-bound-input' ",
                 "type='number' max = 999999 min = 0 step = 1 value=", 
                 (1:input$s_T_timepoints)*3000, ">"))
      out = as.data.frame(t(num.inputs.col1),
                          row.names = "Subgroup proportions")
      names(out) = paste0("T", 0:input$s_T_timepoints)
      out
    }, sanitize.text.function = function(x) x, 
    # rownames = TRUE, 
    width = '100%')
  })
})


s_T_times <- reactive({
  unlist(lapply(0:input$s_T_timepoints, function(i){
    input[[paste0("s_T_time", i)]]
  }))
})

# output$s_T_times <- renderText(s_T_times())






### lambdaMat1 -----------------------------------------------------------------
s_T_grid_Mat1 <- reactive({
  grid_comb. =  expand.grid(1:input$s_T_timepoints, 1:input$s_T_subgr)
  paste0("T_Mat1_s", grid_comb.$Var1, "t", grid_comb.$Var2)
})
observe({
  isolate({
    output$s_T_lambdaMat1 <-renderTable({
      if(is.null(s_T_grid_Mat1())) return(NULL)
      num.inputs.col1 <- 
        paste0("<input id='s_", s_T_grid_Mat1(), 
               "' class='shiny-bound-input' ",
               " type='number' max = 999 min = 0 step = 1 value=", 
                 18, ">")
      out = as.data.frame(matrix(num.inputs.col1, nrow = input$s_T_subgr))
      
      t1=paste0("s_T_time", (1:input$s_T_timepoints)-1)
      t2=paste0("s_T_time", (1:input$s_T_timepoints))
      names(out) = unlist(lapply(1:input$s_T_timepoints, function(i)
        paste0(input[[t1[i]]], "-", input[[t2[i]]])))
      # names(out) = paste0("T", 1:input$s_T_timepoints)
      rownames(out) = paste0("S", 1:input$s_T_subgr)
      out
    }, sanitize.text.function = function(x) x, 
    rownames = TRUE,
    width = 'auto')
  })
})


s_T_lambdaMat1 <- reactive({
  mm = unlist(lapply(1:length(s_T_grid_Mat1()), function(i){
    input[[paste0("s_", s_T_grid_Mat1()[i])]]
  }))
  if(is.null(mm)) return(NULL)
  if(length(mm)%%input$s_T_subgr != 0) return(NULL) # Catch error when mm is not ready to be in the matrix
  matrix(fu(mm), nrow = input$s_T_subgr)
})

output$s_T_lambdaMat1_print <- renderTable({
  if(is.null(s_T_grid_Mat1())) return(NULL)
  s_T_lambdaMat1()
  }, digits = 6)

### lambdaMat2 -----------------------------------------------------------------
s_T_grid_Mat2 <- reactive({
  grid_comb. =  expand.grid(1:input$s_T_timepoints, 1:input$s_T_subgr)
  paste0("T_Mat2_s", grid_comb.$Var1, "t", grid_comb.$Var2)
})
observe({
  isolate({
    output$s_T_lambdaMat2 <-renderTable({
      if(is.null(s_T_grid_Mat2())) return(NULL)
      num.inputs.col1 <- 
        paste0("<input id='s_", s_T_grid_Mat2(), 
               "' class='shiny-bound-input' ",
               " type='number' max = 999 min = 0 step = 1 value=", 
               18, ">")
      out = as.data.frame(matrix(num.inputs.col1, nrow = input$s_T_subgr))
      
      t1=paste0("s_T_time", (1:input$s_T_timepoints)-1)
      t2=paste0("s_T_time", (1:input$s_T_timepoints))
      names(out) = unlist(lapply(1:input$s_T_timepoints, function(i)
        paste0(input[[t1[i]]], "-", input[[t2[i]]])))
      rownames(out) = paste0("S", 1:input$s_T_subgr)
      out
    }, sanitize.text.function = function(x) x, 
    rownames = TRUE,
    width = 'auto')
  })
})


s_T_lambdaMat2 <- reactive({
  mm = unlist(lapply(1:length(s_T_grid_Mat2()), function(i){
    input[[paste0("s_", s_T_grid_Mat2()[i])]]
  }))
  if(is.null(mm)) return(NULL)
  if(length(mm)%%input$s_T_subgr != 0) return(NULL) # Catch error when mm is not ready to be in the matrix
  matrix(fu(mm), nrow = input$s_T_subgr)
})

output$s_T_lambdaMat2_print <- renderTable({
  if(is.null(s_T_grid_Mat2())) return(NULL)
  s_T_lambdaMat2()
  }, digits = 6)


### lambdaProgMat --------------------------------------------------------------
s_T_grid_Prog <- reactive({
  grid_comb. =  expand.grid(1:input$s_T_timepoints, 1:input$s_T_subgr)
  paste0("T_Prog_s", grid_comb.$Var1, "t", grid_comb.$Var2)
})

observe({
  isolate({
    output$s_T_lambdaProgMat <-renderTable({
      if(is.null(s_T_grid_Prog())) return(NULL)
      num.inputs.col1 <- 
        paste0("<input id='s_", s_T_grid_Prog(), 
               "' class='shiny-bound-input' ",
               " type='number' max = 999 min = 0 step = 1 value=", 
               18, ">")
      out = as.data.frame(matrix(num.inputs.col1, nrow = input$s_T_subgr))
      
      t1=paste0("s_T_time", (1:input$s_T_timepoints)-1)
      t2=paste0("s_T_time", (1:input$s_T_timepoints))
      names(out) = unlist(lapply(1:input$s_T_timepoints, function(i)
        paste0(input[[t1[i]]], "-", input[[t2[i]]])))
      rownames(out) = paste0("S", 1:input$s_T_subgr)
      out
    }, sanitize.text.function = function(x) x, 
    rownames = TRUE,
    width = 'auto')
  })
})


s_T_lambdaProgMat <- reactive({
  mm = unlist(lapply(1:length(s_T_grid_Prog()), function(i){
    input[[paste0("s_", s_T_grid_Prog()[i])]]
  }))
  if(is.null(mm)) return(NULL)
  if(length(mm)%%input$s_T_subgr != 0) return(NULL) # Catch error when mm is not ready to be in the matrix
  matrix(fu(mm), nrow = input$s_T_subgr)
})

output$s_T_lambdaProgMat_print <- renderTable({
  if(is.null(s_T_grid_Prog())) return(NULL)
  s_T_lambdaProgMat()
}, digits = 6)




## Control group parameters --------------------------------------------------

# Number of subgroups ----------------------------------------------------------
observe({
  isolate({
    output$s_C_subgr_p <-renderTable({
      if (input$s_C_subgr<1) return("Number of subgroups should be > 0")
      if (input$s_C_subgr>=10) return("Number of subgroups should be < 10")
      num.inputs.col1 <- paste0("<input id='s_C_p", 1:input$s_C_subgr, 
                                "' class='shiny-bound-input' ",
                                "type='number' max = 1 min = 0 step = 0.1 value=", 
                                1/input$s_C_subgr, " />")
      out = as.data.frame(t(num.inputs.col1),
                          row.names = "Subgroup proportions")
      names(out) = paste0("S", 1:input$s_C_subgr)
      out
    },
    sanitize.text.function = function(x) x, 
    # rownames = TRUE, 
    width = '100%')
    
    output$s_C_subgr_labels <-renderTable({
      if (input$s_C_subgr<1) return("Number of subgroups should be > 0")
      if (input$s_C_subgr>=10) return("Number of subgroups should be < 10")
      num.inputs.col1 <- paste0("<input id='s_C_label", 1:input$s_C_subgr, 
                                "' class='shiny-bound-input' ",
                                "type='text' value='Subgr", 
                                1:input$s_C_subgr, "' />")
      out = as.data.frame(t(num.inputs.col1),
                          row.names = "Subgroup labels")
      names(out) = paste0("S", 1:input$s_C_subgr)
      out
    },
    sanitize.text.function = function(x) x, 
    width = '100%')
  })
})


s_C_p <- reactive({
  unlist(lapply(1:input$s_C_subgr, function(i){
    input[[paste0("s_C_p", i)]]
  }))
})
s_C_labels <- reactive({
  unlist(lapply(1:input$s_T_subgr, function(i){
    input[[paste0("s_C_label", i)]]
  }))
})


### lambdaMat1 -----------------------------------------------------------------
s_C_grid_Mat1 <- reactive({
  grid_comb. =  expand.grid(1:input$s_T_timepoints, 1:input$s_C_subgr)
  paste0("C_Mat1_s", grid_comb.$Var1, "t", grid_comb.$Var2)
})
observe({
  isolate({
    output$s_C_lambdaMat1 <-renderTable({
      if(is.null(s_C_grid_Mat1())) return(NULL)
      num.inputs.col1 <- 
        paste0("<input id='s_", s_C_grid_Mat1(), 
               "' class='shiny-bound-input' ",
               " type='number' max = 999 min = 0 step = 1 value=", 
               18, ">")
      out = as.data.frame(matrix(num.inputs.col1, nrow = input$s_C_subgr))
      
      t1=paste0("s_T_time", (1:input$s_T_timepoints)-1)
      t2=paste0("s_T_time", (1:input$s_T_timepoints))
      names(out) = unlist(lapply(1:input$s_T_timepoints, function(i)
        paste0(input[[t1[i]]], "-", input[[t2[i]]])))
      # names(out) = paste0("T", 1:input$s_T_timepoints)
      rownames(out) = paste0("S", 1:input$s_C_subgr)
      out
    }, sanitize.text.function = function(x) x, 
    rownames = TRUE,
    width = 'auto')
  })
})


s_C_lambdaMat1 <- reactive({
  input$s_C_lM1_s1t1
  mm = unlist(lapply(1:length(s_C_grid_Mat1()), function(i){
    input[[paste0("s_", s_C_grid_Mat1()[i])]]
  }))
  if(is.null(mm)) return(NULL)
  if(length(mm)%%input$s_C_subgr != 0) return(NULL) # Catch error when mm is not ready to be in the matrix
  matrix(fu(mm), nrow = input$s_C_subgr)
})

output$s_C_lambdaMat1_print <- renderTable({
  if(is.null(s_C_grid_comb())) return(NULL)
  s_C_lambdaMat1()
})

### lambdaMat2 -----------------------------------------------------------------
s_C_grid_Mat2 <- reactive({
  grid_comb. =  expand.grid(1:input$s_T_timepoints, 1:input$s_C_subgr)
  paste0("C_Mat2_s", grid_comb.$Var1, "t", grid_comb.$Var2)
})
observe({
  isolate({
    output$s_C_lambdaMat2 <-renderTable({
      if(is.null(s_C_grid_Mat2())) return(NULL)
      num.inputs.col1 <- 
        paste0("<input id='s_", s_C_grid_Mat2(), 
               "' class='shiny-bound-input' ",
               " type='number' max = 999 min = 0 step = 1 value=", 
               18, ">")
      out = as.data.frame(matrix(num.inputs.col1, nrow = input$s_C_subgr))
      
      t1=paste0("s_T_time", (1:input$s_T_timepoints)-1)
      t2=paste0("s_T_time", (1:input$s_T_timepoints))
      names(out) = unlist(lapply(1:input$s_T_timepoints, function(i)
        paste0(input[[t1[i]]], "-", input[[t2[i]]])))
      # names(out) = paste0("T", 1:input$s_T_timepoints)
      rownames(out) = paste0("S", 1:input$s_C_subgr)
      out
    }, sanitize.text.function = function(x) x, 
    rownames = TRUE,
    width = 'auto')
  })
})


s_C_lambdaMat2 <- reactive({
  mm = unlist(lapply(1:length(s_C_grid_Mat2()), function(i){
    input[[paste0("s_", s_C_grid_Mat2()[i])]]
  }))
  if(is.null(mm)) return(NULL)
  if(length(mm)%%input$s_C_subgr != 0) return(NULL) # Catch error when mm is not ready to be in the matrix
  matrix(fu(mm), nrow = input$s_C_subgr)
})

output$s_C_lambdaMat2_print <- renderTable({
  if(is.null(s_C_grid_Mat2())) return(NULL)
  s_C_lambdaMat2()
})


### lambdaProgMat --------------------------------------------------------------
s_C_grid_Prog <- reactive({
  grid_comb. =  expand.grid(1:input$s_T_timepoints, 1:input$s_C_subgr)
  paste0("C_Prog_s", grid_comb.$Var1, "t", grid_comb.$Var2)
})

observe({
  isolate({
    output$s_C_lambdaProgMat <-renderTable({
      if(is.null(s_C_grid_Prog())) return(NULL)
      num.inputs.col1 <- 
        paste0("<input id='s_", s_C_grid_Prog(), 
               "' class='shiny-bound-input' ",
               " type='number' max = 999 min = 0 step = 1 value=", 
               18, ">")
      out = as.data.frame(matrix(num.inputs.col1, nrow = input$s_C_subgr))
      
      t1=paste0("s_T_time", (1:input$s_T_timepoints)-1)
      t2=paste0("s_T_time", (1:input$s_T_timepoints))
      names(out) = unlist(lapply(1:input$s_T_timepoints, function(i)
        paste0(input[[t1[i]]], "-", input[[t2[i]]])))
      # names(out) = paste0("T", 1:input$s_T_timepoints)
      rownames(out) = paste0("S", 1:input$s_C_subgr)
      out
    }, sanitize.text.function = function(x) x, 
    rownames = TRUE,
    width = 'auto')
  })
})


s_C_lambdaProgMat <- reactive({
  mm = unlist(lapply(1:length(s_C_grid_Prog()), function(i){
    input[[paste0("s_", s_C_grid_Prog()[i])]]
  }))
  if(is.null(mm)) return(NULL)
  if(length(mm)%%input$s_C_subgr != 0) return(NULL) # Catch error when mm is not ready to be in the matrix
  matrix(fu(mm), nrow = input$s_C_subgr)
})

output$s_C_lambdaProgMat_print <- renderTable({
  if(is.null(s_C_grid_Prog())) return(NULL)
  s_C_lambdaProgMat()
}, digits = 8)



### Analysis -------------------------------------------------------------------


#Calculation of the hazard and survival functions.
#The hazard rates per day before and after PD and the hazard for PD are defined in lambdaMat1, lambdaMat2 and lambdaProgMat
#The matrices contain one row per subgroup and one column per time interval

B5 <- reactive({
  pop_hazVFfun(
    Tint = s_T_times(),
    lambdaMat1 = s_T_lambdaMat1(),
    lambdaMat2 = s_T_lambdaMat2(),
    lambdaProgMat = s_T_lambdaProgMat(),
    p = s_T_p(),
    timezero = FALSE
  )
})

K5 <- reactive({
  pop_hazVFfun(
    Tint = s_T_times(),
    lambdaMat1 = s_C_lambdaMat1(),
    lambdaMat2 = s_C_lambdaMat2(),
    lambdaProgMat = s_C_lambdaProgMat(),
    p = s_C_p(),
  timezero = TRUE
) #Setting timezero=TRUE here means we have the same delayed onset for switchers after progression as we have initially for patients in the treatment group.
})



# Create Plot ------------------------------------------------------------------
react_panplot = eventReactive(input$s_calculate,{
  print(plot_shhr(A = K5(), B = B5()))
})

output$panplot = renderPlot(react_panplot())

output$ctrl_plot = renderPlot({
  plot_diagram(A = K5(),
               B = K5(), 
               s_C_labels(), 
               s_C_labels(), 
               which = "Control")
})
output$trt_plot = renderPlot({
  plot_diagram(A = B5(),
               B = B5(), 
               s_T_labels(), 
               s_T_labels(), 
               which = "Experimental")
})


# Calculate medians ------------------------------------------------------------
k_median = reactive({
  K = K5()
  pos_below = suppressWarnings(min(K$t[(K$S <=  0.5)]) + 1)
  pos_avobe = max(K$t[(K$S >=  0.5)]) + 1
  dd = K$S[pos_avobe] - K$S[pos_below] 
  ((pos_avobe-1)*abs(K$S[pos_below] - 0.5) + (pos_below-1)*abs(K$S[pos_avobe] - 0.5))/dd
})
b_median = reactive({
  K = B5()
  pos_below = suppressWarnings(min(K$t[(K$S <=  0.5)]) + 1)
  pos_avobe = max(K$t[(K$S >=  0.5)]) + 1
  dd = K$S[pos_avobe] - K$S[pos_below] 
  ((pos_avobe-1)*abs(K$S[pos_below] - 0.5) + (pos_below-1)*abs(K$S[pos_avobe] - 0.5))/dd
})

median_table = eventReactive(input$s_calculate,
 data.frame(group = c("Experimental", "Control"),
            median = c(b_median(), k_median()))
)
output$median_table = renderTable({
  median_table()
},
caption = "<h4 style='text-align:center;color:black;'><b>Median survival times by treatment group</b></h4>", 
caption.placement = getOption("xtable.caption.placement", "top"), 
caption.width = getOption("xtable.caption.width", NULL))





# Download code ----------------------------------------------------------------
output$s_downloadCode <- downloadHandler(
  filename = paste('nph-simulate-rcode.r'),
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


# Creat state diagrams
TRT_subgr_labels = %s
CRL_subgr_labels = %s
plot_diagram(A = TRT,
             B = CRL, 
             A_subgr_labels = TRT_subgr_labels,  
             B_subgr_labels = CRL_subgr_labels, 
             which = 'Control')

plot_diagram(A = TRT,
             B = CRL, 
             A_subgr_labels = TRT_subgr_labels,  
             B_subgr_labels = CRL_subgr_labels, 
             which = 'Experimental')


# Create Survival, Hazard and hazard ratio plots
plot_shhr(A = TRT, B = CRL)


", 
    sprintf("c(%s)", paste0(s_T_times(), collapse = ", ")),   # The time intervals
    sprintf("m2r(matrix(c(%s), nrow = %s, ncol = %s))", paste0(m2r(s_T_lambdaMat1()), collapse = ", "), input$s_T_subgr, input$s_T_timepoints),      
    sprintf("m2r(matrix(c(%s), nrow = %s, ncol = %s))", paste0(m2r(s_T_lambdaMat2()), collapse = ", "), input$s_T_subgr, input$s_T_timepoints),      
    sprintf("m2r(matrix(c(%s), nrow = %s, ncol = %s))", paste0(m2r(s_T_lambdaProgMat()), collapse = ", "), input$s_T_subgr, input$s_T_timepoints),      
    sprintf("c(%s)", paste0(s_T_p(), collapse = ", ")),       # The proportions in the subgroups for treatment arm
    sprintf("m2r(matrix(c(%s), nrow = %s, ncol = %s))", paste0(m2r(s_T_lambdaMat1()), collapse = ", "), input$s_T_subgr, input$s_T_timepoints),      
    sprintf("m2r(matrix(c(%s), nrow = %s, ncol = %s))", paste0(m2r(s_T_lambdaMat2()), collapse = ", "), input$s_T_subgr, input$s_T_timepoints),      
    sprintf("m2r(matrix(c(%s), nrow = %s, ncol = %s))", paste0(m2r(s_T_lambdaProgMat()), collapse = ", "), input$s_T_subgr, input$s_T_timepoints),      
    sprintf("c(%s)", paste0(s_C_p(), collapse = ", ")),       # The proportions in the subgroups for control arm
    sprintf("c(%s)", paste0(paste0("'", s_T_labels(), "'"), collapse = ", ")),       # The proportions in the subgroups for control arm
    sprintf("c(%s)", paste0(paste0("'", s_C_labels(), "'"), collapse = ", "))       # The proportions in the subgroups for control arm
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
