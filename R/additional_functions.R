#' Transform median time into rate
#' 
#' This helper function calculates the hazard rate per day of an exponential 
#' distribution from the median given in months.
#' 
#' @param x The median time in months to be transformed into rate
#' 
#' @export
m2r <- function(x){
  12 / 365.25 * log(2) / x
}

#' Draw a state space figure
#' 
#' A figure that shows the states and the possible transitions between them.
#' 
#' @param A An object of class \code{mixpch}, resembling the survival function in treatment group 0
#' @param B An object of class \code{mixpch}, resembling the survival function in treatment group 1
#' @param A_subgr_labels A character vector with the same length as A$p. It indicates names for the subgroups in A
#' @param B_subgr_labels A character vector with the same length as B$p. It indicates names for the subgroups in B
#' @param which Which treatment arm should be shown? One of "Both", "Experimental", "Control".
#' @param treatment_labels A character vector of length 2 indicating the treatment labels.
#' @param colors Either a vector of length two with colors for A and B, or "default".
#' @param show.rate A logical indicating whether the rate should be shown in the diagram
#' 
#' @import ggplot2
#' @export
plot_diagram = function(A,
                        B,
                        A_subgr_labels = "",
                        B_subgr_labels = "",
                        which = c("Both", "Experimental", "Control"),
                        treatment_labels = c("Experimental", "Control"),
                        colors = "default", show.rate = TRUE){
  which = match.arg(which)
  
  times = A$Tint
  t_pr = A$p
  t_label = A_subgr_labels
  t_lambdaMat1 = A$lambdaMat1
  t_lambdaMat2 = A$lambdaMat2
  t_lambdaProgMat = A$lambdaProgMat
  c_pr = B$p
  c_label = B_subgr_labels
  c_lambdaMat1 = B$lambdaMat1
  c_lambdaMat2 = B$lambdaMat2
  c_lambdaProgMat = B$lambdaProgMat
  
  if (length(times) > 3) stop("Only two intervals are supported")
  
  t_group = rep(t_label, each = 3)
  c_group = rep(c_label, each = 3)
  t_group_id = rep(1:length(t_pr), each = 3)
  c_group_id = rep(1:length(c_pr), each = 3)
  
  dat = data.frame(treatment = c(rep("Experimental", 3*length(t_pr)),
                                 rep("Control", 3*length(c_pr))),
                   treatment_labels = c(rep(treatment_labels[1], 3*length(t_pr)),
                                        rep(treatment_labels[2], 3*length(c_pr))),
                   group = c(t_group, c_group),
                   group_id = c(t_group_id, c_group_id))
  dat$x = c(-1, 0, 1) * 4
  dat$y = c( 1, 0, 1) * 2
  dat$x_begin = c(-1,        0 + 1/4,  1 - 1/4) * 4
  dat$y_begin = c( 1 - .5/2, 0,        1) * 2
  dat$x_end   = c( 0 -  1/4, 1      , -1 + 1/4) * 4
  dat$y_end   = c( 0,        1 -.5/2,  1) * 2
  dat$x_mid   = dat$x_begin + dat$x_end
  dat$y_mid   = dat$y_begin + dat$y_end
  
  dat$ends = c("last","last","first")
  dat$label = c("Alive\nNo progression","Progression","Death")
  dat$lambdas = c("lambda[P](t)",
                  "lambda[PD](t)",
                  "lambda[D](t)")
  dat$nudge_x = c(-.4,  0, 0)
  dat$nudge_y = c(-.4,-.4,.4)
  
  t_lambdas = data.frame()
  for (i in 1:length(t_pr)){
    t_lambdas = rbind(t_lambdas,
                      t_lambdaMat1[i,],
                      t_lambdaMat2[i,],
                      t_lambdaProgMat[i,])
  }
  c_lambdas = data.frame()
  for (i in 1:length(c_pr)){
    c_lambdas = rbind(c_lambdas,
                      c_lambdaMat1[i,],
                      c_lambdaMat2[i,],
                      c_lambdaProgMat[i,])
  }
  colnames(t_lambdas) = colnames(c_lambdas) = paste0("T", 1:(length(times)-1))
  dat = cbind(dat, rbind(t_lambdas, c_lambdas))
  if(show.rate){
    if (length(times) == 2){
      dat$lambda_label = paste0(dat$lambdas, ": ", round(dat$T1, 4))
    } else if (length(times) == 3){
      dat$lambdast1 = c(sprintf("lambda[P]^paste('[%s-%s]')*(t)",  0, times[2]),
                        sprintf("lambda[PD]^paste('[%s-%s]')*(t)", 0, times[2]),
                        sprintf("lambda[D]^paste('[%s-%s]')*(t)",  0, times[2]))
      dat$lambdast2 = c(sprintf("lambda[P]^paste('[%s-%s]')*(t)",  times[2], times[3]),
                        sprintf("lambda[PD]^paste('[%s-%s]')*(t)", times[2], times[3]),
                        sprintf("lambda[D]^paste('[%s-%s]')*(t)",  times[2], times[3]))
      dat$lambda_label = paste0("atop(",dat$lambdast1, ": ", round(dat$T1, 4), ",",
                                dat$lambdast2, ": ", round(dat$T2, 4), ")")
    } 
  } else{
    dat$lambda_label = dat$lambdas
  }
  
  if (which == "Experimental"){
    dat = dat[which(dat$treatment == "Experimental"), ]
  } else if (which == "Control"){
    dat = dat[which(dat$treatment == "Control"), ]      
  }
  
  if(colors == "default"){
    col_trt = "pink" 
    col_ctrl = "lightblue" 
  }
  ggplot(dat) + 
    facet_grid(group ~ treatment_labels) + 
    geom_rect(aes_string(xmin = "x - 1", 
                  xmax = "x + 1",
                  ymin = "y - 0.5", 
                  ymax = "y + 0.5",
                  fill = "treatment"))  +
    geom_text(aes_string(x = "x", 
                         y = "y", label = "label"))  +
    geom_text(aes_string(x = "(x_begin + x_end)/2 + nudge_x", 
                         y = "(y_begin + y_end)/2 + nudge_y",
                         label = "lambda_label"), 
              parse = TRUE)  +
    geom_segment(aes_string(x    = "x_begin", 
                            xend = "x_end" ,
                            y    = "y_begin" , 
                            yend = "y_end"),
                 arrow = arrow(type = "closed", 
                               length = unit(10, "points"),
                               ends = dat$ends)) +
    coord_equal(ylim = c(-0.5, 3.5), xlim = c(-5, 5)) + 
    theme_void() + 
    theme(legend.position = "none", 
          text = element_text(size = 20)) +
    scale_fill_manual(values = c("Experimental" = col_trt, 
                                 "Control" = col_ctrl)) -> p1
  
  p1
}


#' Plot of survival, hazard and hazard ratio of two groups as a function of time
#' 
#' A convenience function that uses the generic plot function in the nph package to
#' plot the three functions in a layout of 3 columns and 1 row.
#' 
#' @param A An object of class \code{mixpch}, resembling the survival function in treatment group 0
#' @param B An object of class \code{mixpch}, resembling the survival function in treatment group 1
#' @param main An overall title for the plot
#' @param xmax A maximum value for the x-axis. The plot is drawn using xlim = c(0, xmax)
#' @param ymax_haz A maximum value for the y-axis for the hazards plot. The plot is drawn using ylim = c(0, ymax_haz)
#' @param ymax_hr  A maximum value for the y-axis for the hazards ratio plot. The plot is drawn using ylim = c(0, ymax_hr)
#' 
#' @import graphics
#' @export
plot_shhr <- function(A, B, main = "", xmax = NULL, ymax_haz = NULL, ymax_hr = NULL) {
  if (is.null(xmax)) xmax = length(A$t)
  CX <- 2
  CX2 <- 1
  dd <- 4
  # Set graphical parameters. Save previous setup in old_par
  old_par = par(
    lwd = dd,
    font = 2,
    font.axis = 2,
    font.lab = 2,
    cex = CX,
    cex.lab = CX,
    cex.axis = CX,
    cex.main = CX * 1.5,
    mfrow = c(1, 3),
    oma = c(0, 0, 2, 0),
    mar = c(5, 5, 3, 2)
  )
  on.exit(par(old_par))
  
  # Plot the survival functions
  plot(A,
       ylim = c(0, 1),
       xlim = c(0,xmax),
       ylab = "Survival",
       xlab = "Days")
  plot(B, add = TRUE, lty = 2, col = 3)
  
  # Plot the hazard functions
  if (is.null(ymax_haz)) ymax_haz = max(B$haz, A$haz, na.rm = TRUE) + 0.001
  plot(A,
       "haz",
       xlim = c(0, xmax),
       ylim = c(0, ymax_haz),
       ylab = "Hazard",
       xlab = "Days")
  plot(B,
       "haz",
       add = TRUE,
       lty = 2,
       col = 3)
  
  # Plot the hazard ratio
  if (is.null(ymax_hr)) ymax_hr = max(B$haz / A$haz, na.rm = TRUE) + 1
  plot(
    A$t,
    B$haz / A$haz,
    xlim = c(0, xmax),
    ylim = c(0, ymax_hr),
    ylab = "Hazard ratio",
    xlab = "Days",
    type = "l"
  )
  
  # Add title
  mtext(
    outer = TRUE,
    side = 3,
    text = main, cex = 2,
    line = -2
  )

}


# Generic print function for the weighted log-rank tests
#' @export
print.wlogrank = function (x, digits = max(options()$digits - 4, 3), ...) {
  saveopt <- options(digits = digits)
  on.exit(options(saveopt))
  
  otmp <- x$obs
  etmp <- x$exp
  df <- (sum(1 * (etmp > 0))) - 1
  temp <- cbind(x$n, otmp, etmp,
                ((otmp - etmp)^2)/etmp, 
                ((otmp - etmp)^2)/x$var)
  dimnames(temp) <- list(names(x$n), c("N", "Observed", 
                                       "Expected", "(O-E)^2/E", "(O-E)^2/V"))
  # Print everything
  cat("Call:\n")
  dput(x$call)
  cat("\n")
  print(data.frame(temp, check.names = FALSE))
  cat("\n Chisq=", format(round(x$test$Chisq, 1)), " on", df, 
      "degrees of freedom, p=", format.pval(pchisq(x$test$Chisq, 
                                                   df, lower.tail = FALSE)))
  cat("\n rho   = ", x$test$rho,
      "gamma = ", x$test$gamma, "\n")
  invisible(x)
}


# Generic print function for the maximum of weighted log-rank tests
#' @export
print.wlogrank_max = function (x, digits = max(options()$digits - 4, 3), ...) {
  saveopt <- options(digits = digits)
  on.exit(options(saveopt))
  
  # Print everything -
  cat("Call:\n")
  dput(x$call)
  cat("\n Two sided p-value =", format(x$pmult), 
      paste0("(Bonferroni corrected: ", format(x$p.Bonf), ")\n"))
  cat("\n Individual weighted log-rank tests:\n") 
  print(x$tests)
  invisible(x)
}



#' Launch a GUI (shiny app) for the nph package
#'
#' @details The packages \code{shinycssloaders}, \code{formatR} and \code{styler} are required for correct display of the GUI.
#' The package \code{rmarkdown} with access to pandoc is required to save reports.
#' @export
nph_gui <- function(){
  appDir <- system.file("shiny", package = "nph")
  shiny::runApp(appDir, display.mode = "normal");
}

#' Draw a population composition plot
#' 
#' A figure that shows the composition of the population under study though time
#' 
#' @param A An object of class \code{mixpch}, resembling the survival function in treatment group 0
#' @param B An object of class \code{mixpch}, resembling the survival function in treatment group 1
#' @param colors Either a vector of length four with colors for A and B and subgroup 1 and 2, or "default".
#' @param max_time the maximum value for the x-axis.
#' @param position Either "stack" or "fill". By default (stack), the total population decreases through time. If position="fill", the size of the population is rescaled to show conditional percentages.
#' @param title The text for the title.
#' 
#' @author Robin Ristl, \email{robin.ristl@@meduniwien.ac.at}, Nicolas Ballarini
#' @seealso \code{\link{pop_pchaz}}
#' @examples
#' A <- pop_pchaz(Tint = c(0, 90, 365),
#'   lambdaMat1 = matrix(c(0.2, 0.1, 0.4, 0.1), 2, 2) / 365,
#'	 lambdaMat2 = matrix(c(0.5, 0.2, 0.6, 0.2), 2, 2) / 365,
#'	 lambdaProg = matrix(c(0.5, 0.5, 0.4, 0.4), 2, 2) / 365,
#'	 p = c(0.8, 0.2), 
#'	 timezero = FALSE, discrete_approximation = TRUE)
#' B <- pop_pchaz(Tint = c(0, 90, 365),
#'   lambdaMat1 = matrix(c(0.2, 0.1, 0.4, 0.1), 2, 2) / 365,
#'	 lambdaMat2 = matrix(c(0.5, 0.1, 0.6, 0.1), 2, 2) / 365,
#'	 lambdaProg = matrix(c(0.5, 0.5, 0.04, 0.04), 2, 2) / 365,
#'	 p = c(0.8, 0.2), 
#'	 timezero = FALSE, discrete_approximation = TRUE)
#' plot_subgroups(A, B, title = "position='stack'")
#' plot_subgroups(A, B, position='fill', title = "position='fill'")
#'
#'
#' @import ggplot2
#' @export
plot_subgroups = function(A,
                          B,
                          colors = "default",
                          max_time = max(A$Tint), 
                          position = c("stack", "fill"),
                          title = ""){

  position = match.arg(position)
  if(colors == "default"){
    cols_new = c("#FAAB18", "#fccc74", "#1380A1", "#71b2c6")
  } else{
    if(length(colors)!=4) stop("colors should be of length 4")
    cols_new = colors
  }
  
  A_n_subgr = length(A$p)
  B_n_subgr = length(B$p)
  
  t = A$t
  if(any(A_n_subgr!=2, B_n_subgr!=2)) stop("This function is only implemented for two subgroups in each treatment")
  
  nsubpop = 4
  lambdaProgMat = rbind(A$lambdaProgMat, B$lambdaProgMat)
  lambdaMat1 = rbind(A$lambdaMat1, B$lambdaMat1)
  lambdaMat2 = rbind(A$lambdaMat2, B$lambdaMat2)
  timezero = c(A$timezero, B$timezero)
  p = c(A$p, B$p)/2
  Tint = A$Tint
  S = S_scaled = list()
  # A$
  for(i in 1:nsubpop) {
    lambdaProg = lambdaProgMat[i,]
    lambda1 = lambdaMat1[i,]
    lambda2 = lambdaMat2[i,]
    timezero[i]
    if(all(lambdaProg == Inf)){
      funs = hazVFfun(Tint, lambda2)
    } else if(all(lambdaProg == 0)){
      funs = hazVFfun(Tint, lambda1)
    } else {
      funs = subpop_hazVFfun(Tint, lambda1, lambda2, lambdaProg, timezero[i])
    }
    S[[i]] = 1-funs$F
    S_scaled[[i]] = (1-funs$F)*p[i]
  }
  
  Smat = do.call(rbind,S)
  Smat_scaled = do.call(rbind,S_scaled)
  Smat[, 1:10]
  Smat_scaled[, 1:10]
  
  dat = do.call(rbind,lapply(1:ncol(Smat_scaled), function(x){
    values = Smat_scaled[, x]
    if(position=="fill") values = values / sum(values)
    data.frame(time = x,
               treatment = rep(c(1,2), each = 2),
               subgroup =  rep(c(1,2), 2),
               value =  values)
  }))
  comb = paste0(dat$treatment, dat$subgroup)
  dat$group = factor(comb, levels = rev(unique(comb)))
  names(cols_new) = levels(dat$group)
  
  # Calculate percentages for labels
  y1p = Smat_scaled[, 1]
  y1b1 = c(0,cumsum(y1p))
  y1b = c(0,cumsum(y1p))
  y1b = y1b[-length(y1b)] + y1p/2
  y1p = y1p*100
  
  y2p = Smat_scaled[, ncol(Smat_scaled)]
  if(position=="fill") y2p = y2p/sum(y2p)
  y2b1 = c(0,cumsum(y2p))
  y2b  = c(0,cumsum(y2p))
  y2b = y2b[-length(y2b)] + y2p/2
  y2p = y2p*100
  
  ggplot(dat) +
    geom_area(aes_string(x = "time", fill = "group", y = "value")) + 
    scale_fill_manual(values = cols_new, 
                      labels = c("CTRL - Subgr 2", "CTRL - Subgr 1", 
                                 "TRT - Subgr 2", "TRT - Subgr 1")) +
    guides(fill = guide_legend(title = " ", reverse=TRUE)) +
    coord_cartesian(expand = FALSE, clip="off") + 
    labs(x = "Time", y = "Percent", title = title) +
    theme(panel.background = element_rect(fill = "gray50"),
          panel.grid = element_blank(),
          legend.position = "bottom",
          axis.text.y=element_text(margin = margin(r = 10)),
          axis.text.y.right = element_text(margin = margin(l = 10)),
          axis.line.x.bottom = element_line(size = 1)) +
    scale_x_continuous(limits = c(0,max_time), breaks = seq(0, max_time, 500)) + 
    scale_y_continuous(breaks = y1b,
                       labels =  sprintf("%.1f%%", y1p), 
                       sec.axis = sec_axis(~., 
                                           breaks = y2b, 
                                           labels = sprintf("%.1f%%", y2p))) -> pppp
  pppp + theme(axis.ticks.length = unit(0.5, "cm"),
               axis.ticks.x = element_blank(),             
               axis.ticks.y = element_line(colour = "gray40"),
               axis.text.y = element_text(size = 11),
               axis.text.x = element_text(size = 12, color = "black", margin = margin(t = -10)))
    
}



