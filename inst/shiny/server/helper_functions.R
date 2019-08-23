# This helper function calculates the hazard rate per day of an exponential 
# distribution from the median given in months.
fu <- function(x){
  12 / 365.25 * log(2) / x
}
#A function to show survival, hazard and hazard ratio of two groups as function of time
panplot <- function(K,
                    B,
                    main = "",
                    savepdf = NA,
                    savepng = NA) {
  if (!is.na(savepdf))
    pdf(paste(savepdf, ".pdf", sep = ""),
        width = 10,
        height = 4)
  if (!is.na(savepng))
    jpeg(
      paste(savepng, ".png", sep = ""),
      width = 25,
      height = 8,
      unit = "cm",
      res = 300
    )
  CX <- 2
  CX2 <- 1
  dd <- 4
  par(
    lwd = dd,
    font = 2,
    font.axis = 2,
    font.lab = 2,
    cex = CX,
    cex.lab = CX,
    cex.axis = CX,
    cex.main = CX * 1.5
  )
  par(
    mfrow = c(1, 3),
    oma = c(0, 0, 2, 0),
    mar = c(5, 5, 3, 2)
  )
  plot(K,
       ylim = c(0, 1),
       ylab = "Survival",
       xlab = "Days")
  plot(B, add = TRUE, lty = 2, col = 3)
  
  ymax = max(B$haz, K$haz, na.rm = TRUE)
  cat(ymax)
  plot(K,
       "haz",
       ylim = c(0, ymax + 0.001),
       ylab = "Hazard",
       xlab = "Days")
  plot(B,
       "haz",
       add = TRUE,
       lty = 2,
       col = 3)
  
  plot(
    K$t,
    B$haz / K$haz,
    ylim = c(0, max(B$haz / K$haz, na.rm = TRUE) + 1),
    ylab = "Hazard ratio",
    xlab = "Days",
    type = "l"
  )
  mtext(
    outer = TRUE,
    side = 3,
    text = main,
    line = -2
  )#,cex=CX2,adj=0)
  if (!is.na(savepng))
    dev.off()
  if (!is.na(savepdf))
    dev.off()
}




#A function to show survival, hazard and hazard ratio of two groups as function of time
intplot <- function(dat, 
                    K,
                    B,
                    main = "",
                    savepdf = NA,
                    savepng = NA) {
  if (!is.na(savepdf))
    pdf(paste(savepdf, ".pdf", sep = ""),
        width = 10,
        height = 4)
  if (!is.na(savepng))
    jpeg(
      paste(savepng, ".png", sep = ""),
      width = 25,
      height = 8,
      unit = "cm",
      res = 300
    )
  CX <- 1
  CX2 <- 1
  dd <- 4
  par(
    lwd = dd,
    font = 2,
    font.axis = 2,
    font.lab = 2,
    cex = CX,
    cex.lab = CX,
    cex.axis = CX,
    cex.main = CX * 1.5
  )
  par(
    # mfrow = c(1, 3),
    oma = c(0, 0, 2, 0),
    mar = c(5, 5, 3, 2)
  )
  km.by.trt <- survival::survfit(survival::Surv(y, event) ~ group, data = dat)
  plot(km.by.trt, col = c(1,3),
       lwd = dd,
       ylab = "Survival",
       xlab = "Days")
  plot(K, add = TRUE, 
       lwd = dd/2,
       ylim = c(0, 1),
       ylab = "Survival",
       xlab = "Days")
  plot(B, 
       lwd = dd/2,
       add = TRUE, lty = 2, col = 3)
  
  mtext(
    outer = TRUE,
    side = 3,
    text = main,
    line = -2
  )#,cex=CX2,adj=0)
  if (!is.na(savepng))
    dev.off()
  if (!is.na(savepdf))
    dev.off()
}





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
#' @param colors Either a vector of length two with colors for A and B, of "default".
#' @param show.rate A logical indicating whether the rate should be shown in the diagram
#' 
#' @export
plot_diagram = function(A,
                        B,
                        A_subgr_labels,
                        B_subgr_labels,
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
    dat %>% 
      filter(treatment == "Experimental") -> dat
  } else if (which == "Control"){
    dat %>% 
      filter(treatment == "Control") -> dat
  }
  
  if(colors == "default"){
    col_trt = "pink" 
    col_ctrl = "lightblue" 
  }
  ggplot(dat) + 
    facet_grid(group ~ treatment_labels) + 
    geom_rect(aes(xmin = x - 1, 
                  xmax = x + 1,
                  ymin = y - 0.5, 
                  ymax = y + 0.5,
                  fill = treatment))  +
    geom_text(aes(x = x, 
                  y = y, label = label))  +
    geom_text(aes(x = (x_begin + x_end)/2 + nudge_x, 
                  y = (y_begin + y_end)/2 + nudge_y,
                  label = lambda_label), 
              parse = TRUE)  +
    geom_segment(aes(x    = x_begin, 
                     xend = x_end ,
                     y    = y_begin , 
                     yend = y_end),
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
#' @param xmax A maximum value for the x-axis. The plots is drawn using xlim = c(0, xmax)
#' 
#' @export
plot_shhr <- function(A, B, main = "", xmax = NULL) {
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
  
  # Plot the survival functions
  plot(A,
       ylim = c(0, 1),
       xlim = c(0,xmax),
       ylab = "Survival",
       xlab = "Days")
  plot(B, add = TRUE, lty = 2, col = 3)
  
  # Plot the hazard functions
  ymax = max(B$haz, A$haz, na.rm = TRUE)
  plot(A,
       "haz",
       xlim = c(0,xmax),
       ylim = c(0, ymax + 0.001),
       ylab = "Hazard",
       xlab = "Days")
  plot(B,
       "haz",
       add = TRUE,
       lty = 2,
       col = 3)
  
  # Plot the hazard ratio
  plot(
    A$t,
    B$haz / A$haz,
    xlim = c(0,xmax),
    ylim = c(0, max(B$haz / A$haz, na.rm = TRUE) + 1),
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
  par(old_par)
}
