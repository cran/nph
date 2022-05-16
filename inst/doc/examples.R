## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(cache = FALSE,
  collapse = TRUE,
  comment = "#>"
)
options()

## ----fig.show='hold'----------------------------------------------------------
library(nph)

## ----001, fig.show='hold', fig.width=12.5, warning=FALSE, out.width='100%'----
times <- c(0, 5 * 365)   # Time interval boundaries, in days

# Treatment group
t_resp <- 1              # There are no subgroups
B5 <- pop_pchaz(
  T = times,
  lambdaMat1    = m2r(matrix(18,nrow = 1)),
  lambdaMat2    = m2r(matrix(11,nrow = 1)),
  lambdaProgMat = m2r(matrix(9, nrow = 1)),
  p = t_resp,
  timezero = FALSE, discrete_approximation = TRUE
)

# Control group
c_resp <- 1              # There are no subgroups
K5  <- pop_pchaz(
  T = times,
  lambdaMat1    = m2r(matrix(11,nrow = 1)),
  lambdaMat2    = m2r(matrix(9, nrow = 1)),
  lambdaProgMat = m2r(matrix(5, nrow = 1)),
  p = c_resp,
  timezero = TRUE, discrete_approximation = TRUE
)

pp = plot_diagram(B5, K5)
pp

## ----fig.show='hold', fig.height=3.5, fig.width=10, warning=FALSE, out.width='100%'----
plot_shhr(K5, B5)

## ----002,  fig.show='hold', fig.width=12.5, fig.height=10.5, warning=FALSE, out.width='100%'----
times <- c(0, 5 * 365)   # Time interval boundaries, in days

# Treatment group
t_resp <-  c(.2, .8)
B5 <- pop_pchaz(
  T = times,
  lambdaMat1    = m2r(matrix(c(30,
                               18),nrow = 2)),
  lambdaMat2    = m2r(matrix(c(20,
                               11),nrow = 2)),
  lambdaProgMat = m2r(matrix(c(15, 
                               9), nrow = 2)),
  p = t_resp,
  timezero = FALSE, discrete_approximation = TRUE
)

# Control group
c_resp <- 1 
K5  <- pop_pchaz(
  T = times,
  lambdaMat1    = m2r(matrix(11,nrow = 1)),
  lambdaMat2    = m2r(matrix(9, nrow = 1)),
  lambdaProgMat = m2r(matrix(5, nrow = 1)),
  p = c_resp,
  timezero = TRUE, discrete_approximation = TRUE
)

pp = plot_diagram(B5, K5, A_subgr_labels = c("S1", "S2"))
pp

## ----fig.show='hold', fig.height=3.5, fig.width=10, warning=FALSE, out.width='100%'----
plot_shhr(K5, B5)

## ----003,  fig.show='hold', fig.width=11, warning=FALSE, out.width='100%'-----
times <- c(0, 5 * 365)   # Time interval boundaries, in days

# Treatment group
t_resp <-  c(.5, .5)
B5 <- pop_pchaz(
  T = times,
  lambdaMat1    = m2r(matrix(c(30,
                               18),nrow = 2)),
  lambdaMat2    = m2r(matrix(c(20,
                               11),nrow = 2)),
  lambdaProgMat = m2r(matrix(c(15, 
                               9), nrow = 2)),
  p = t_resp,
  timezero = FALSE, discrete_approximation = TRUE
)

# Control group
c_resp <- 1 
K5  <- pop_pchaz(
  T = times,
  lambdaMat1    = m2r(matrix(11,nrow = 1)),
  lambdaMat2    = m2r(matrix(9, nrow = 1)),
  lambdaProgMat = m2r(matrix(5, nrow = 1)),
  p = c_resp,
  timezero = TRUE, discrete_approximation = TRUE
)

## ----fig.show='hold', fig.height=3.5, fig.width=10, warning=FALSE, out.width='100%'----
plot_shhr(K5, B5)

## ----004,  fig.show='hold', fig.width=12.5, fig.height=10.5, warning=FALSE, out.width='100%'----
times <- c(0, 5 * 365)   # Time interval boundaries, in days

# Treatment group
t_resp <-  c(.5, .5)
B5 <- pop_pchaz(
  T = times,
  lambdaMat1    = m2r(matrix(c(30,
                               18), nrow = 2)),
  lambdaMat2    = m2r(matrix(c(30,
                               11), nrow = 2)),
  lambdaProgMat = m2r(matrix(c(15, 
                                9), nrow = 2)),
  p = t_resp,
  timezero = FALSE, discrete_approximation = TRUE
)

# Control group
c_resp <- 1 
K5  <- pop_pchaz(
  T = times,
  lambdaMat1    = m2r(matrix(11, nrow = 1)),
  lambdaMat2    = m2r(matrix( 9, nrow = 1)),
  lambdaProgMat = m2r(matrix( 5, nrow = 1)),
  p = c_resp,
  timezero = TRUE, discrete_approximation = TRUE
)


pp = plot_diagram(B5, K5, A_subgr_labels = c("S1", "S2"))
pp

## ----fig.show='hold', fig.height=3.5, fig.width=10, warning=FALSE, out.width='100%'----
plot_shhr(K5, B5)

## ----005, fig.show='hold', fig.width=12.5, warning=FALSE, out.width='100%'----
times <- c(0, 100, 5 * 365)   # Time interval boundaries, in days

# Treatment group
t_resp <- 1               # There are no subgroups
B5 <- pop_pchaz(
  T = times,
  lambdaMat1    = m2r(matrix(c(11, 18), nrow = 1)),
  lambdaMat2    = m2r(matrix(c( 9, 11), nrow = 1)),
  lambdaProgMat = m2r(matrix(c( 5,  9), nrow = 1)),
  p = t_resp,
  timezero = FALSE, discrete_approximation = TRUE
)

# Control group
c_resp <- 1 
K5  <- pop_pchaz(
  T = times,
  lambdaMat1    = m2r(matrix(c(11, 11), nrow = 1)),
  lambdaMat2    = m2r(matrix(c( 9,  9), nrow = 1)),
  lambdaProgMat = m2r(matrix(c( 5,  5), nrow = 1)),
  p = c_resp,
  timezero = TRUE, discrete_approximation = TRUE
)


pp = plot_diagram(B5, K5)
pp

## ----fig.show='hold', fig.height=3.5, fig.width=10, warning=FALSE, out.width='100%'----
plot_shhr(K5, B5)

## ----006, fig.show='hold', fig.width=12.5, fig.height=10.5, warning=FALSE, out.width='100%'----
times <- c(0, 100, 5 * 365)   # Time interval boundaries, in days

# Treatment group
t_resp <- c(0.5, 0.5) #Proportion of subgroups
B5 <- pop_pchaz(
  T = times,
    lambdaMat1    = m2r(matrix(c(11, 30,
                                 11, 18), byrow = TRUE, nrow = 2)),
    lambdaMat2    = m2r(matrix(c( 9, 20,
                                  9, 11), byrow = TRUE, nrow = 2)),
    lambdaProgMat = m2r(matrix(c( 5, 15,
                                  5,  9), byrow = TRUE, nrow = 2)),
  p = t_resp,
  timezero = FALSE, discrete_approximation = TRUE
)

# Control group
c_resp <- 1 
K5  <- pop_pchaz(
  T = times,
  lambdaMat1    = m2r(matrix(c(11, 11), nrow = 1, ncol = 2)),
  lambdaMat2    = m2r(matrix(c( 9,  9), nrow = 1, ncol = 2)),
  lambdaProgMat = m2r(matrix(c( 5,  5), nrow = 1, ncol = 2)),
  p = c_resp,
  timezero = TRUE, discrete_approximation = TRUE
)

pp = plot_diagram(B5, K5, A_subgr_labels = c("S1", "S2"))
pp

## ----fig.show='hold', fig.height=3.5, fig.width=10, warning=FALSE, out.width='100%'----
plot_shhr(K5, B5)

## ----007, fig.show='hold', fig.width=12.5, fig.height=10.5, warning=FALSE, out.width='100%'----
times <- c(0, 5 * 365)   # Time interval boundaries, in days

# Treatment group
t_resp <- 1               # There are no subgroups
B5 <- pop_pchaz(
  T = times,
  lambdaMat1    = m2r(matrix(18, nrow = 1)),
  lambdaMat2    = m2r(matrix(11, nrow = 1)),
  lambdaProgMat = m2r(matrix( 9, nrow = 1)),
  p = t_resp,
  timezero = FALSE, discrete_approximation = TRUE
)

# Control group
c_resp <- c(1/3, 2/3) #non-switcher, switcher reponse
K5  <- pop_pchaz(
  T = times, 
  lambdaMat1    = m2r(matrix(c(11, 
                               11), nrow = 2, ncol = 1)),
  lambdaMat2    = m2r(matrix(c( 9, 
                               11), nrow = 2, ncol = 1)),
  lambdaProgMat = m2r(matrix(c( 5,
                                5), nrow = 2, ncol = 1)),
  p = c_resp,
  timezero = TRUE, discrete_approximation = TRUE
)


pp = plot_diagram(B5, K5, B_subgr_labels = c("Non-switcher","Switcher"))
                  
pp

## ----fig.show='hold', fig.height=3.5, fig.width=10, warning=FALSE, out.width='100%'----
plot_shhr(K5, B5)

## ----007b, fig.show='hold', fig.width=12.5, fig.height=10.5, warning=FALSE, out.width='100%'----
times <- c(0, 5 * 365)   # Time interval boundaries, in days

# Treatment group
t_resp <- 1               # There are no subgroups
B5 <- pop_pchaz(
  T = times,
  lambdaMat1    = m2r(matrix(18, nrow = 1)),
  lambdaMat2    = m2r(matrix(11, nrow = 1)),
  lambdaProgMat = m2r(matrix( 9, nrow = 1)),
  p = t_resp,
  timezero = FALSE, discrete_approximation = TRUE
)

# Control group
c_resp <- c(1/3, 2/3) #non-switcher, switcher reponse
K5  <- pop_pchaz(
  T = times, 
  lambdaMat1    = m2r(matrix(c(11, 
                               11), nrow = 2, ncol = 1)),
  lambdaMat2    = m2r(matrix(c( 9, 
                               14), nrow = 2, ncol = 1)),
  lambdaProgMat = m2r(matrix(c( 5,
                                5), nrow = 2, ncol = 1)),
  p = c_resp,
  timezero = TRUE, discrete_approximation = TRUE
)

pp = plot_diagram(B5, K5, B_subgr_labels = c("Non-switcher","Switcher"))
pp

## ----fig.show='hold', fig.height=3.5, fig.width=10, warning=FALSE, out.width='100%'----
plot_shhr(K5, B5)

## ----008,  fig.show='hold', fig.width=12.5, fig.height=10.5, warning=FALSE, out.width='100%'----
times <- c(0, 5 * 365)   # Time interval boundaries, in days

# Treatment group
t_resp <- 1               # There are no subgroups
B5 <- pop_pchaz(
  T = times,
  lambdaMat1    = m2r(matrix(18,nrow = 1)),
  lambdaMat2    = m2r(matrix(11,nrow = 1)),
  lambdaProgMat = m2r(matrix(9, nrow = 1)),
  p = t_resp,
  timezero = FALSE, discrete_approximation = TRUE
)

# Control group
c_resp <- c(2/3, 1/3) #non-switcher, switcher reponse
K5  <- pop_pchaz(
  T = times, 
  lambdaMat1    = m2r(matrix(c(11, 
                               11), nrow = 2, ncol = 1)),
  lambdaMat2    = m2r(matrix(c(9, 
                               11), nrow = 2, ncol = 1)),
  lambdaProgMat = m2r(matrix(c(5,
                               5),  nrow = 2, ncol = 1)),
  p = c_resp,
  timezero = TRUE, discrete_approximation = TRUE
)

pp = plot_diagram(B5, K5, B_subgr_labels = c("Non-switcher","Switcher"))
pp

## ----fig.show='hold', fig.height=3.5, fig.width=10, warning=FALSE, out.width='100%'----
plot_shhr(K5, B5)

## ----008b,  fig.show='hold', fig.width=12.5, fig.height=10.5, warning=FALSE, out.width='100%'----
times <- c(0, 5 * 365)   # Time interval boundaries, in days

# Treatment group
t_resp <- 1               # There are no subgroups
B5 <- pop_pchaz(
  T = times,
  lambdaMat1    = m2r(matrix(18, nrow = 1)),
  lambdaMat2    = m2r(matrix(11, nrow = 1)),
  lambdaProgMat = m2r(matrix( 9, nrow = 1)),
  p = t_resp,
  timezero = FALSE, discrete_approximation = TRUE
)

# Control group
c_resp <- c(2/3, 1/3) #non-switcher, switcher reponse
K5  <- pop_pchaz(
  T = times, 
  lambdaMat1    = m2r(matrix(c(11, 
                               11), nrow = 2, ncol = 1)),
  lambdaMat2    = m2r(matrix(c( 9, 
                               14), nrow = 2, ncol = 1)),
  lambdaProgMat = m2r(matrix(c( 5,
                                5), nrow = 2, ncol = 1)),
  p = c_resp,
  timezero = TRUE, discrete_approximation = TRUE
)

pp = plot_diagram(B5, K5, B_subgr_labels = c("Non-switcher","Switcher"))
pp

## ----fig.show='hold', fig.height=3.5, fig.width=10, warning=FALSE, out.width='100%'----
plot_shhr(K5, B5)

## ----009,  fig.show='hold', fig.width=12.5, fig.height=13, warning=FALSE, out.width='100%'----
times <- c(0, 5 * 365)   # Time interval boundaries, in days

# Treatment group
t_resp <- c(.5, .5)
B5 <- pop_pchaz(
  T = times,
    lambdaMat1    = m2r(matrix(c( 2,
                                 18), nrow = 2)),
    lambdaMat2    = m2r(matrix(c( 2,
                                 18), nrow = 2)),
    lambdaProgMat = m2r(matrix(c(19, 
                                 19), nrow = 2)),
  p = t_resp,
  timezero = FALSE, discrete_approximation = TRUE
)

# Control group
c_resp <- c(1) 
K5  <- pop_pchaz(
  T = times,
  lambdaMat1    = m2r(matrix(c(11), nrow = 1, ncol = 1)),
  lambdaMat2    = m2r(matrix(c( 9), nrow = 1, ncol = 1)),
  lambdaProgMat = m2r(matrix(c( 5), nrow = 1, ncol = 1)),
  p = c_resp,
  timezero = TRUE, discrete_approximation = TRUE
)

pp = plot_diagram(B5, K5,
                  A_subgr_labels = c("S1", "S2"))
pp

## ----fig.show='hold', fig.height=3.5, fig.width=10, warning=FALSE, out.width='100%'----
plot_shhr(K5, B5)

## ----010, fig.show='hold', fig.width=12.5, fig.height=13.5, warning=FALSE, out.width='100%'----
times <- c(0, 5 * 365)   # Time interval boundaries, in days

# Treatment group
t_resp <- c(.2, .8)
B5 <- pop_pchaz(
  T = times,
  lambdaMat1    = m2r(matrix(c(30,
                               18), nrow = 2)),
  lambdaMat2    = m2r(matrix(c(30,
                               11), nrow = 2)),
  lambdaProgMat = m2r(matrix(c(15, 
                                9), nrow = 2)),
  p = t_resp,
  timezero = FALSE, discrete_approximation = TRUE
)

# Control group
c_resp <- c(2/3, 1/3*t_resp[1], 1/3*t_resp[2])
K5  <- pop_pchaz(
  T = times, 
  lambdaMat1    = m2r(matrix(c(11), nrow = 3, ncol = 1)),   
  lambdaMat2    = m2r(matrix(c( 9, 
                               25,
                               14), nrow = 3)), 
  lambdaProgMat = m2r(matrix(c( 5), nrow = 3,ncol = 1)),
  p = c_resp,
  timezero = TRUE, discrete_approximation = TRUE
)

pp = plot_diagram(B5, K5,
                  A_subgr_labels = c("S1", "S2"),
                  B_subgr_labels = c("Non-switcher","Switcher-S1","Switcher-S2"))
pp

## ----fig.show='hold', fig.height=3.5, fig.width=10, warning=FALSE, out.width='100%'----
plot_shhr(K5, B5)

## ----010b, fig.show='hold', fig.width=12.5, fig.height=13.5, warning=FALSE, out.width='100%'----
times <- c(0, 5 * 365)   # Time interval boundaries, in days

# Treatment group
t_resp <- c(.5, .5)
B5 <- pop_pchaz(
  T = times,
  lambdaMat1    = m2r(matrix(c(30,
                               18), nrow = 2)),
  lambdaMat2    = m2r(matrix(c(30,
                               11), nrow = 2)),
  lambdaProgMat = m2r(matrix(c(15, 
                                9), nrow = 2)),
  p = t_resp,
  timezero = FALSE, discrete_approximation = TRUE
)

# Control group
c_resp <- c(2/3, 1/3*t_resp[1], 1/3*t_resp[2])
K5  <- pop_pchaz(
  T = times, 
  lambdaMat1    = m2r(matrix(c(11), nrow = 3, ncol = 1)),   
  lambdaMat2    = m2r(matrix(c( 9, 
                               25,
                               14), nrow = 3)), 
  lambdaProgMat = m2r(matrix(c( 5), nrow = 3,ncol = 1)),
  p = c_resp,
  timezero = TRUE, discrete_approximation = TRUE
)

pp = plot_diagram(B5, K5,
                  A_subgr_labels = c("S1", "S2"),
                  B_subgr_labels = c("Non-switcher","Switcher-S1","Switcher-S2"))
pp

## ----fig.show='hold', fig.height=3.5, fig.width=10, warning=FALSE, out.width='100%'----
plot_shhr(K5, B5)

## ----010c, fig.show='hold', fig.width=12.5, fig.height=13.5, warning=FALSE, out.width='100%'----
times <- c(0, 5 * 365)   # Time interval boundaries, in days

# Treatment group
t_resp <- c(.5, .5)
B5 <- pop_pchaz(
  T = times,
  lambdaMat1    = m2r(matrix(c(30,
                               18), nrow = 2)),
  lambdaMat2    = m2r(matrix(c(20,
                               11), nrow = 2)),
  lambdaProgMat = m2r(matrix(c(15, 
                                9), nrow = 2)),
  p = t_resp,
  timezero = FALSE, discrete_approximation = TRUE
)

# Control group
c_resp <- c(1/3, 2/3*t_resp[1], 2/3*t_resp[2])
K5  <- pop_pchaz(
  T = times, 
  lambdaMat1    = m2r(matrix(c(11), nrow = 3, ncol = 1)),   
  lambdaMat2    = m2r(matrix(c( 9, 
                               25,
                               14), nrow = 3)), 
  lambdaProgMat = m2r(matrix(c( 5), nrow = 3,ncol = 1)),
  p = c_resp,
  timezero = TRUE, discrete_approximation = TRUE
)


pp = plot_diagram(B5, K5,
                  A_subgr_labels = c("S1", "S2"),
                  B_subgr_labels = c("Non-switcher","Switcher-S1","Switcher-S2"))
pp

## ----fig.show='hold', fig.height=3.5, fig.width=10, warning=FALSE, out.width='100%'----
plot_shhr(K5, B5)

