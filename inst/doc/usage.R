## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE, fig.show='hold'-----------------------------------------
#  install.packages("nph")
#  
#  # For dev version
#  # install.packages("devtools")
#  devtools::install_github("repo/nph")

## ----echo=FALSE, fig.height=3.5, fig.show='hold', fig.width=10, warning=FALSE, out.width='100%'----
library(nph)
times <- c(0, 5 * 365)   # Time interval boundaries, in days

# Treatment group
t_resp <- 1              # There are no subgroups
B5 <- pop_pchaz(
  T = times,
  lambdaMat1    = m2r(matrix(18, nrow = 1)),
  lambdaMat2    = m2r(matrix(11, nrow = 1)),
  lambdaProgMat = m2r(matrix( 9, nrow = 1)),
  p = t_resp,
  timezero = FALSE, discrete_approximation = TRUE
)

# Control group
c_resp <- 1              # There are no subgroups
K5  <- pop_pchaz(
  T = times,
  lambdaMat1    = m2r(matrix(11, nrow = 1)),
  lambdaMat2    = m2r(matrix( 9, nrow = 1)),
  lambdaProgMat = m2r(matrix( 5, nrow = 1)),
  p = c_resp,
  timezero = TRUE, discrete_approximation = TRUE
)
plot_shhr(A = K5, B = B5, main = "Different disease progression by treatment")

## ----echo=FALSE, fig.height=3.5, fig.show='hold', fig.width=10, warning=FALSE, out.width='100%'----
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

plot_shhr(K5, B5, main = "Different effect by time intervals")

## ----echo=FALSE, fig.height=3.5, fig.show='hold', fig.width=10, warning=FALSE, out.width='100%'----
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

plot_shhr(K5, B5, main = "Presence of subgroups with differential treatment effect")

## ----echo=FALSE, fig.height=4, fig.show='hold', fig.width=7, message=FALSE, warning=FALSE, out.width='75%', paged.print=TRUE----
times <- c(0, 5 * 365)   # Time interval boundaries, in days

# Treatment group
t_resp <- 1              # There are no subgroups
B5 <- pop_pchaz(
  T = times,
  lambdaMat1    = m2r(matrix(18, nrow = 1)),
  lambdaMat2    = m2r(matrix(11, nrow = 1)),
  lambdaProgMat = m2r(matrix( 9, nrow = 1)),
  p = t_resp,
  timezero = FALSE, discrete_approximation = TRUE
)

# Control group
c_resp <- 1              # There are no subgroups
K5  <- pop_pchaz(
  T = times,
  lambdaMat1    = m2r(matrix(11, nrow = 1)),
  lambdaMat2    = m2r(matrix( 9, nrow = 1)),
  lambdaProgMat = m2r(matrix( 5, nrow = 1)),
  p = c_resp,
  timezero = TRUE, discrete_approximation = TRUE
)
plot_diagram(K5, B5, which = "Control")

## ----echo=FALSE, fig.align='center', fig.show='hold', fig.width=12.5, warning=FALSE, out.width='50%'----
knitr::include_graphics("lambdamat.png")

## ---- echo=TRUE, fig.show='hold', fig.width=12.5, warning=FALSE, out.width='100%'----
times <- c(0, 100, 5 * 365)   # Time interval boundaries, in days
t_resp <- c(0.2, 0.8) #Proportion of subgroups
B5 <- pop_pchaz(
  T = times,
  lambdaMat1    = m2r(matrix(c(11, 30,
                               11, 18), byrow = TRUE, nrow = 2)),
  lambdaMat2    = m2r(matrix(c( 9, 20,
                                9, 11), byrow = TRUE, nrow = 2)),
  lambdaProgMat = m2r(matrix(c( 5, 15,
                                5,  9), byrow = TRUE, nrow = 2)),
  p = t_resp, discrete_approximation = TRUE
)

## ---- echo=TRUE, fig.show='hold', fig.width=6, warning=FALSE, out.width='33%'----
plot(B5, main = "Survival function")
plot(B5, fun = "haz", main = "Hazard function")
plot(B5, fun = "cumhaz", main = "Cumulative Hazard function")

## ------------------------------------------------------------------------
times <- c(0, 100, 5 * 365)   # Time interval boundaries, in days
# Treatment group
B5 <- pop_pchaz(T = times,
                   lambdaMat1    = m2r(matrix(c(11, 30,
                                                11, 18), byrow = TRUE, nrow = 2)),
                   lambdaMat2    = m2r(matrix(c( 9, 20,
                                                 9, 11), byrow = TRUE, nrow = 2)),
                   lambdaProgMat = m2r(matrix(c( 5, 15,
                                                 5,  9), byrow = TRUE, nrow = 2)),
                   p = c(0.2, 0.8),#Proportion of subgroups
                discrete_approximation = TRUE 
)
# Control group
K5  <- pop_pchaz(T = times,
                    lambdaMat1    = m2r(matrix(c(11, 11), nrow = 1)),
                    lambdaMat2    = m2r(matrix(c( 9,  9), nrow = 1)),
                    lambdaProgMat = m2r(matrix(c( 5,  5), nrow = 1)),
                    p = 1, discrete_approximation = TRUE
)

## ------------------------------------------------------------------------
# Study set up and Simulation of a data set until interim analysis at 150 events
set.seed(15657)
dat <- sample_fun(K5, B5,
                  r0 = 0.5,                     # Allocation ratio
                  eventEnd = 450,               # maximal number of events
                  lambdaRecr = 300 / 365,       # recruitment rate per day (Poisson assumption)
                  lambdaCens = 0.013 / 365,     # censoring rate per day  (Exponential assumption)
                  maxRecrCalendarTime = 3 * 365,# Maximal duration of recruitment
                  maxCalendar = 4 * 365.25)     # Maximal study duration
head(dat)
tail(dat)

## ------------------------------------------------------------------------
logrank.test(time  = dat$y,
             event = dat$event,
             group = dat$group,
             # alternative = "greater",
             rho   = 1,
             gamma = 0) 
# survival::survdiff(formula = survival::Surv(time  = dat$y, event = dat$event) ~ dat$group)

## ------------------------------------------------------------------------
lrmt = logrank.maxtest(
      time  = dat$y,
      event = dat$event,
      group = dat$group,
      rho   = c(0, 0, 1, 1),
      gamma = c(0, 1, 0, 1)
)
lrmt

## ------------------------------------------------------------------------
lrmt$logrank.test[[1]]
lrmt$logrank.test[[2]]
lrmt$logrank.test[[3]]
lrmt$logrank.test[[4]]

