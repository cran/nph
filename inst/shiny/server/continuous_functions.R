## haz_fun is a function that generates a hazard function
# haz_fun <- function(Tint, lambda){
#   hf = approxfun(x = Tint, y = c(lambda),
#                  method = "constant", 
#                  yleft = 0, yright = lambda[length(lambda)], 
#                  f = 0, ties = "ordered")
#   class(hf) = "stepfun"
#   hf
# }
haz_fun <- function(Tint, lambda){
  function(v) {
    mapply(function(vv){
      Tp0 = c(-Inf, Tint)
      lambda0 = c(0, lambda)
      index = max(which(vv >= Tp0))
      lambda0[index]
    }, v)}
}
# hf = haz_fun(Tint = c(0, 40), lambda = c(.02,.05))
# hf
# # plot(hf, verticals = FALSE, pch = 1)
# hf(0)
# hf(0:100)

# Cumulative hazard function -
## cumhaz_fun is a function that generates a Cumulative hazard function
cumhaz_fun <- function(Tint, lambda){
  function(v) {
    mapply(function(vv){
      Tp1 = c(Tint[-1], Inf)
      Tp1 = ifelse(vv > Tp1, Tp1, vv)
      d = (Tp1-Tint)
      sum(d*lambda*(d>0))
    }, v)}
}
# chf = cumhaz_fun(Tint = c(0, 40), lambda = c(.02,.05))
# chf
# plot(chf, pch = 1)
# chf(0)
# chf(100)
# plot(chf, pch = 1, from = 0, to = 100)
# lines(1:100, chf(1:100), pch = 1, col = 2)

# Survival function -
## surv_fun is a function that generates a survival function
surv_fun <- function(Tint, lambda){
  chf = function(v) {
    mapply(function(vv){
      Tp1 = c(Tint[-1], Inf)
      Tp1 = ifelse(vv > Tp1, Tp1, vv)
      d = (Tp1-Tint)
      d*lambda*(d>0)
      exp(-sum(d*lambda*(d>0)))}, v)}
  # class(chf) = "stepfun"
  chf
}

# Cumulative distribution function -
cdf_fun <- function(Tint, lambda){
  chf = function(v) {
    mapply(function(vv){
      Tp1 = c(Tint[-1], Inf)
      Tp1 = ifelse(vv > Tp1, Tp1, vv)
      d = (Tp1-Tint)
      d*lambda*(d>0)
      1-exp(-sum(d*lambda*(d>0)))}, v)}
  # class(chf) = "stepfun"
  chf
}

# distribution function -
pdf_fun <- function(Tint, lambda){
  hf = haz_fun(Tint = Tint, lambda = lambda)
  chf = function(v) {
    mapply(function(vv){
      Tp1 = c(Tint[-1], Inf)
      Tp1 = ifelse(vv > Tp1, Tp1, vv)
      d = (Tp1-Tint)
      d*lambda*(d>0)
      ee = sum(d*lambda*(d>0))
      exp(-ee)*hf(vv)}, # h(t) = f(t) / (1-F(t)), then f(t) = (1-F(t)) * h(t)
      v)}
  # class(chf) = "stepfun"
  chf
}

## Internal function to calculate all functions at once
internal_pchaz <- function(Tint, lambda){
  funs_T   = Tint[-length(Tint)]          # For creating the functions, we dont need the last value of T as we assume the hazard is constant until \infty. 
  haz_f    = haz_fun(Tint = funs_T, lambda = lambda)     # Creates the hazard function
  cumhaz_f = cumhaz_fun(Tint = funs_T, lambda =lambda)  # Creates the cumhazard function
  surv_f   = surv_fun(Tint = funs_T, lambda = lambda)    # Creates the survival function
  cdf_f    = cdf_fun(Tint = funs_T, lambda = lambda)     # Creates the cdf 
  pdf_f    = pdf_fun(Tint = funs_T, lambda = lambda) # Creates the pdf 
  out <- list(haz_f    = haz_f,  cumhaz_f = cumhaz_f, 
              surv_f   = surv_f, cdf_f    = cdf_f, 
              pdf_f    = pdf_f)
  out
}

#' @title Calculate survival for piecewise constant hazard
#'
#' @description Calculates hazard, cumulative hazard, survival and distribution function
#' 	based on hazards that are constant over pre-specified time-intervals. 
#' 	
#' @param Tint vector of length \eqn{k+1}, for the boundaries of \eqn{k} time intervals (presumably in days) with piecewise constant hazard. The boundaries should be increasing and the first one should
#'	be \code{0}, the last one should be larger than the assumed trial duration.
#' @param lambda vector of length \eqn{k} with the piecewise constant hazards for the intervals specified via \code{Tint}. 
#' 
#' @details Given \eqn{k} time intervals \eqn{[t_{j-1},t_j), j=1,\dots,k} with 
#'  \eqn{0 =  t_0 < t_1 \dots < t_k}, the function assume constant hazards \eqn{\lambda_{j}} at each interval.
#'  The resulting hazard function is
#'  \eqn{\lambda(t) =\sum_{j=1}^k \lambda_{j} {1}_{t \in [t_{j-1},t_j)}},
#'  the cumulative hazard function is\\ 
#'  \eqn{\Lambda(t) = \int_0^t \lambda(s) ds =\sum_{j=1}^k \left( (t_j-t_{j-1})\lambda_{j} {1}_{t > t_j} + (t-t_{j-1}) \lambda_{j} {1}_{t \in [t_{j-1},t_j) } \right)} 
#'  and the survival function \eqn{S(t) = e^{-\Lambda(t)}}.
#'  The output includes the functions values calculated for all integer time points
#'  between 0 and the maximum of \code{Tint}.
#'  Additionally, a list with functions is also given to calculate the values at any arbitrary point \eqn{t}.
#' 
#' @return A list with class \code{mixpch} containing the following components:
#' \describe{
#'	\item{\code{haz}}{Values of the hazard function over discrete times t.}
#'	\item{\code{cumhaz}}{Values of the cumulative hazard function over discrete times t.}
#'	\item{\code{S}}{Values of the survival function over discrete times t.}
#'	\item{\code{F}}{Values of the distribution function over discrete times t.}
#'	\item{\code{t}}{Time points for which the values of the different functions are calculated.}
#'	\item{\code{Tint}}{Input vector of boundaries of time intervals.}
#'	\item{\code{lambda}}{Input vector of piecewise constant hazards.}
#'	\item{\code{funs}}{A list with functions to calculate the hazard, cumulative hazard, survival, pdf and cdf over arbitrary continuous times.}
#' }
#' 
#' @author Robin Ristl, \email{robin.ristl@@meduniwien.ac.at}, Nicolas Ballarini
#' @seealso \code{\link{subpop_pchaz}}, \code{\link{pop_pchaz}}, \code{\link{plot.mixpch}}
#' @examples
#' pchaz(Tint = c(0, 40, 100), lambda=c(.02, .05))
#Das heiszt ab inklusive Tag 0 gilt lambda=0.02, ab inklusive Tag 40 gilt lambda 0.05
#Wenn das Event am Tag 40 mit dem Hazard 0.05 eintritt, wird es am Tag 41 registriert.
#' @export
pchaz <- function(Tint, lambda){
  if (length(Tint) != (length(lambda) + 1)){
    stop("The length of Tint should be equal to the length of lambda + 1")
  }
  if (any(lambda < 0)) {
    stop("Only positive values are allowed for the hazard rates")
  }
  if (any(diff(Tint) <= 0)) {
    stop("Tint should be non-decreasing")
  }
  
  call <- match.call()
  t        = seq(0, max(Tint) - 1, 1)  # We define a vector to discretize time 
  funs = internal_pchaz(Tint, lambda)
  
  haz      = funs$haz_f(t)
  cumhaz   = funs$cumhaz_f(t)
  surv     = funs$surv_f(t)
  cdf      = funs$cdf_f(t)
  pdf      = funs$pdf_f(t)
  
  out<-list(haz = haz, cumhaz = cumhaz, S = surv, F = cdf, f = pdf, 
            t=t, Tint = Tint, lambda = lambda, 
            call = call,
            funs = funs)
  class(out)<-"mixpch" #damit wir die plot.mixpch Funktion verwenden koennen
  out
}

# hvfffun = pchaz(Tint = c(0, 40, 100), lambda = c(.02,.05))
# library(nph)
# hvfffun2 = hazVFfun(Tint = c(0, 40, 100), lambda = c(.02,.05))
# plot(hvfffun2, fun = "haz")
# plot(hvfffun,  fun = "haz", add = TRUE, col = "red")
# plot(hvfffun$haz, hvfffun2$haz)
# plot(hvfffun2, fun = "cumhaz")
# plot(hvfffun,  fun = "cumhaz", add = TRUE, col = "red")
# plot(hvfffun2, fun = "S")
# plot(hvfffun,  fun = "S", add = TRUE, col = "red")
# plot(hvfffun2, fun = "F")
# plot(hvfffun,  fun = "F", add = TRUE, col = "red")


## Internal function to calculate all functions at once
internal_subpop_pchaz<-function(Tint, lambda1, lambda2, lambdaProg, timezero = FALSE,
                                int_control) {
  hVF_Tod <- internal_pchaz(Tint, lambda1)
  hVF_TodnachSwitch <- internal_pchaz(Tint, lambda2)
  hVF_Switch <- internal_pchaz(Tint, lambdaProg)
  
  # We first create haz, cumhaz, surv, cdf and pdf functions
  if(!timezero){
    int1 = function(s, t) {
      t1 = hVF_Tod$cumhaz_f(min(s, t))
      t2 = (hVF_TodnachSwitch$cumhaz_f(t) - hVF_TodnachSwitch$cumhaz_f(s)) * (s < t)
      exp(-(t1 + t2)) * hVF_Switch$pdf_f(s)
    }
  } else{
    int1 = function(s, t) {
      t1 = hVF_Tod$cumhaz_f(min(s, t))
      t2 = (hVF_TodnachSwitch$cumhaz_f(t - s)) * (s < t)
      exp(-(t1 + t2)) * hVF_Switch$pdf_f(s)
    }
  }
  int1_vec = Vectorize(int1)
  
  # Mix of Survival functions
  Smix = function(v) {
    mapply(function(x){
      Ss = integrate(f = int1_vec, lower = 0, upper = Inf, t = x, 
                     rel.tol = int_control$rel.tol, abs.tol = int_control$abs.tol)$value
      Ss
    }, v)
  }
  
  # Mix of CDF functions
  Fmix = function(v) {
    mapply(function(x){
      Ss = integrate(f = int1_vec, lower = 0, upper = Inf, t = x, 
                     rel.tol = int_control$rel.tol, abs.tol = int_control$abs.tol)$value
      1 - (Ss)}, v)
  }
  
  # Mix of Cumulative hazards functions
  cummixhaz = function(v) {
    mapply(function(x){
      Ss = integrate(f = int1_vec, lower = 0, upper = Inf, t = x, 
                     rel.tol = int_control$rel.tol, abs.tol = int_control$abs.tol)$value
      -log(Ss)}, v)
  }
  
  if(!timezero){
    int2 = function(s, t) {
      t1 = hVF_Tod$cumhaz_f(min(s, t))
      t2 = (hVF_TodnachSwitch$cumhaz_f(t) - hVF_TodnachSwitch$cumhaz_f(s)) * (s < t)
      exp(-(t1 + t2)) * hVF_Switch$pdf_f(s) * 
        (hVF_Tod$haz_f(t) * (s>t) + hVF_TodnachSwitch$haz_f(t) * (s<t))
    }
  } else {
    int2 = function(s, t) {
      t1 = hVF_Tod$cumhaz_f(min(s, t))
      t2 = (hVF_TodnachSwitch$cumhaz_f(t - s)) * (s < t)
      exp(-(t1 + t2)) * hVF_Switch$pdf_f(s) * 
        (hVF_Tod$haz_f(t) * (s>t) + hVF_TodnachSwitch$haz_f(t-s) * (s<t))
    }
  }
  int2_vec = Vectorize(int2)
  
 # Mix of hazard functions
 hmix = function(v) {
    mapply(function(x){
      invS = (1/Smix(x))
      Ss = integrate(f = int2_vec, lower = 0, upper = Inf, t = x, 
                     rel.tol = int_control$rel.tol, abs.tol = int_control$abs.tol)$value
      invS * Ss
    }, v)
  }
  
  out<-list(haz_f    = hmix, cumhaz_f = cummixhaz, 
            surv_f   = Smix, cdf_f    = Fmix)
  out
}



#' @title Calculate survival for piecewise constant hazards with change after random time
#'
#' @description Calculates hazard, cumulative hazard, survival and distribution function
#' 	based on hazards that are constant over pre-specified time-intervals
#' 	
#' @param Tint vector of length \eqn{k+1}, for the boundaries of \eqn{k} time intervals (presumably in days) with piecewise constant hazard. The boundaries should be increasing and the first one should
#'	be \code{0}, the last one should be larger than the assumed trial duration.
#' @param lambda1 vector of length \eqn{k} for piecewise constant hazards before the changing event happens, for the intervals specified via \code{T}. 
#' @param lambda2 vector of length \eqn{k} for piecewise constant hazards after the changing event has happened, for the intervals specified via \code{T}. 
#' @param lambdaProg vector of length \eqn{k} for piecewise constant hazards for the changing event, for the intervals specified via \code{T}. 
#' @param timezero logical, indicating whether after the changeing event the timecount, governing which interval in \code{Tint} and which according value in 
#'   \code{lambda2} is used, should restart at zero.
#' @param int_control A list with the \code{rel.tol} and \code{abs.tol} paramaters to be passed to the  \code{\link{integrate}} function. 
#' @param discrete_approximation if TRUE, the function uses an approximation based on discretizing the time, instead of integrating. This speeds up the calculations
#'   
#' @details We assume that the time to disease progression \eqn{T_{PD}} is governed
#'  by a separate process with hazard function \eqn{\eta(t)}, 
#'  which does not depend on the hazard function for death \eqn{\lambda(t)}. 
#'  \eqn{\eta(t)}, too, may be modelled as piecewise constant or, for simplicity, 
#'  as constant over time. We define \eqn{\lambda_{prePD}(t)} and \eqn{\lambda_{postPD}(t)} 
#'  as the hazard functions for death before and after disease progression. 
#'  Conditional on \eqn{T_{PD}=s}, the hazard function for death is 
#'  \eqn{\lambda(t|T_{PD}=s)=\lambda_{prePD}(t){I}_{t\leq s}+\lambda_{postPD}(t){I}_{t>s}}.
#'  The conditional survival function is 
#'  \eqn{S(t|T_{PD}=s)=\exp(-\int_0^t \lambda(t|T_{PD}=s)ds)}. 
#'  The unconditional survival function results from integration over all 
#'  possible progression times as \eqn{S(t)=\int_0^t S(t|T_{PD}=s)dP(T_{PD}=s)}. 
#'  The output includes the function values calculated for all integer time points
#'  between 0 and the maximum of \code{Tint}.
#'  Additionally, a list with functions is also given to calculate the values at any arbitrary point \eqn{t}.
#' 
#' @return  A list with class \code{mixpch} containing the following components:
#' \describe{
#'	\item{\code{haz}}{Values of the hazard function.}
#'	\item{\code{cumhaz}}{Values of the cumulative hazard function.}
#'	\item{\code{S}}{Values of the survival function.}
#'	\item{\code{F}}{Values of the distribution function.}
#'	\item{\code{t}}{Time points for which the values of the different functions are calculated.}
#'	\item{\code{Tint}}{Input vector of boundaries of time intervals.}
#'	\item{\code{lambda1}}{Input vector of piecewise constant hazards before the changing event happen.}
#'	\item{\code{lambda2}}{Input vector of piecewise constant hazards after the changing event happen.}
#'	\item{\code{lambdaProg}}{Input vector of piecewise constant hazards for the changing event.}
#'	\item{\code{funs}}{A list with functions to calculate the hazard, cumulative hazard, survival, and cdf over arbitrary continuous times.}
#' }
#' @author Robin Ristl, \email{robin.ristl@@meduniwien.ac.at}, Nicolas Ballarini
#' @seealso \code{\link{pchaz}}, \code{\link{pop_pchaz}}, \code{\link{plot.mixpch}}
#' @examples
#' subpop_pchaz(Tint = c(0, 40, 100), lambda1 = c(0.2, 0.4), lambda2 = c(0.1, 0.01),
#'	lambdaProg = c(0.5, 0.4),timezero = FALSE, discrete_approximation = TRUE)
#' subpop_pchaz(Tint = c(0, 40, 100), lambda1 = c(0.2, 0.4), lambda2 = c(0.1, 0.01),
#'	lambdaProg = c(0.5, 0.4), timezero = TRUE, discrete_approximation = TRUE)
#' @export
subpop_pchaz<-function(Tint, lambda1, lambda2, lambdaProg, timezero = FALSE, 
                       int_control = list(rel.tol = .Machine$double.eps^.4, 
                                          abs.tol = 1e-9),
                       discrete_approximation = FALSE) {
  call <- match.call()
  
  if (any(diff(Tint) <= 0)) {
    stop("Tint should be non-decreasing")
  }
  if (discrete_approximation){
    out = subpop_hazVFfun(Tint, lambda1, lambda2, lambdaProg, timezero)
    return(out)
  }
  
  t <- seq(0, max(Tint) - 1, 1)
  if(!all(c("rel.tol", "abs.tol") %in% names(int_control))){
    int_control = list(rel.tol = .Machine$double.eps^.4, 
                       abs.tol = 1e-9)
  }
  # Handle special cases where lambda progr is 0 or Inf (progression is not considered)
  if(all(lambdaProg == Inf)){
    funs = internal_pchaz(Tint, lambda2)
  } else if(all(lambdaProg == 0)){
    funs = internal_pchaz(Tint, lambda1)
  } else {
    funs = internal_subpop_pchaz(Tint, lambda1, lambda2, lambdaProg, timezero, int_control)
  }
  
  haz = funs$haz_f(t)
  cumhaz = funs$cumhaz_f(t)
  S = funs$surv_f(t)
  F = funs$cdf_f(t)
  
  out <- list(haz = haz, cumhaz = cumhaz, S = S, F = F, 
              t=t, Tint = Tint, 
              lambda1 = lambda1, lambda2 = lambda2, lambdaProg = lambdaProg,  timezero = timezero,
              call = call,
              funs = funs)
  
  class(out) <- c("mixpch", "subpop") #mix of piecewise constantant hazard
  out
}
# a1 = subpop_hazVFfun(Tint=c(0,40,100),
#                      lambda1=c(0.2,0.4),
#                      lambda2=c(0.1,0.01),
#                      lambdaProg=c(0.5,0.4),timezero=FALSE)
# a2 = subpop_pchaz(Tint=c(0,40,100),
#                   lambda1=c(0.2,0.4),
#                   lambda2=c(0.1,0.01),
#                   lambdaProg=c(0.5,0.4),
#                   timezero=FALSE)
# a2$haz
# plot(a1)
# plot(a2, add = TRUE, col = 2)
# str(a2)


#' @title Calculate survival for piecewise constant hazards with change after random time and mixture of subpopulations
#'
#' @description Calculates hazard, cumulative hazard, survival and distribution function
#' 	based on hazards that are constant over pre-specified time-intervals
#' 	
#' @param Tint vector of length \eqn{k+1}, for the boundaries of \eqn{k} time intervals (presumably in days) with piecewise constant hazard. The boundaries should be increasing and the first one should
#'	be \code{0}, the last one should be larger than the assumed trial duration.
#' @param lambdaMat1 matrix of dimension \eqn{m}-by-\eqn{k}, each row contains the vector of piecewise constant hazards for one subpopulation before the changeing event happens, for the intervals speciefied via \code{Tint}. 
#' @param lambdaMat2 matrix of dimension \eqn{m}-by-\eqn{k}, each row contains the vector piecewise constant hazards for one subpopulation after the changeing event has happened, for the intervals speciefied via \code{Tint}. 
#' @param lambdaProgMat matrix of dimension \eqn{m}-by-\eqn{k}, each row contains the vector of piecewise constant hazards for one subpopulation for the changeing event, for the intervals speciefied via \code{Tint}. 
#' @param p vector of length \eqn{m} for relative sizes (proportions) of the subpopulations. They should sum up to 1.
#' @param timezero logical, indicating whether after the changing event the timecount, governing which interval in \code{Tint} and which according value in 
#'   \code{lambda2} is used, should restart at zero. This argument is either of length 1 (applying the same to all subgroups) or the same length as the number of subgroups.
#' @param int_control A list with additional paramaters to be passed to the  \code{\link{integrate}} function. 
#' @param discrete_approximation if TRUE, the function uses an approximation based on discretizing the time, instead of integrating. This speeds up the calculations
#'   
#' @details Given \eqn{m} subgroups with relative sizes \eqn{p_1, \dots, p_m} and 
#' subgroup-specific survival functions \eqn{S{l}(t)}, 
#' the marginal survival function is the mixture \eqn{S(t)=\sum_{l=1}^m p_l S_{l}(t)}.
#' Note that the respective hazard function is not a linear combination of the 
#' subgroup-specific hazard functions. 
#' It may be calculated by the general relation \eqn{\lambda(t)=-\frac{dS(t)}{dt}\frac{1}{S(t)}}.
#' In each subgroup, the hazard is modelled as a piecewise constant hazard, with
#' the possibility to also model disease progression. 
#' Therefore, each row of the hazard rates is used in \code{\link{subpop_pchaz}}. 
#' See \code{\link{pchaz}} and \code{\link{subpop_pchaz}}
#' for more details.
#' The output includes the function values calculated for all integer time points
#'  between 0 and the maximum of \code{Tint}.
#' 
#' Note: this function may be very slow in cases where many time points need to be calculated. If this happens, use
#' \code{discrete_approximation = TRUE}. 
#' 
#' @return A list with class \code{mixpch} containing the following components:
#' \describe{
#'	\item{\code{haz}}{Values of the hazard function.}
#'	\item{\code{cumhaz}}{Values of the cumulative hazard function.}
#'	\item{\code{S}}{Values of the survival function.}
#'	\item{\code{F}}{Values of the distribution function.}
#'	\item{\code{t}}{Time points for which the values of the different functions are calculated.}
#' }
#' @author Robin Ristl, \email{robin.ristl@@meduniwien.ac.at}, Nicolas Ballarini
#' @seealso \code{\link{pchaz}}, \code{\link{subpop_pchaz}}, \code{\link{plot.mixpch}}
#' @examples
#' pop_pchaz(Tint = c(0, 40, 100),
#'   lambdaMat1 = matrix(c(0.2, 0.1, 0.4, 0.1), 2, 2),
#'	 lambdaMat2 = matrix(c(0.5, 0.2, 0.6, 0.2), 2, 2),
#'	 lambdaProg = matrix(c(0.5, 0.5, 0.4, 0.4), 2, 2),
#'	 p = c(0.8, 0.2),
#'	 timezero = FALSE, discrete_approximation = TRUE)
#' @importFrom stats integrate
#' @export
pop_pchaz <- function(Tint, lambdaMat1, lambdaMat2, lambdaProgMat, p, timezero=FALSE, 
                      int_control = list(rel.tol = .Machine$double.eps^.4, 
                                         abs.tol = 1e-9),
                      discrete_approximation = FALSE) {
  call <- match.call()
  if (!identical(dim(as.matrix(lambdaMat1)), dim(as.matrix(lambdaMat2)), dim(as.matrix(lambdaProgMat)))) {
    stop("Dimensions of all lambda matrices should be the same")
  }
  if (any(diff(Tint) <= 0)) {
    stop("Tint should be non-decreasing")
  }
  if (any(c(lambdaMat1, lambdaMat2, lambdaProgMat) < 0)) {
    stop("Only positive values are allowed for the hazard rates")
  }
  # Now we are sure lambdaMat1,lambdaMat2,lambdaProgMat have the same dimensions
  # We can compare only one to Tint and p
  if ((length(Tint)-1) != ncol(as.matrix(lambdaMat1))){
    stop("The length of Tint should be equal to the number of column of the lambda matrices + 1")
  }
  if (length(p) != nrow(as.matrix(lambdaMat1))){
    stop("The length of p should be equal to the number of column of the lambda matrices")
  }  
  if (abs(sum(p)-1)>.Machine$double.eps^.5){
    stop("The prevalences in p should sum up to 1")
  }  
  if (length(timezero) == 1){
    timezero = rep(timezero, length(p))
  }
  if (length(timezero) != length(p)){
    stop("timezero should be of length 1 or the same length as p")
  }
  
  if(discrete_approximation){
    out = pop_hazVFfun(Tint, lambdaMat1, lambdaMat2, lambdaProgMat, p, timezero)
    
    return(out)
  }
  
  if(!all(c("rel.tol", "abs.tol") %in% names(int_control))){
    int_control = list(rel.tol = .Machine$double.eps^.4, 
                       abs.tol = 1e-9)
  }
  t <- seq(0, max(Tint) - 1, 1)
  
  #wichtig: damit cumsum fuer kumulativen Hazard genau stimmt:
  nsubpop <- length(p)
  Fmix <- rep(0,length(t))
  for(i in 1:nsubpop) {
    lambdaProg = lambdaProgMat[i,]
    lambda1 = lambdaMat1[i,]
    lambda2 = lambdaMat2[i,]
    if(all(lambdaProg == Inf)){
      funs = internal_pchaz(Tint, lambda2)
    } else if(all(lambdaProg == 0)){
      funs = internal_pchaz(Tint, lambda1)
    } else {
      funs = internal_subpop_pchaz(Tint, lambda1, lambda2, lambdaProg, timezero[i], int_control)
    }
    Fmix <- Fmix + p[i]*funs$cdf_f(t)
  }
  Smix <- 1 - Fmix
  cummixhaz <- -log(Smix)
  mixhaz <- c(diff(cummixhaz), NA)
  #hazmix<-fmix/Smix #Achtung, wenn Smix nahe an 0 ist, wird das numerisch ungenau
  out <- list(haz = mixhaz, cumhaz = cummixhaz, S = Smix, F = Fmix,
              t=t, Tint = Tint, 
              lambdaMat1 = lambdaMat1, lambdaMat2 = lambdaMat2, 
              lambdaProgMat = lambdaProgMat, p = p, timezero = timezero,
              call = call)
  class(out) <- c("mixpch", "pop") #mix of piecewise constantant hazard
  out
}
# a1 = pop_hazVFfun(Tint=c(0,40,100),lambdaMat1=matrix(c(0.2,0.1,0.4,0.1),2,2),
#              lambdaMat2=matrix(c(0.5,0.2,0.6,0.2),2,2),
#              lambdaProg=matrix(c(0.5,0.5,0.4,0.4),2,2),p=c(0.8,0.2),timezero=FALSE)
# a2 = pop_pchaz(Tint=c(0,40,100),lambdaMat1=matrix(c(0.2,0.1,0.4,0.1),2,2),
#              lambdaMat2=matrix(c(0.5,0.2,0.6,0.2),2,2),
#              lambdaProg=matrix(c(0.5,0.5,0.4,0.4),2,2),p=c(0.8,0.2),timezero=FALSE)
# 
# 
# plot(a1)
# plot(a1, add = TRUE, col = 2)
# plot(a1, fun = "cumhaz")
# plot(a1, fun = "cumhaz", add = TRUE, col = 2)

#' @export
print.mixpch = function(x, ...){
  m1 = "pchaz object: piecewise constant hazard"
  if(class(x)[2] %in% "subpop") m1 = sprintf("%s, with change after random time", m1)
  if(class(x)[2] %in% "pop")    m1 = sprintf("%s, with change after random time and mixture of subpopulations", m1)
  
  tt = c("Time intervals: ",
         paste0(rep("|", length(x$Tint)), collapse = "-----"), 
         "\n           ",
         paste0(sprintf("%6.f", x$Tint), collapse = ""), "\n\n")
  
  ss = ""
  if(class(x)[2] %in% "pop"){
    ss = paste0("Number of subgroups: ", length(x$p))
    if(length(x$p)>1) ss = sprintf("%s, with prevalences: %s", 
                                   ss,
                                   paste0(x$p, collapse = " - "))
  }
  
  overall_median = which.min(abs(x$F - 0.5))+1
  if (overall_median == max(x$Tint)) overall_median = NA
  mm = c("Median survival time: ", overall_median)
  
  cat(m1, "\n\n")
  cat(tt)
  cat(ss, "\n")
  cat(mm, "\n")
}
