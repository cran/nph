#R Paket Survival non-proportional hazards

# @title Calculate survival for piecewise constant hazard
#
# @description Calculates hazard, cumulative hazard, survival and distribution function
# 	based on hazards that are constant over pre-specified time-intervals
# 	
# @param Tint vector of boundaries of time intervals (presumably in days) with piecewise constant hazard. The boundaries should be increasing and the first one should
#	be \code{0}, the last one should be larger than the assumed trial duration.
# @param lambda vector of piecewise constant hazards for the intervals speciefied via \code{Tint}. 
# 
# @details Function values are calculated for all integer time points between 0 and the maximum of \code{Tint}.
# @return A list with class \code{mixpch} containing the following components:
# \describe{
#	\item{\code{haz}}{Values of the hazard function.}
#	\item{\code{cumhaz}}{Values of the cumulative hazard function.}
#	\item{\code{S}}{Values of the survival function.}
#	\item{\code{F}}{Values of the distribution function.}
#	\item{\code{t}}{Time points for which the values of the different functions are calculated.}
# }
# @author Robin Ristl, \email{robin.ristl@@meduniwien.ac.at}
# @seealso \code{\link{subpop_hazVFfun}}, \code{\link{pop_hazVFfun}}, \code{\link{plot.mixpch}}
# @examples
# hazVFfun(Tint=c(0,40,100),lambda=c(.02,.05))
#Das heißt ab inklusive Tag 0 gilt lambda=0.02, ab inklusive Tag 40 gilt lambda 0.05
#Wenn das Event am Tag 40 mit dem Hazard 0.05 eintritt, wird es am Tag 41 registriert.
hazVFfun<-function(Tint,lambda) {
  if (length(Tint) != (length(lambda) + 1)) {
    stop("The length of Tint should be equal to the length of lambda + 1")
  }
  t<-seq(0,max(Tint)-1,1)
  haz<-rep(0,length(t))
  cumhaz<-rep(0,length(t))
  for(i in 1:(length(Tint)-1)) {
    haz<-haz+lambda[i]*(Tint[i] <= t & t < Tint[i+1])
    cumhaz<-cumhaz+lambda[i]*(Tint[i+1]-Tint[i])*(t>=Tint[i+1])+ lambda[i]*(t-Tint[i])*(Tint[i] <= t & t < Tint[i+1])
  }
  #haz
  #cumhaz<-c(0,cumsum(haz[-length(t)]))
  #cumhaz-c(0,cumsum(haz[-length(t)]))
  S<-exp(-cumhaz)
  F<- 1-S
  out<-list(haz=haz,cumhaz=cumhaz,S=S,F=F,t=t,
            Tint = Tint, lambda = lambda)
  class(out)<-"mixpch" #damit wir die plot.mixpch Ffunktion verwenden koennen
  out
}

#lambdaTod = lambda1
#lambdaTodnachSwitch = lambda2
#lambdaSwitch = lambdaProg
#nachswithZeit0 = timezero

# @title Calculate survival for piecewise constant hazards with change after random time
#
# @description Calculates hazard, cumulative hazard, survival and distribution function
# 	based on hazards that are constant over pre-specified time-intervals
# 	
# @param Tint vector of boundaries of time intervals (presumably in days) with piecewise constant hazard. The boundaries should be increasing and the first one should
#	be \code{0}, the last one should be larger than the assumed trial duration.
# @param lambda1 vector of piecewise constant hazards before the changeing event happens, for the intervals speciefied via \code{Tint}. 
# @param lambda2 vector of piecewise constant hazards after the changeing event has happened, for the intervals speciefied via \code{Tint}. 
# @param lambdaProg vector of piecewise constant hazards for the changeing event, for the intervals speciefied via \code{Tint}. 
# @param timezero logical, indicating whether after the changeing event the timecount, governing which interval in \code{Tint} and which according value in 
#   \code{lambda2} is used, should restart at zero.
# 
# @details Function values are calculated for all integer time points between 0 and the maximum of \code{Tint}.
# @return A list with class \code{mixpch} containing the following components:
# \describe{
#	\item{\code{haz}}{Values of the hazard function.}
#	\item{\code{cumhaz}}{Values of the cumulative hazard function.}
#	\item{\code{S}}{Values of the survival function.}
#	\item{\code{F}}{Values of the distribution function.}
#	\item{\code{t}}{Time points for which the values of the different functions are calculated.}
# }
# @author Robin Ristl, \email{robin.ristl@@meduniwien.ac.at}
# @seealso \code{\link{hazVFfun}}, \code{\link{pop_hazVFfun}}, \code{\link{plot.mixpch}}
# @examples
# subpop_hazVFfun(Tint=c(0,40,100),lambda1=c(0.2,0.4),lambda2=c(0.1,0.01),
#	lambdaProg=c(0.5,0.4),timezero=FALSE)
# subpop_hazVFfun(Tint=c(0,40,100),lambda1=c(0.2,0.4),lambda2=c(0.1,0.01),
#	lambdaProg=c(0.5,0.4),timezero=TRUE)
subpop_hazVFfun<-function(Tint,lambda1,lambda2,lambdaProg,timezero) {
  t<-seq(0,max(Tint)-1,1)
  
  hVF_Tod<-hazVFfun(Tint,lambda1)
  hVF_TodnachSwitch<-hazVFfun(Tint,lambda2)
  hVF_Switch<-hazVFfun(Tint,lambdaProg)
  #Mischen:
  #gehe alle moeglichen Switchzeitpunkt t durch,
  #fuer jeden Zeitpunkt berechen F bedingt auf diesen Switchzeitpunkt.
  #Mische diese F mit der WS fuer die Switchzeitpunkt
  nt<-length(t)
  
  # Calculate assuming that theprogression is at the end of the day
  Smat<-matrix(NA,nrow=nt,ncol=nt+1) #+1 weil auch die Moeglichkeit gar nicht zu switchen bis max(t) vorkommen kann
  #Zeilen: t, Spalten bedingte Survfunktion (P(Y>t|switch=t_spalte)) fuer switchzeit t_spalte
  prob_t<-c(hVF_Switch$F[1],diff(hVF_Switch$F),hVF_Switch$S[nt]) #WS fuer switch im Intervall (t_i-1,t]
  #if(sum(prob_t)!=1) warning("Prob. ist nicht 1")
  #S=Smat%*%prob_t
  set<-rep(TRUE,nt)	
  #delta<-diff(t)
  for(i in 1:(nt+1)) {
    #Es gibt zwei Versionen: 1) Zeit nach Switch geht fuer hazard_nachswitch(t) einfach weiter:
    if(!timezero) haz<-c(hVF_Tod$haz[!set],hVF_TodnachSwitch$haz[set])
    #2) Zeit nach Switch beginnt fuer hazard_nachswitch wieder bei 0.
    #bzw. nach switch zu Zeit s ist der hazard hazard_nachswitch(t-s):
    if(timezero) {
      #haz<-c(hVF_Tod$haz[!set],hVF_TodnachSwitch$haz[1:(nt-sum(!set))])
      haz<-c(hVF_Tod$haz[!set],hVF_TodnachSwitch$haz[rev(set)])
      #set ist immer von der Struktur FALSE,...,FALSE,TRUE,..,TRUE.
      #rev(set) bewirkt, dass die ersten sum(set) Eintraege von 
      #VF_TodnachSwitch$haz benutzt werden,
      #beachte: c(x,vektor[c(FALSE,...,FALSE)])  ergibt x
    }
    #cumhaz<-c(0,cumsum(haz[-nt]*delta)) #delta ist jetzt immer 1
    cumhaz<-c(0,cumsum(haz[-nt]))
    Smat[,i]<-exp(-cumhaz)
    set[i]<-!set[i] #Switchzeitpunkt um eine Position verschieben
  }
  #plot(t,hVF_TodnachSwitch$haz,col=1,ylim=c(0,3))
  #lines(t,hVF_Tod$haz,col=2)
  
  Smix1<-as.numeric(Smat%*%prob_t)
  
  
  # Calculate assuming that the progression is at the beginning of the day
  Smat<-matrix(NA,nrow=nt,ncol=nt) #+1 weil auch die Moeglichkeit gar nicht zu switchen bis max(t) vorkommen kann
  #Zeilen: t, Spalten bedingte Survfunktion (P(Y>t|switch=t_spalte)) fuer switchzeit t_spalte
  prob_t<-c(diff(hVF_Switch$F),hVF_Switch$S[nt]) #WS fuer switch im Intervall (t_i-1,t]
  #if(sum(prob_t)!=1) warning("Prob. ist nicht 1")
  #S=Smat%*%prob_t
  set<-rep(TRUE,nt)	
  #delta<-diff(t)
  for(i in 1:(nt)) {
    #Es gibt zwei Versionen: 1) Zeit nach Switch geht fuer hazard_nachswitch(t) einfach weiter:
    if(!timezero) haz<-c(hVF_Tod$haz[!set],hVF_TodnachSwitch$haz[set])
    #2) Zeit nach Switch beginnt fuer hazard_nachswitch wieder bei 0.
    #bzw. nach switch zu Zeit s ist der hazard hazard_nachswitch(t-s):
    if(timezero) {
      #haz<-c(hVF_Tod$haz[!set],hVF_TodnachSwitch$haz[1:(nt-sum(!set))])
      haz<-c(hVF_Tod$haz[!set],hVF_TodnachSwitch$haz[rev(set)])
      #set ist immer von der Struktur FALSE,...,FALSE,TRUE,..,TRUE.
      #rev(set) bewirkt, dass die ersten sum(set) Eintraege von 
      #VF_TodnachSwitch$haz benutzt werden,
      #beachte: c(x,vektor[c(FALSE,...,FALSE)])  ergibt x
    }
    #cumhaz<-c(0,cumsum(haz[-nt]*delta)) #delta ist jetzt immer 1
    cumhaz<-c(0,cumsum(haz[-nt]))
    Smat[,i]<-exp(-cumhaz)
    set[i]<-!set[i] #Switchzeitpunkt um eine Position verschieben
  }
  #plot(t,hVF_TodnachSwitch$haz,col=1,ylim=c(0,3))
  #lines(t,hVF_Tod$haz,col=2)
  
  Smix2<-as.numeric(Smat%*%prob_t)
  
  Smix = (Smix1+Smix2)/2
  
  #lines(t,Smix,col=5)
  #lines(t,Smat[,nt+1],col=4)
  
  Fmix<-1-Smix
  cummixhaz<- -log(Smix)
  mixhaz<- c(diff(cummixhaz),NA)
  #fmix<-c(diff(Fmix),NA) 
  
  
  
  
  #cummixhaz-cumsum(c(0,mixhaz[1:(length(mixhaz)-1)]))
  #zwei moegliche numerische Approximationen, besser die, die keine Division erforden, sonst Probleme bei Nennerwerten nahe 0.
  #fmix<-c(diff(Fmix)/delta,NA) 
  #mixhaz<-fmix/Smix
  #mixhaz2<- c(diff(cummixhaz)/delta,NA)
  #fmix2<-mixhaz2*Smix
  
  #plot(mixhaz,mixhaz2)
  #plot(fmix,fmix2)
  #abline(0,1)
  
  out<-list(haz=mixhaz,cumhaz=cummixhaz,S=Smix,F=Fmix,
            t=t, Tint = Tint, 
            lambda1 = lambda1, lambda2 = lambda2, lambdaProg = lambdaProg, timezero = timezero)
  class(out)<-c("mixpch", "subpop")  #mix of piecewise constantant hazard
  out
}


# @title Calculate survival for piecewise constant hazards with change after random time and mixture of subpopulations
#
# @description Calculates hazard, cumulative hazard, survival and distribution function
# 	based on hazards that are constant over pre-specified time-intervals
# 
# @param Tint vector of boundaries of time intervals (presumably in days) with piecewise constant hazard. The boundaries should be increasing and the first one should
#	be \code{0}, the last one should be larger than the assumed trial duration.
# @param lambdaMat1 matrix, each row contains the vector of piecewise constant hazards for one subpopulation before the changeing event happens, for the intervals speciefied via \code{Tint}. 
# @param lambdaMat2 matrix, each row contains the vector piecewise constant hazards for one subpopulation after the changeing event has happened, for the intervals speciefied via \code{Tint}. 
# @param lambdaProgMat matrix, each row contains the vector of piecewise constant hazards for one subpopulation for the changeing event, for the intervals speciefied via \code{Tint}. 
# @param p vector of relative sizes (proportions) of the subpopulations.
# @param timezero logical, indicating whether after the changing event the timecount, governing which interval in \code{Tint} and which according value in 
#   \code{lambda2} is used, should restart at zero. This argument is either of length 1 (applying the same to all subgroups) or the same length as the number of subgroups.
#
# @details Function values are calculated for all integer time points between 0 and the maximum of \code{Tint}.
# @return A list with class \code{mixpch} containing the following components:
# \describe{
#	\item{\code{haz}}{Values of the hazard function.}
#	\item{\code{cumhaz}}{Values of the cumulative hazard function.}
#	\item{\code{S}}{Values of the survival function.}
#	\item{\code{F}}{Values of the distribution function.}
#	\item{\code{t}}{Time points for which the values of the different functions are calculated.}
# }
# @author Robin Ristl, \email{robin.ristl@@meduniwien.ac.at}
# @seealso \code{\link{hazVFfun}}, \code{\link{subpop_hazVFfun}}, \code{\link{plot.mixpch}}
# @examples
# pop_hazVFfun(Tint=c(0,40,100),lambdaMat1=matrix(c(0.2,0.1,0.4,0.1),2,2),
#	lambdaMat2=matrix(c(0.5,0.2,0.6,0.2),2,2),
#	lambdaProg=matrix(c(0.5,0.5,0.4,0.4),2,2),p=c(0.8,0.2),timezero=FALSE)
pop_hazVFfun<-function(Tint,lambdaMat1,lambdaMat2,lambdaProgMat,p,timezero=TRUE) {
  if (!identical(dim(as.matrix(lambdaMat1)), dim(as.matrix(lambdaMat2)), dim(as.matrix(lambdaProgMat)))) {
    stop("Dimensions of all lambda matrices should be the same")
  }
  # Now we are sure lambdaMat1,lambdaMat2,lambdaProgMat have the same dimensions
  # We can compare only one to Tint and p
  if ((length(Tint)-1) != ncol(as.matrix(lambdaMat1))){
    stop("The length of Tint should be equal to the number of column of the lambda matrices + 1")
  }
  if (length(p) != nrow(as.matrix(lambdaMat1))){
    stop("The length of p should be equal to the number of column of the lambda matrices")
  }
  if (length(timezero) == 1){
    timezero = rep(timezero, length(p))
  }
  if (length(timezero) != length(p)){
    stop("timezero should be of length 1 or the same length as p")
  }
  
	t<-seq(0,max(Tint)-1,1)

	#wichtig: damit cumsum fuer kumulativen Hazard genau stimmt:
	nsubpop<-length(p)
	Fmix<-rep(0,length(t))
	#fmix<-rep(0,length(t))
	for(i in 1:nsubpop) {
		hazVF_pop_i<-subpop_hazVFfun(Tint,lambdaMat1[i,],lambdaMat2[i,],lambdaProgMat[i,],timezero=timezero[i])
		Fmix<-Fmix + p[i]*hazVF_pop_i$F
		#fmix<-fmix + p[i]*hazVF_pop_i$f
	}
	Smix<-1-Fmix
	cummixhaz<- -log(Smix)
	mixhaz<- c(diff(cummixhaz),NA)
	#hazmix<-fmix/Smix #Achtung, wenn Smix nahe an 0 ist, wird das numerisch ungenau
	out<-list(haz=mixhaz,cumhaz=cummixhaz,S=Smix,F=Fmix,
	          t=t, Tint = Tint, 
	          lambdaMat1 = lambdaMat1, lambdaMat2 = lambdaMat2, 
	          lambdaProgMat = lambdaProgMat, p = p, timezero = timezero)
	class(out)<-c("mixpch", "pop") #mix of piecewise constantant hazard
	out
}

#' @title Plot mixpch Objects
#'
#' @description Plots survival and other functions stored in \code{mixpch} objects versus time.
#' 
#' @param x an object of class \code{mixpch}.
#' @param fun character string in \code{c("S","F","haz","cumhaz")} indicating which function to plot. Select \code{"S"} for the survival function,
#'	\code{"F"} for the distribution functin, \code{"haz"} for the hazard function or \code{"cumhaz"} for the cumulative hazard function.
#' @param add logical, indicates if the drawing should be added to an existing plot.
#' @param xlab label of the x-axis
#' @param ylab label of the y-axis
#' @param ... further arguments passed to the plotting functions
#' 
#' @author Robin Ristl, \email{robin.ristl@@meduniwien.ac.at}
#' 
#' @seealso \code{\link{pchaz}}, \code{\link{subpop_pchaz}}, \code{\link{pop_pchaz}}
#' 
#' @references
#' Robin Ristl, Nicolas Ballarini, Heiko Götte, Armin Schüler, Martin Posch, Franz König. Delayed treatment effects, treatment switching and 
#' heterogeneous patient populations: How to design and analyze RCTs in oncology. Pharmaceutical statistics. 2021; 20(1):129-145. 
#'
#' @examples
#' A <- pop_pchaz(Tint = c(0, 90, 1500),
#'   lambdaMat1 = matrix(c(0.2, 0.1, 0.4, 0.1), 2, 2) / 365,
#'	 lambdaMat2 = matrix(c(0.5, 0.2, 0.6, 0.2), 2, 2) / 365,
#'	 lambdaProg = matrix(c(0.5, 0.5, 0.4, 0.4), 2, 2) / 365,
#'	 p = c(0.8, 0.2), 
#'	 timezero = FALSE, discrete_approximation = TRUE)
#' plot(A)
#' plot(A, "haz", add = TRUE)
#' 
#' @import graphics
#' @export
plot.mixpch<-function(x,fun=c("S","F","haz","cumhaz"),add=FALSE,ylab=fun,xlab="Time",...) {
	fun<-match.arg(fun,c("S","F","haz","cumhaz"))
	if(!add) plot(x$t,x[[fun]],type="l",ylab=ylab,xlab=xlab,...)
	if(add) lines(x$t,x[[fun]],...)
}


#' @title Draw random survival times from mixpch object.
#'
#' @description Draws independent random survival times from \code{mixpch} objects.
#' 
#' @param n Number of random draws
#' @param x An object of class \code{mixpch}
#' 
#' @return A vector of random survival times.
#' 
#' @details The mixpch object stores the survival function up to some time T. For random times equal or larger T, the value T is returned.
#'
#' @author Robin Ristl, \email{robin.ristl@@meduniwien.ac.at}
#' 
#' @references
#' Robin Ristl, Nicolas Ballarini, Heiko Götte, Armin Schüler, Martin Posch, Franz König. Delayed treatment effects, treatment switching and 
#' heterogeneous patient populations: How to design and analyze RCTs in oncology. Pharmaceutical statistics. 2021; 20(1):129-145. 
#'
#' @seealso \code{\link{rSurv_conditional_fun}}, \code{\link{sample_fun}}, \code{\link{sample_conditional_fun}}
#' 
#' @examples
#' A <- pop_pchaz(Tint = c(0, 90, 1500),
#'   lambdaMat1 = matrix(c(0.2, 0.1, 0.4, 0.1), 2, 2) / 365,
#'	 lambdaMat2 = matrix(c(0.5, 0.2, 0.6, 0.2), 2, 2) / 365,
#'	 lambdaProg = matrix(c(0.5, 0.5, 0.4, 0.4), 2, 2) / 365,
#'	 p = c(0.8, 0.2), 
#'	 timezero = FALSE, discrete_approximation = TRUE)
#' rSurv_fun(n = 10, x = A)
#' 
#' @export
rSurv_fun<-function(n,x) {
	#max(t) muss >= maximale Beobachtungsdauer sein. Denn fuer grosze u wird einfach max(t) erzeugt.
	#x pop_hazVFfun output
	u<-runif(n)
	Y<-rep(0,n)
	for(i in 1:n) {
		#ind<-sum(x$F<=u[i])
		#Wir geben nur diskrete Zeiten auf einen Tag genau aus.
		#Auf diese Art (siehe unten) wird auf den naechsten ganzen Tag aufgerundet.
		#(Die mixpch Objekte kennen dazwischen die Verteilungsfunktion ja auch gar nicht genau.)
		#So kann auch keine Zeit = 0 gezogen werden.
		ind<-sum(x$F<=u[i])
		Y[i]<-x$t[ind]+1
		#ind<-sum(x$F<u[i])+1 #das ist das selbe, aber kann zu Index auszerhalb des Bereichs fuehren
		#Y[i]<-x$t[ind]
	}
	Y
}

#' @title Draw conditional random survival times from mixpch object.
#'
#' @description Draws independent random survival times from \code{mixpch} objects conditional on 
#' 	observed time.
#' @param x An object of class \code{mixpch}
#' @param y A vector of observed right censored times
#' @return A vector of random survival times, conditional on the observed censored times.
#' @details 
#' 	Note that the mixpch object stores the survival function up to some time T. For random times equal or larger T, the value T is returned.
#' @author Robin Ristl, \email{robin.ristl@@meduniwien.ac.at}
#' @seealso \code{\link{rSurv_fun}}, \code{\link{sample_fun}}, \code{\link{sample_conditional_fun}}
#'
#' @references
#' Robin Ristl, Nicolas Ballarini, Heiko Götte, Armin Schüler, Martin Posch, Franz König. Delayed treatment effects, treatment switching and 
#' heterogeneous patient populations: How to design and analyze RCTs in oncology. Pharmaceutical statistics. 2021; 20(1):129-145. 
#'
#' @examples
#' A <- pop_pchaz(Tint = c(0, 90, 1500),
#'   lambdaMat1 = matrix(c(0.2, 0.1, 0.4, 0.1), 2, 2) / 365,
#'	 lambdaMat2 = matrix(c(0.5, 0.2, 0.6, 0.2), 2, 2) / 365,
#'	 lambdaProg = matrix(c(0.5, 0.5, 0.4, 0.4), 2, 2) / 365,
#'	 p = c(0.8, 0.2), 
#'	 timezero = FALSE, discrete_approximation = TRUE)
#' rSurv_conditional_fun(x = A, y = c(10,15,9,2,1))
#' @export
rSurv_conditional_fun<-function(x,y) {
	#max(t) muss >= maximale Beobachtungsdauer sein. Denn fuer grosze u wird einfach max(t) erzeugt.
	#x pop_hazVFfun output
	#y Beobachtungszeit bisher
	n<-length(y)
	u<-runif(n)
	Y<-rep(0,n)
	for(i in 1:n) {
		ind_y<-sum(x$t<=y[i])
		gr_y<-x$t>=y[i]
		F_bedingt<-((x$F - x$F[ind_y] )/ x$S[ind_y])[gr_y]
		#ind<-sum(F_bedingt<=u[i])
		#Wir geben nur diskrete Zeiten auf einen Tag genau aus.
		#Auf diese Art (siehe unten) wird auf den naechsten ganzen Tag aufgerundet.
		#(Die mixpch Objekte kennen dazwischen die Verteilungsfunktion ja auch gar nicht genau.)
		#ind<-sum(F_bedingt<u[i])+1 
		#Y[i]<-(x$t[gr_y])[ind]
		ind<-sum(F_bedingt<=u[i])
		Y[i]<-(x$t[gr_y])[ind]+1
	}
	Y
}

#bedingt auf 0 ergibt das selbe wie unbedingt (muss so sein, auch im R Code ist dann alles gleich)
#bed<-rSurv_conditional_fun(A,rep(0,1000))
#ran<-rSurv_fun(1000,A)
#t.test(bed,ran)
#qqplot(bed,ran)
#table(bed)
#table(ran)
###

#' @title Draw survival times based on study settings
#'
#' @description Simulates data for a randomized controlled survival study.
#' @param A An object of class \code{mixpch}, resembling the survival function in treatment group 0
#' @param B An object of class \code{mixpch}, resembling the survival function in treatment group 1
#' @param r0 Allocation ratio to group 0 (must be a number between 0 and 1)
#' @param eventEnd Number of events, after which the study stops
#' @param lambdaRecr Rate per day for recruiting patients, assuming recruitung follows a Poisson process
#' @param lambdaCens Rate per day for random censoring, assuming censoring times are exponential
#' @param maxRecrCalendarTime Maximal duration of recruitment in days
#' @param maxCalendar Maximal total study duration in days, after which the study stops
#' 
#' @return A data frame with each line representing data for one patient and the following columns:
#' \describe{
#'	\item{\code{group}}{Treatment group}
#'	\item{\code{inclusion}}{Start of observation in terms of calendar time}
#'	\item{\code{y}}{Observed survival/censored time}
#'	\item{\code{yCalendar}}{End of observation in terms of calendar time.}
#'	\item{\code{event}}{logical, \code{TRUE} indicates the observation ended with an event, \code{FALSE} corresponds to censored times}
#'	\item{\code{adminCens}}{logical, \code{True} indicates that the observation is subject to administrative censoring, i.e. the subject was observed until the 
#'		end of the study without an event.}
#'	\item{\code{cumEvents}}{Cumulative number of events over calendar time of end of observation}
#' }
#' The data frame is ordered by \code{yCalendar}
#' @details For simulating the data, patients are allocated randomly to either group (unrestricted randomization).
#' @author Robin Ristl, \email{robin.ristl@@meduniwien.ac.at}
#' @seealso \code{\link{rSurv_fun}}, \code{\link{rSurv_conditional_fun}}, \code{\link{sample_conditional_fun}}
#'
#' @references
#' Robin Ristl, Nicolas Ballarini, Heiko Götte, Armin Schüler, Martin Posch, Franz König. Delayed treatment effects, treatment switching and 
#' heterogeneous patient populations: How to design and analyze RCTs in oncology. Pharmaceutical statistics. 2021; 20(1):129-145. 
#'
#' @examples
#' A <- pop_pchaz(Tint = c(0, 90, 1500),
#'   lambdaMat1 = matrix(c(0.2, 0.1, 0.4, 0.1), 2, 2) / 365,
#'	 lambdaMat2 = matrix(c(0.5, 0.2, 0.6, 0.2), 2, 2) / 365,
#'	 lambdaProg = matrix(c(0.5, 0.5, 0.4, 0.4), 2, 2) / 365,
#'	 p = c(0.8, 0.2), 
#'	 timezero = FALSE, discrete_approximation = TRUE)
#' B <- pop_pchaz(Tint = c(0, 90, 1500),
#'   lambdaMat1 = matrix(c(0.2, 0.1, 0.4, 0.1), 2, 2) / 365,
#'	 lambdaMat2 = matrix(c(0.5, 0.1, 0.6, 0.1), 2, 2) / 365,
#'	 lambdaProg = matrix(c(0.5, 0.5, 0.04, 0.04), 2, 2) / 365,
#'	 p = c(0.8, 0.2), 
#'	 timezero = FALSE, discrete_approximation = TRUE)
#' plot(A)
#' plot(B, add = TRUE)
#' dat <- sample_fun(A, B, r0 = 0.5, eventEnd = 30, lambdaRecr = 0.5,
#'   lambdaCens = 0.25 / 365, maxRecrCalendarTime = 2 * 365,
#'   maxCalendar = 4 * 365)
#' @export
sample_fun<-function(A,B,r0=0.5,eventEnd,lambdaRecr,lambdaCens,maxRecrCalendarTime,maxCalendar) {
  #A, B pop_hazVFfun Objekte
  if(maxCalendar<maxRecrCalendarTime | max(A$t)<maxCalendar | max(B$t)<maxCalendar) stop("Zeiten falsch spezifiziert.")
  
  #alle Raten pro Tag
  #r0<-0.5 #Randomisierungswahrsch fuer Gruppe 0
  
  #Rekrutierungsrate (Poissonprozess)
  #lambdaRecr<-0.5 #pro Tag, alles in Tagen
  #lambdaCens<-0.25/365
  n<-0 #Anzahl Pat
  e<-0 #Anzahl events
  kal<- 0 #kalenderzeit
  
  #eventEnd (frueher e_interim) #Eventzahl, bei der (Interim)analyse ausgefuehrt wird.
  #maxRecrCalendarTime<-2*365
  #maxCalendar<-4*365
  dat<-NULL
  nmax<-rpois(1,lambda=lambdaRecr * maxRecrCalendarTime)
  inclusion<-runif(nmax,0, maxRecrCalendarTime)
  inclusion<-round(inclusion)
  group<-sample(0:1,nmax,TRUE,c(r0,1-r0))
  n0<-sum(group==0)
  n1<-sum(group==1)
  y<-rep(0,nmax)
  y[group==0] <- rSurv_fun(n0,A)
  y[group==1] <- rSurv_fun(n1,B)
  if(lambdaCens>0) zens<-rexp(nmax,lambdaCens) else zens<-rep(Inf,nmax)
  event<-y<zens
  y<-ifelse(event,y,zens)
  y<-round(y)
  yCalendar<-inclusion+y
  #events nach maxCalendar zensieren
  zensEnde<-yCalendar>maxCalendar
  event[zensEnde]<-FALSE
  y[zensEnde]<-maxCalendar-inclusion[zensEnde]
  yCalendar[zensEnde]<-maxCalendar
  
  #suche die Zeit, bis zu der eventEnd events aufgetreten sind (falls das passiert ist).	
  dat<-data.frame(group,inclusion,y,yCalendar,event,adminCens=zensEnde)
  dat<-dat[order(dat$yCalendar),]
  #head(dat)
  dat$cumEvents<-cumsum(dat$event)
  n_events<-sum(dat$event)
  if(n_events>=eventEnd) {
    kalenderEnde<-dat$yCalendar[sum(dat$cumEvents<eventEnd)+1]
    #min(which(dat$cumEvents>=eventEnd))
    kalenderEnde
    #neu zensieren:
    dat<-dat[dat$inclusion<=kalenderEnde,]
    zensEnde2<-dat$yCalendar>kalenderEnde
    dat$adminCens<-zensEnde2
    dat$event[zensEnde2]<-FALSE
    dat$y[zensEnde2]<-kalenderEnde-dat$inclusion[zensEnde2]
    dat$yCalendar[zensEnde2]<-kalenderEnde
    
  } else {
    warning("Geplante Eventzahl nicht erreicht.")
  }
  dat$cumEvents<-cumsum(dat$event)
  #while(e<eventEnd) {
  #	kal<-kal+rexp(1,lambdaRecr)
  #	kal<-round(kal) #Zeiten auf ganze Tage runden:
  #	group<-sample(0:1,1,TRUE,c(p_0,1-p_0))
  #	if(group==0) zeit<-rSurv_fun(1,A) #kontrolle
  #	if(group==1) zeit<-rSurv_fun(1,B) #treatment
  #	zens<-rexp(1,lambdaCens)
  #	y<-min(zeit,zens)	
  #	y<-round(y) #Zeiten auf ganze Tage runden:
  #	yCalendar<-kal+y  #event/zens Zeit in Kalenderzeit
  #	event<-zeit<zens
  #	dat<-rbind(dat,data.frame(group,inclusion=kal,y,yCalendar,event))
  #	e<-sum(dat$yCalendar[dat$event]<=kal)
  #}
  #e
  #sum(dat$event)
  #dim(dat)
  #kal
  #hist(dat$inclusion)
  #dat_obs<-dat
  #head(dat)
  #dat_obs$event[dat_obs$yCalendar>kal]<-FALSE
  #sum(dat_obs$event)
  #dat_obs$yCalendar[dat_obs$yCalendar>kal]<-kal
  #dat_obs$y<-dat_obs$yCalendar-dat_obs$inclusion
  #head(dat_obs)
  #dat_obs
  #class(dat)<-"simsurv" #nicht sinnvoll, weil dadurch dat kein Dataframe mehr ist
  dat
}


##
#' @title Draw conditional survival times based on study settings
#'
#' @description Simulates data for a randomized controlled survival study conditional on observed interim data.
#' 
#' @param dat A data frame with the same structure and column names as the output of \code{\link{sample_fun}}, containing the data to condition on
#' @param A An object of class \code{mixpch}, resembling the survival function in treatment group 0
#' @param B An object of class \code{mixpch}, resembling the survival function in treatment group 1
#' @param r0 Allocation ratio to group 1 (must be a number between 0 and 1)
#' @param eventEnd Number of events, after which the study stops
#' @param lambdaRecr Rate per day for recruiting patients, assuming recruitung follows a Poisson process
#' @param lambdaCens Rate per day for random censoring, assuming censoring times are exponential
#' @param maxRecrCalendarTime Maximal duration of recruitment in days
#' @param maxCalendar Maximal total study duration in days, after which the study stops
#' 
#' @return A data frame with each line representing data for one patient and the following columns:
#' \describe{
#'	\item{\code{group}}{Treatment group}
#'	\item{\code{inclusion}}{Start of observation in terms of calendar time}
#'	\item{\code{y}}{Observed survival/censored time}
#'	\item{\code{yCalendar}}{End of observation in terms of calendar time.}
#'	\item{\code{event}}{logical, \code{TRUE} indicates the observation ended with an event, \code{FALSE} corresponds to censored times}
#'	\item{\code{adminCens}}{logical, \code{True} indicates that the observation is subject to administrative censoring, i.e. the subject was observed until the 
#'		end of the study without an event.}
#'	\item{\code{cumEvents}}{Cumulative number of events over calendar time of end of observation}
#' }
#' The data frame is ordered by \code{yCalendar}
#' @details For simulating the data, patients are allocated randomly to either group (unrestricted randomization).
#' @author Robin Ristl, \email{robin.ristl@@meduniwien.ac.at}
#' @seealso \code{\link{rSurv_fun}}, \code{\link{rSurv_conditional_fun}}, \code{\link{sample_fun}}
#'
#' @references
#' Robin Ristl, Nicolas Ballarini, Heiko Götte, Armin Schüler, Martin Posch, Franz König. Delayed treatment effects, treatment switching and 
#' heterogeneous patient populations: How to design and analyze RCTs in oncology. Pharmaceutical statistics. 2021; 20(1):129-145. 
#'
#'
#' @examples
#' A <- pop_pchaz(Tint = c(0, 90, 1500),
#'   lambdaMat1 = matrix(c(0.2, 0.1, 0.4, 0.1), 2, 2) / 365,
#'	 lambdaMat2 = matrix(c(0.5, 0.2, 0.6, 0.2), 2, 2) / 365,
#'	 lambdaProg = matrix(c(0.5, 0.5, 0.4, 0.4), 2, 2) / 365,
#'	 p = c(0.8, 0.2), 
#'	 timezero = FALSE, discrete_approximation = TRUE)
#' B <- pop_pchaz(Tint = c(0, 90, 1500),
#'   lambdaMat1 = matrix(c(0.2, 0.1, 0.4, 0.1), 2, 2) / 365,
#'	 lambdaMat2 = matrix(c(0.5, 0.1, 0.6, 0.1), 2, 2) / 365,
#'	 lambdaProg = matrix(c(0.5, 0.5, 0.04, 0.04), 2, 2) / 365,
#'	 p = c(0.8, 0.2), 
#'	 timezero = FALSE, discrete_approximation = TRUE)
#' datinterim <- sample_fun(A, B, r0 = 0.5, eventEnd = 30, lambdaRecr = 1,
#'   lambdaCens = 0.25 / 365,
#'	 maxRecrCalendarTime = 3 * 365,
#'	 maxCalendar = 4 * 365)
#' datcond <- sample_conditional_fun(datinterim, A, B, r0 = 0.5, eventEnd = 60,
#'   lambdaRecr = 1, lambdaCens = 0.25 / 365, maxRecrCalendarTime = 3 * 365, 
#'   maxCalendar = 4 * 365)
#' @export
sample_conditional_fun<-function(dat,A,B,r0=0.5,eventEnd,lambdaRecr,lambdaCens,maxRecrCalendarTime,maxCalendar) {
  #Kalenderzeit kal, zu der Stichprobe beendet wurde
  kal<-max(dat$yCalendar)
  if(kal>maxRecrCalendarTime) stop("Rekrutierungsende erreicht.")
  if(maxCalendar<maxRecrCalendarTime | max(A$t)<maxCalendar | max(B$t)<maxCalendar) stop("Zeiten falsch spezifiziert.")
  
  
  #1) Fuer alle Zensierten eventzeiten und zenszeiten ziehen 
  dat_sim<-dat
  #setA<-dat$event==FALSE & dat$group==0
  #setB<-dat$event==FALSE & dat$group==1
  setA<-dat$adminCens==TRUE & dat$group==0
  setB<-dat$adminCens==TRUE & dat$group==1
  
  #zuerst Zensierungszeiten berechnen, denn dafuer brauchen wir das bisherige y
  setAB<-setA | setB
  if(lambdaCens>0) {
    zens_sim<-dat$y[setAB] + rexp(sum(setAB),lambdaCens) #Achtung, das geht nur, wenn die Zensierungszeiten exponentialverteilt sind,
    #weil nur dann die bedingte Verteilung wieder die selbe Exponentialverteilung ist.
  } else {
    zens_sim<-rep(Inf,sum(setAB))
  }
  
  #setA<-dat_sim$event==FALSE & dat$group==0
  #setB<-dat_sim$event==FALSE & dat$group==1
  
  dat_sim$y[setA]<-rSurv_conditional_fun(A,y=dat$y[setA])
  dat_sim$y[setB]<-rSurv_conditional_fun(B,y=dat$y[setB])
  #zensInd<-dat_sim$y[setAB]<zens_sim
  dat_sim$event[setAB]<- dat_sim$y[setAB]<zens_sim
  dat_sim$y[setAB & dat_sim$event==FALSE ]<-zens_sim[!(dat_sim$y[setAB] < zens_sim)]
  
  dat_sim$y[setAB]<-round(dat_sim$y[setAB])
  dat_sim$yCalendar[setAB]<-dat_sim$inclusion[setAB]+dat_sim$y[setAB]
  
  #head(dat_obs,10)
  #head(dat_sim,10)
  
  #2) hypothetische neue Patienten rekrutieren bis zum Ende der Rekrutierungszeit
  #Anzahl ist Poissonverteilt mit Rate lambdaRecr * (maxRecrCalendarTime -kal)
  nneu<-rpois(1,lambda=lambdaRecr * (maxRecrCalendarTime -kal))
  #wann sie rekrutiert werden ist gleichverteilt auf (kal,maxRecrCalendarTime]
  inclusionneu<-runif(nneu,kal,maxRecrCalendarTime)
  inclusionneu<-round(inclusionneu)
  #Gruppe zufaellig zuordnen
  groupneu<-sample(0:1,nneu,TRUE,c(r0,1-r0))
  #Eventzeiten ziehen
  yneu<-rep(NA,nneu)
  yneu[groupneu==0]<-rSurv_fun(sum(groupneu==0),A)
  yneu[groupneu==1]<-rSurv_fun(sum(groupneu==1),B)
  #Zensierung ziehen
  if(lambdaCens>0) zensneu<-rexp(nneu,lambdaCens) else zensneu<-rep(Inf,nneu)
  eventneu<-yneu<zensneu
  yneu[eventneu==FALSE]<-zensneu[eventneu==FALSE]
  yneu<-round(yneu)
  yCalendarneu<-inclusionneu+yneu  #event/zens Zeit in Kalenderzeit
  dat_simneu<-data.frame(group=groupneu,inclusion=inclusionneu,y=yneu,yCalendar=yCalendarneu,event=eventneu,adminCens=FALSE)
  dat_simneu$cumEvents<-cumsum(dat_simneu$event)
  
  
  #3) Daten zusammenfuegen und alle Zeiten die nach Studienende liegen zensieren
  dat_sim_alle<-rbind(dat_sim,dat_simneu)
  #suche die Zeit, bis zu der eventEnd events aufgetreten sind (falls das passiert ist).	
  dat_sim_alle<-dat_sim_alle[order(dat_sim_alle$yCalendar),]
  dat_sim_alle$cumEvents<-cumsum(dat_sim_alle$event)
  n_events<-sum(dat_sim_alle$event)
  if(n_events>=eventEnd) {
    kalenderEnde<-dat_sim_alle$yCalendar[sum(dat_sim_alle$cumEvents<eventEnd)+1]
  } else {
    warning("Geplante Eventzahl nicht erreicht.")
    # in diesem Fall alle Events zensieren, die nach Studienende liegen
    kalenderEnde<-maxCalendar
  }
  dat_sim_alle<-dat_sim_alle[dat_sim_alle$inclusion<=kalenderEnde,]
  zensEnde2<-dat_sim_alle$yCalendar>kalenderEnde
  dat_sim_alle$event[zensEnde2]<-FALSE
  dat_sim_alle$adminCens<-zensEnde2
  dat_sim_alle$y[zensEnde2]<-kalenderEnde-dat_sim_alle$inclusion[zensEnde2]
  dat_sim_alle$yCalendar[zensEnde2]<-kalenderEnde
  dat_sim_alle$cumEvents<-cumsum(dat_sim_alle$event)
  dat_sim_alle
}

###TESTS
#' @title Weighted log-rank test
#'
#' @aliases logrank.test print.wlogrank
#' @description Calculates a weighted log-rank test for the comparison of two groups.
#' 
#' @param time Vector of observed event/censored times
#' @param event logical vector or numeric vector with entries 0 or 1, indicating if an event was observed (TRUE or 1) or the time is censored (FALSE or 0)
#' @param group Vector of group allocations
#' @param alternative Either of \code{"two.sided"},\code{"less"} or \code{"greater"}, specifies if two-sided or respective
#'	one-sided p-values are calculated. In any case the z test statistic of each included weighted log-rank test 
#'	is based on the (weighted) sum of expected minus observed
#'	events in the group corresponding to the first factor level of \code{group}. Hence a small value of the test statistic corresponds to a 
#'	lower (weighted average) hazard rate in the first group.
#' @param rho Parameter to calculate weights in the rho-gamma family
#' @param gamma Parameter to calculate weights in the rho-gamma family
#' @param event_time_weights Optional vector of user defined weights. This weight vector needs to have one entry per event time (not per event, as multiple events may
#'	occur at the same time) and must be sorted by increasing event time.
#' 
#' @return A list with elements:
#' \describe{
#'	\item{\code{D}}{A data frame event numbers, numbers at risk and expected number of events for each event time}
#'	\item{\code{test}}{A data frame containing the z and chi-squared statistic for the one-sided and two-sided test, respectively,
#'		of the null hypothesis of equal hazard functions in both groups and the p-value for the one-sided test.
#'	}
#'	\item{\code{var}}{The estimated variance of the sum of expected minus observed events in the first group.}
#' }
#' 
#' @details 
#' For a given sample, let \eqn{\mathcal{D}} be the set of unique event times.
#' For a time-point \eqn{t \in \mathcal{D}}, let \eqn{n_{t,ctr}} and \eqn{n_{t,trt}} be 
#' the number of patients at risk in the control and treatment group and let
#' \eqn{d_{t,ctr}} and \eqn{d_{t,trt}} be the respective number of events. 
#' The expected number of events in the control group is calculated under the
#' least favorable configuration in \eqn{H_0}, 
#' \eqn{\lambda_{ctr}(t) = \lambda_{trt}(t)}, as \eqn{e_{t,ctr}=(d_{t,ctr}+d_{t,trt}) 
#' \frac{n_{t0}}{n_{t0}+n_{t1}}}. The conditional variance of \eqn{d_{t,ctr}}
#' is calculated from a hypergeometric distribution as 
#' \eqn{var(d_{t,ctr})=\frac{n_{t0} n_{t1} (d_{t0}+d_{t1}) (n_{t0}+n_{t1} - d_{t0} - d_{t1})}{(n_{t0}+n_{t1})^2 (n_{t0}+n_{t1}-1)}}.
#' Further define a weighting function \eqn{w(t)}.
#' The weighted logrank test statistic for a comparison of two groups is 
#' \deqn{z=\sum_{t \in \mathcal{D}} w(t) (d_{t,ctr}-e_{t,ctr}) / \sqrt{\sum_{t \in \mathcal{D}} w(t)^2 var(d_{t,ctr})}}
#' 
#' Under the the least favorable configuration in \eqn{H_0}, 
#' the test statistic is asymptotically standard normally distributed and large 
#' values of \eqn{z} are in favor of the alternative.
#' 
#' The function consider particular weights in the Fleming-Harrington \eqn{\rho-\gamma} 
#' family \eqn{w(t)=\hat S(t-)^\rho (1-\hat S(t-))^\gamma}. 
#' Here, \eqn{\hat{S}(t)=\prod_{s \in \mathcal{D}: s \leq t} 1-\frac{d_{t,ctr}+d_{t,trt}}{n_{t,ctr}+n_{t,trt}}}
#' is the pooled sample Kaplan-Meier estimator.
#' (Note: Prior to package version 2.1, \eqn{S(t)} was used in the definition of \eqn{\rho-\gamma} weights,
#' this was changed to \eqn{S(t-)} with version 2.1.)
#' Weights \eqn{\rho=0, \gamma=0} correspond to the standard logrank test with 
#' constant weights \eqn{w(t)=1}. Choosing \eqn{\rho=0, \gamma=1} puts more weight on 
#' late events, \eqn{\rho=1, \gamma=0} puts more weight on early events and 
#' \eqn{\rho=1, \gamma=1} puts most weight on events at intermediate time points.
#' 
#' @author Robin Ristl, \email{robin.ristl@@meduniwien.ac.at}
#' @seealso \code{\link{logrank.maxtest}}
#' 
#' @references 
#' Robin Ristl, Nicolas Ballarini, Heiko Götte, Armin Schüler, Martin Posch, Franz König. Delayed treatment effects, treatment switching and 
#' heterogeneous patient populations: How to design and analyze RCTs in oncology. Pharmaceutical statistics. 2021; 20(1):129-145. 
#'
#' Thomas R Fleming and David P Harrington. Counting processes and survival analysis. John Wiley & Sons, 2011
#'
#' @examples
#' A <- pop_pchaz(Tint = c(0, 90, 1500),
#'   lambdaMat1 = matrix(c(0.2, 0.1, 0.4, 0.1), 2, 2) / 365,
#'	 lambdaMat2 = matrix(c(0.5, 0.2, 0.6, 0.2), 2, 2) / 365,
#'	 lambdaProg = matrix(c(0.5, 0.5, 0.4, 0.4), 2, 2) / 365,
#'	 p = c(0.8, 0.2), 
#'	 timezero = FALSE, discrete_approximation = TRUE)
#' B <- pop_pchaz(Tint = c(0, 90, 1500),
#'   lambdaMat1 = matrix(c(0.2, 0.1, 0.4, 0.1), 2, 2) / 365,
#'	 lambdaMat2 = matrix(c(0.5, 0.1, 0.6, 0.1), 2, 2) / 365,
#'	 lambdaProg = matrix(c(0.5, 0.5, 0.04, 0.04), 2, 2) / 365,
#'	 p = c(0.8, 0.2), 
#'	 timezero = FALSE, discrete_approximation = TRUE)
#' dat <- sample_fun(A, B, r0 = 0.5, eventEnd = 30,
#'   lambdaRecr = 0.5, lambdaCens = 0.25 / 365,
#'	 maxRecrCalendarTime = 2 * 365,
#'	 maxCalendar = 4 * 365)
#' logrank.test(dat$y, dat$event, dat$group)
#' 
#' 
#' @import stats
#' @export
logrank.test<-function(time,event,group,alternative=c("two.sided", "less", "greater"),
                       rho=0,gamma=0,event_time_weights=NULL) {
  call <- match.call()
  alternative = match.arg(alternative)
  n<-length(time)
  ng<-table(group)
  #ord<-order(time)
  #time<-time[ord]
  #event<-event[ord]
  #group<-factor(group[ord])
  group<-factor(group)
  #if(is.null(weights)) weights<-weights[ord]
  
  Ag<-aggregate(event,by=list(time=time,group=group),FUN=sum,drop=FALSE)
  Ag$x <- ifelse(is.na(Ag$x), 0, Ag$x) #Nico's fix, because in new R version, drop=FALSE creates NA instead of 0 for empty combinations
  #aggregate(event~gruppe*y,data=dat_sim_alle,FUN=sum,drop=FALSE)
  #aggregate(data.frame(event),by=list(time=time,group=group),FUN=sum,drop=FALSE)
  #Die drop=FALSE Einstellung ergibt diese Warnmeldung:
  #Warnmeldung:
  #In `levels<-`(`*tmp*`, value = if (nl == nL) as.character(labels) else paste0(labels,  :
  #  duplicated levels in factors are deprecated
  #tritt nicht auf, wenn y ganzzahlig ist
  
  #str(Ag)
  #Ag
  #cbind(Ag[Ag$group==levels(group)[1],],	Ag[Ag$group==levels(group)[2],])
  tab<-data.frame(time=Ag$time[Ag$group==levels(group)[1]],event1=Ag$x[Ag$group==levels(group)[1]],event2=Ag$x[Ag$group==levels(group)[2]])
  tab$atrisk1<-NA
  tab$atrisk2<-NA
  for(i in 1:dim(tab)[1]) {
    tab$atrisk1[i]<-sum(time[group==levels(group)[1]]>=tab$time[i])
    tab$atrisk2[i]<-sum(time[group==levels(group)[2]]>=tab$time[i])
  }
  
  nz<-dim(tab)[1]
  #tab$atrisk1<-ng[1]-c(0, cumsum(tab$event1)[1:(nz-1)])
  #tab$atrisk2<-ng[2]-c(0, cumsum(tab$event2)[1:(nz-1)])
  tab$atrisk<-tab$atrisk1+tab$atrisk2
  tab$event<-tab$event1+tab$event2
  #tab<-data.frame(time,event,group)#,atrisk,atrisk1,atrisk2,group)
  #tab
  #sf<-survfit(Surv(time=time,event=event)~group)
  #str(sf)
  #as.data.frame(sf)
  #sftab<-data.frame(time=sf$time,n.risk=sf$n.risk,n.event=sf$n.event,gruppe=rep(c(0,1),times=sf$strata))
  #sf0<-sftab[sftab$gruppe==0,]
  #sf1<-sftab[sftab$gruppe==1,]
  #sfall<-merge(sf0,sf1,by="time",all=TRUE)
  
  D<-tab[tab$event>0,]
  D$expected1<-D$event*D$atrisk1/D$atrisk
  D$expected2<-D$event*D$atrisk2/D$atrisk
  #D$expected1_mode<-(D$event+1)*(D$atrisk1+1)/(D$atrisk+2)
  
  #D$observed1<-as.numeric(D$group==levels(group)[1])
  D$diff1<-D$event1-D$expected1
  #D$var<-D$atrisk1*D$event*(D$atrisk-D$atrisk1)*(D$atrisk-D$event)/( D$atrisk^2*(D$atrisk-1) )
  D$var<-D$atrisk1*D$event*D$atrisk2*(D$atrisk-D$event)/( D$atrisk^2*(D$atrisk-1) )
  
  D$S<-	cumprod((D$atrisk-D$event)/D$atrisk) #Kaplan-Meier estimator
  D$Sminus<-c(1,D$S[-length(D$S)])

  #wenn weights bereitgestellt werden:
  if (is.null(event_time_weights)) {
	#rho-gamma weights now use S(t-) instead of S(t)
	D$w <- D$Sminus^rho * (1 - D$Sminus)^gamma 
  } else {
	#user defined weights now need one entry per event time, sorted by event time
	D$w <- event_time_weights/sum(event_time_weights)
	rho<-NA
	gamma<-NA
  }
  #if(is.null(weights)) D$w<-D$S^rho*(1-D$S)^gamma else D$w<-weights[D$time+1] #+1, denn Zeiten starten im weights Vektor mit t=0
  #weights muss Gewichte fuer alle Zeiten t=0,1,... bis zur maximalen Eventzeit enthalten
  #D$w<-D$S^rho*(1-D$S)^gamma
  #D
  Dvoll<-D #aufheben, falls wir das ganze D im output haben wollen
  D<-D[D$atrisk1>0 & D$atrisk2>0,] #sobald es nur mehr eine Gruppe gibt, gibt es keinen Beitrag mehr zur Teststatistik bzw. Varianz wird NaN.
  #Gewichtete Logrank Statistik
  z<-sum(D$w*D$diff1)/sqrt(sum(D$w^2*D$var)) #z-Statistik
  Chisq<-(sum(D$w*D$diff1))^2/sum(D$w^2*D$var)	#Chi-Quadrat Statisitk
  df<-nlevels(group)-1
  if(alternative=="two.sided")	p<-1-pchisq(Chisq,df=df)
  if(alternative=="less")	p<-pnorm(z)
  if(alternative=="greater") p<-1-pnorm(z)
  #S<-Surv(time=time,event=event)
  #lt<-survdiff(S~group) #ok
  #str(lt)
  #colSums(D)
  out = list(D=D,test=data.frame(rho,gamma,z,Chisq,df,p,alternative),var=sum(D$w^2*D$var),
             obs = c(sum(D$event1), sum(D$event2)),
             exp = c(sum(D$expected1), sum(D$expected2)),
             n   =  as.numeric(ng),
             call = call)
  class(out) = "wlogrank"
  out
}


###


#' @title Maximum combination (MaxCombo) log-rank test
#'
#' @aliases logrank.maxtest print.wlogrank_max
#' 
#' @description Calculates a MaxCombo test for the comparison of two groups based on the maximum of test statistics of a set of weighted log-rank tests
#' 
#' @param time Vector of observed event/censored times
#' @param event logical vector or numeric vector with entries 0 or 1, indicating if an event was observed (TRUE or 1) or the time is censored (FALSE or 0)
#' @param group Vector of group allocations
#' @param alternative Either of \code{"two.sided"},\code{"less"} or \code{"greater"}, specifies if two-sided or respective
#'	one-sided p-values are calculated. In any case the z test statistic of each included weighted log-rank test 
#'	is based on the (weighted) sum of expected minus observed
#'	events in the group corresponding to the first factor level of \code{group}. Hence a small value of the test statistic corresponds to a 
#'	lower (weighted average) hazard rate in the first group.
#' @param rho Vector of parameter values rho for a set of weighting functions in the rho-gamma family
#' @param gamma Vector of parameter values gamma for a set of weighting functions in the rho-gamma family
#' @param event_time_weights Optional matrix, each column containing a different weighting vector for the event times.
#' 	These weight vectors need to have one entry per event time (not per event, as multiple events may
#'	occur at the same time) and must be sorted by increasing event time.
#' @param algorithm algorithm for the multivariate normal integration to be used in \code{\link[mvtnorm]{pmvnorm}}.
#' 
#' @return A list with elements:
#' \describe{
#'	\item{\code{pmult}}{The two sided p-value for the null hypothesis of equal hazard functions in both groups, based on the multivariate 
#'	normal approximation for the z-statistics of differently weighted log-rank tests.}
#'	\item{\code{p.Bonf}}{The two sided p-value for the null hypothesis of equal hazard functions in both groups, based on a
#'		Bonferroni multiplicity adjustment for differently weighted log-rank tests.}
#'	\item{\code{tests}}{Data frame with z-statistics and two-sided unadjusted p-values of the individual weighted log-rank tests}
#'	\item{\code{korr}}{Estimated correlation matrix for the z-statistics of the differently weighted log-rank tests.}
#' }
#' 
#' @details 
#' To perform a maximum-type combination test, a set of \eqn{m} different weight 
#' functions \eqn{w_1(t), \dots, w_m(t)} is specified and the correspondingly
#' weighted logrank statistics \eqn{z_1,\dots,z_m} are calculated. The maximum
#' test statistic is \eqn{z_{max}=\max_{i=1,\dots,m} z_i}. If at least one of
#' the selected weight functions results in high power, we may expect a large
#' value of \eqn{z_{max}}. 
#' Under the least favorable configuration in \eqn{H_0}, approximately
#' \eqn{(Z_1,\dots,Z_m)\sim N_m({0},{\Sigma})}. The p-value of the maximum 
#' test, \eqn{P_{H_0}(Z_{max}>z_{max})=1-P(Z_1 \leq z_{max},\dots,Z_m \leq z_{max})},
#' is calculated based on this multivariate normal approximation via numeric integration.
#' The integration is done using \code{\link[mvtnorm]{pmvnorm}}. The default settings in
#' \code{logrank.maxtest} correspond to greater precision than the original default of
#' \code{\link[mvtnorm]{pmvnorm}}. Precision can be set via the argument \code{algorithm}.
#' Lower precision settings may speed up caclulation.
#' 
#' The multivariate normal approach automatically corrects for multiple testing with
#' different weights and does so efficiently since the correlation between the different 
#' tests is incorporated in \eqn{{\Sigma}}. For actual calculations, \eqn{{\Sigma}} is
#'  replaced by an estimate.
#  Note that \eqn{cov (w_i(t)d_{t,ctr},w_j(t)d_{t,ctr})=w_i(t)w_j(t) var(d_{t,ctr})},
#  at least approximately assuming weights are converging in probability to 
#  a non-random function. Thus the \eqn{i,j}-the element of \eqn{{\Sigma}} is estimated as
#  \deqn{\hat{cov}(Z_i,Z_j)=\sum_{t \in \mathcal{D}}  w_i(t)w_j(t)var(d_{t,ctr})/ \sqrt{\sum_{t \in \mathcal{D}}  w_i^2(t) var(d_{t,ctr})   \sum_{t \in \mathcal{D}} w_j^2(t) var(d_{t,ctr})}}
#'  
#' 
#' @author Robin Ristl, \email{robin.ristl@@meduniwien.ac.at}
#' @seealso \code{\link{logrank.test}}
#'
#' @references
#' Robin Ristl, Nicolas Ballarini, Heiko Götte, Armin Schüler, Martin Posch, Franz König. Delayed treatment effects, treatment switching and 
#' heterogeneous patient populations: How to design and analyze RCTs in oncology. Pharmaceutical statistics. 2021; 20(1):129-145. 
#'
#' Pranab Ghosh, Robin Ristl, Franz König, Martin Posch, Christopher Jennison, Heiko Götte, Armin Schüler, Cyrus Mehta. Robust group sequential
#' designs for trials with survival endpoints and delayed response. Biometrical Journal. First published online: 21 December 2021
#'
#' Tarone RE. On the distribution of the maximum of the logrank statistic and the modified wilcoxon statistic. Biometrics. 1981; 37:79-85.
#'
#' Lee S-H. On the versatility of the combination of the weighted log-rank statistics. Comput Stat Data Anal. 2007; 51(12):6557-6564. 
#'
#' Karrison TG et al. Versatile tests for comparing survival curves based on weighted log-rank statistics. Stata J. 2016; 16(3):678-690. 
#'
#' @examples
#' A <- pop_pchaz(Tint = c(0, 90, 1500),
#'   lambdaMat1 = matrix(c(0.2, 0.1, 0.4, 0.1), 2, 2) / 365,
#'	 lambdaMat2 = matrix(c(0.5, 0.2, 0.6, 0.2), 2, 2) / 365,
#'	 lambdaProg = matrix(c(0.5, 0.5, 0.4, 0.4), 2, 2) / 365,
#'	 p = c(0.8, 0.2), 
#'	 timezero = FALSE, discrete_approximation = TRUE)
#' B <- pop_pchaz(Tint = c(0, 90, 1500),
#'   lambdaMat1 = matrix(c(0.2, 0.1, 0.4, 0.1), 2, 2) / 365,
#'	 lambdaMat2 = matrix(c(0.5, 0.1, 0.6, 0.1), 2, 2) / 365,
#'	 lambdaProg = matrix(c(0.5, 0.5, 0.04, 0.04), 2, 2) / 365,
#'	 p = c(0.8, 0.2), 
#'	 timezero = FALSE, discrete_approximation = TRUE)
#' dat <- sample_fun(A, B, r0 = 0.5, eventEnd = 30,
#'   lambdaRecr = 0.5, lambdaCens = 0.25 / 365,
#'	 maxRecrCalendarTime = 2 * 365,
#'	 maxCalendar = 4 * 365)
#' logrank.maxtest(dat$y, dat$event, dat$group)
#' 
#' @importFrom mvtnorm pmvnorm
#' @export
logrank.maxtest<-function(time,event,group,alternative=c("two.sided", "less","greater"),
                          rho=c(0,0,1),gamma=c(0,1,0),event_time_weights=NULL,algorithm = mvtnorm::GenzBretz(maxpts = 50000, abseps = 0.00001, releps = 0)) {  #event_time_weights ist Matrix, pro Spalte ein Gewichtsvektor
  #require(mvtnorm) #das muss fuer das Paket aus der Funktion in die namespace datei verschoben werden.
  call <- match.call()
  alternative = match.arg(alternative)
  if(!is.null(rho))	mrg<-length(rho) else mrg<-0
  if(!is.null(event_time_weights)) {
    if(!is.matrix(event_time_weights)) event_time_weights<-matrix(event_time_weights,ncol=1)
    mw<-dim(event_time_weights)[2]
  } else {
    mw<-0
  }
  m<-mrg+mw
  
  testListe<-vector(mode="list",length=m)
  z<-rep(0,m)
  p<-rep(0,m)
  for(i in 1:m) {
    if(i<=mrg) testListe[[i]]<-logrank.test(time=time,event=event,group=group,alternative=alternative,rho=rho[i],gamma=gamma[i])
    if(i>mrg) testListe[[i]]<-logrank.test(time=time,event=event,group=group,alternative=alternative,event_time_weights=event_time_weights[,i-mrg])
    z[i]<-testListe[[i]]$test$z
    p[i]<-testListe[[i]]$test$p
  }
  V<-diag(1,m)
  for(i in 1:(m-1)) {
    for(j in 2:m) {
      kov<-sum(testListe[[i]]$D$w*testListe[[j]]$D$w*testListe[[i]]$D$var)  #test1$D$var==test2$D$var
      kor<-kov/sqrt(testListe[[i]]$var)/sqrt(testListe[[j]]$var)
      V[i,j]<-V[j,i]<-kor
    }
  }
  #V
  #z
  if(alternative=="two.sided") {
    maxz<-max(abs(z))
    low=-rep(maxz,m)
    up=rep(maxz,m)
    pmult<-1-mvtnorm::pmvnorm(lower=low,upper=up,corr=V,algorithm=algorithm)
  }
  
  if(alternative=="less") {
    minz<-min(z)
    low=rep(minz,m)
    up=rep(Inf,m)
    pmult<-1-mvtnorm::pmvnorm(lower=low,upper=up,corr=V,algorithm=algorithm)
  }
  if(alternative=="greater") {
    maxz<-max(z)
    low=rep(-Inf,m)
    up=rep(maxz,m)
    pmult<-1-mvtnorm::pmvnorm(lower=low,upper=up,corr=V,algorithm=algorithm)
  }
  
  #pmult
  #p
  #2*(1-pnorm(maxz))
  p.Bonf=p.adjust(p,"bonferroni")
  #tests=data.frame(rho,gamma,z,p)
  out = list(pmult=pmult[[1]],p.Bonf=min(p.Bonf),tests=data.frame(Test=1:m,z,p),korr=V,alternative,
             logrank.test = testListe,
             call = call)
  class(out) = "wlogrank_max"
  out
}


