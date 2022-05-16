#package Funktionen
#Basis: 
#"Funktionen_nph_params_1Dez2021_1-20_mit_perturb_scen5_neu.R"
#"perturb_18Okt2021_1-1.R"

####Funktionen simultane Inferenz
indfun<-function(t,km) sum(km$time<=t)

est_fun<-function(t,ev) {
  KMest<-survfit(Surv(time=t,event=ev)~1)
  NAest<-cumsum(KMest$n.event/KMest$n.risk)
  est<-exp(-NAest)
  data.frame(t=KMest$time,ev=KMest$n.event,NAsurv=est,KMsurv=KMest$surv,obs=1,atrisk=KMest$n.risk)
}

fillo<-function(x) {
  for(i in 1:length(x)) {
    if(is.na(x[i])) {
      if(i==1) x[i]<-1 else x[i]<-x[i-1]
    }
  }
  x
}

fillo_at_risk<-function(x,ev) {
  for(i in length(x):1) {
    if(is.na(x[i])) {
     	if(i==length(x)) x[i]<-0 else x[i]<-x[i+1]
    }
  }
  x
}

closest<-function(x,vek) {
  dist<-abs(vek-x)
  which(dist==min(dist))[1]
}

getHval<-function(i,est,param,est01,gruppe,t0,t1,ev0,ev1,haz_method) {
	if(param$type[i]=="S") {
		pos<-sum(est$t<=param$par[i])
		val<-est$NAsurv[pos]
		H<-ifelse(est$t<=param$par[i],-val,0)
	}
	if(param$type[i]=="logS") {
		pos<-sum(est$t<=param$par[i])
		val<-log(est$NAsurv[pos])
		H<-ifelse(est$t<=param$par[i],-1,0)
	}
	if(param$type[i]=="cloglogS" | param$type[i]=="cumhaz") {
		pos<-sum(est$t<=param$par[i])
		kumhaz<- -log(est$NAsurv[pos])
		val<-log(kumhaz)
		H<-ifelse(est$t<=param$par[i],1/kumhaz,0)
	}

	if(param$type[i]%in%c("Q","logQ")) {
		pos<-sum(est$NAsurv>(1-param$par[i]))+1
		#pos<-sum(est$NAsurv>param$par[i])+1
		val<-est$t[pos]
		if(haz_method=="local") {
			#lokaler Hazard	
			#simple Methode:
			#m<-dim(est)[1]
			est$personyears<-diff(c(0,est$t))*est$atrisk
		
			set<-est$ev>0
			ind<-(1:dim(est)[1])[set]
			ind_val<-sum((est$ev>0)[1:pos])
			m<-sum(est$ev)
			delta<-ceiling(sqrt(m))*2
			low<-ind[max(c(1,ind_val-delta))]
			up<-ind[min(c(dim(est)[1],ind_val+delta))]
			#persontimes<-sum(diff(est$t[low:up])*est$atrisk[(low+1):up])
			#events<-sum(est$ev[(low+1):up])
			haz<-sum(est$ev[low:up])/sum(est$personyears[low:up])

		}
		if(haz_method=="muhaz") {
			#muhaz
			if(gruppe==0) mhaz<-muhaz(times=t0,delta=ev0) #,bw.method="global",bw.grid=BW)
	      	if(gruppe==1) mhaz<-muhaz(times=t1,delta=ev1)
			pos_q<-closest(val,vek=mhaz$est.grid)
      		haz<-mhaz$haz.est[pos_q]
		}

		if(param$type[i]=="Q") {	
			H<-ifelse(est$t<=val,-1/haz,0)
		}
		if(param$type[i]=="logQ") {
			H<-ifelse(est$t<=val,-1/haz/val,0)
			val<-log(val)
		}
	}
	if(param$type[i]=="score") {
		if(is.na(param$par[i])) par_score<-Inf else par_score<-param$par[i]
		t_leq_pa<-as.numeric(est01$t<=par_score)
		if(gruppe==0) group_ind<-est01$ev0>=1
		if(gruppe==1) group_ind<-est01$ev1>=1
		#if(gruppe==0) H<-est01$atrisk1/(est01$atrisk0+est01$atrisk1)
		#if(gruppe==1) H<-(1-est01$atrisk1/(est01$atrisk0+est01$atrisk1))
		#if(gruppe==1) H<-(est01$atrisk0/(est01$atrisk0+est01$atrisk1))
		#Bei Gruppe 0 umgekehrtes Vorzeichen, weil val_1 - val_0 gerechnet wird
		
		#Wenn wir noch für das allgemeine Schema dN_i/Y_i^2 mit Y_i multiplizieren, wird es einheitlich:
		H<-(est01$atrisk0*est01$atrisk1/(est01$atrisk0+est01$atrisk1)*t_leq_pa)[group_ind] #*dN_gruppe
		#H0<-H
		#H1<-H
		#sum(H0/est0$atrisk)^2+sum(H1/est1$atrisk)^2
		H/est01$atrisk0[group_ind]
		H0<-(0-est01$atrisk1/(est01$atrisk0+est01$atrisk1))
		
		H1<-(1-est01$atrisk1/(est01$atrisk0+est01$atrisk1))
		U0<- H0*est01$ev0*t_leq_pa
		U1<- H1*est01$ev1*t_leq_pa
		if(gruppe==0) val = -sum(U0)
		if(gruppe==1) val = sum(U1)		
		#score<-sum(U1)+sum(U0) #aber in kov2fun wird val_1-val_0 gerechnet
		N<-length(t0)+length(t1)
		val<-val/N #mean score wird als wahrer Paremeter benutzt.
		#daher H auch durch N dividieren
		H<-H/N
		#meanU<-stat/N#(sum(ev0)+sum(ev1))
		#v<-sum(H0^2*est01$ev0)+sum(H1^2*est01$ev1)
		#stat/sqrt(v)
		#stat^2/v
		#n_ev<-sum(est01$ev1)+sum(est01$ev0)


	}

	if(param$type[i]=="HR") {

		gr<-rep(0:1,times=c(length(t0),length(t1)))
		if(is.na(param$par[i])) par_cox<-Inf else par_cox<-param$par[i]
		t_leq_pa<-as.numeric(est01$t<=par_cox)

		time_cox<-c(t0,t1)
		event_cox<-c(ev0,ev1)
		event_cox[time_cox>par_cox]<-FALSE

		cox<-coxph(Surv(time=time_cox,event=event_cox)~gr,ties="breslow")
		logHR_est<-as.numeric(coef(cox))
		hr<-exp(logHR_est)
		z<-est01$atrisk1*hr/(est01$atrisk0+est01$atrisk1*hr)
		#H0<- -z
		#H1<- 1-z
		#U0<- H0*est01$ev0*t_leq_pa
		#U1<- H1*est01$ev1*t_leq_pa
		#sum(c(U0,U1)) #=0
		if(gruppe==0) {
			group_ind<-est01$ev0>=1
			H<- (z*est01$atrisk0*t_leq_pa)[group_ind] #ohne minus Vorzeichen, weil val1-val0 gerechnet wird
							#mal atrisk0, weil im allgemeinen Schema mal dN0/atrisk0 dazukommt
			val<-0
		}

		if(gruppe==1) {
			group_ind<-est01$ev1>=1
			H<- ((1-z)*est01$atrisk1*t_leq_pa)[group_ind] #
							#mal atrisk0, weil im allgemeinen Schema mal dN1/atrisk1 dazukommt
			val<-logHR_est
		}
		invHesse<-as.numeric(vcov(cox))
		H<-H*invHesse  #normale Taylorentwicklung, um Varianz von logHR_est zu bekommen.

	}
	if(param$type[i]=="avgHR") {
		pos<-sum(est$t<=param$par[i])
		weights01_predictable<-(c(1,est01$NAsurv0)*c(1,est01$NAsurv1))[-(1+dim(est01)[1])]
		weights01<-weights01_predictable
  		#weights0<-weights01[!is.na(est01$obs0)]
  		#weights1<-weights01[!is.na(est01$obs1)]
 		if(gruppe==0) weights<-weights01[est01$ev0>=1]
 		if(gruppe==1) weights<-weights01[est01$ev1>=1]

		avgHR_contr<-sum((weights*est$ev/est$atrisk)[est$t<=param$par[i]])
		H<- ifelse(est$t<=param$par[i],weights/avgHR_contr,0)
		val<-log(avgHR_contr)
	}
	if(param$type[i]=="RMST") {
		pos<-sum(est$t<=param$par[i])
		times<-est$t[est$t<=param$par[i]]
		
		S<-est$NAsurv[est$t<=param$par[i]]
		#if(times[1]>0) { kann man immer machen, erstes dt kann ja auch 0 sein.
			times<-c(0,times)
			S<-c(1,S)
			#plot(times,S,type="s")
		#}
		m<-length(times)
		#times[m]<-param$par[i]
		#if(times[m]<param$par[i]) { #kann man immer machen, letztes dt kann ja auch 0 sein.
			times<-c(times,param$par[i])
			#S<-c(S,S[m])
		#}
		dt<-diff(times)
		Sint<-S #[-length(S)]
		val<-sum(Sint*dt)
	
		W<-rev(cumsum(rev(Sint[-1]*dt[-1])))
		W<-c(W,rep(0,length(est$ev)-length(W)))
		H<- -W  #est$ev_sqrt/est$atrisk
	}
	list(H=H,val=val)
}


dec2bin<-function (x, digits = 8) {
    b <- rep(0, digits)
    for (k in 1:digits) {
        b[k] <- x%%2
        x <- (x - b[k])/2
    }
    b
}

clo_matrix_fun<-function(m) {
	n_tests_ct<-2^m-1
	clo_matrix<-matrix(NA,nrow=n_tests_ct,ncol=m)
	for(i in 1:n_tests_ct) {
		clo_matrix[i,]<-as.logical(dec2bin(i,digits=m))
	}
	clo_matrix
}

#' @title Simultaneous Inference For Parameters Quantifying Differences Between Two Survival Functions
#'
#' @description Hypothesis tests with parametric multiple testing adjustment and simultaneous confidence intervals
#' for a set of parameters, which quantify differences between two survival functions. Eligible parameters are differences in survival probabilities, log survival probabilities,
#' complementary log log (cloglog) transformed survival probabilities, quantiles of the survival functions,
#' log transformed quantiles, restricted mean survival times, as well as an average hazard ratio, the Cox model score statistic
#' (logrank statistic), and the Cox-model hazard ratio.
#'
#' @param time vector of observed event/censored times.
#' @param event Vector with entries 0 or 1 (or FALSE/TRUE) indicating if an event was observed (1) or the time is censored (0).
#' @param group group indicator, must be a vector with entries 0 or 1 indicating the allocation of a subject to one of two groups.
#'	Group 0 is regarded as reference group when calculating parameters.
#' @param data an optional data frame containing the time, event and group data.
#' @param param_type character vector defining the set of parameters that should be analysed. Possible entries are "S","logS","cloglogS",
#' 	"Q","logQ","RMST","avgHR","score" and "HR", representing differences in survival probabilities, log survival probabilities,
#'	complementary log log (cloglog) transformed survival probabilities, quantiles of the survival functions,
#'	log transformed quantiles, restricted mean survival times, as well as an average hazard ratio, the Cox model score statistic (logrank statistic),
#'	and the Cox-model hazard ratio.
#' @param param_par numeric vector which contains the time points at which the requested parameters are evaluated (e.g. x-year survival or RMST after x-years),
#'	or, in case of analysing quantiles, the according probability. May be \code{NA} for parameter types "RMST","avgHR","score" or "HR". In this case,
#'	the minimum of the largest event times of the two groups is used. Also, times greater than this minimum are replaced by this minumum for "RMST","avgHR","score" or "HR".
#' @param param_alternative optional character vector with entries "less" or "greater", defining the alternative for each parameter.
#'	Only required if one-sided tests or one-sided confidence intervals are requested.
#'	Note that group 0 is regarded as reference group when calculating parameters and therefore whether "greater" or "less" corresponds 
#'	to a benefit may depend on the type of parameter. In general, to show
#'	larger survival in group 1 compared to group 0, alternatives would be "greater" for parameter types "S", "logS", "Q", "logQ" and "RMST"
#'	and would be "less" for parameters types "cloglogS",
#'	"avgHR","HR", and "score". 
#'	(The score test is defined here such that alternative "less" corresponds to smaller hazard (and better survival) in group 1 compared to group 0.)
#' @param lvl Confidence level. Applies to, both, unadjusted and multiplicity adjusted (simultaneous) confidence intervals.
#' @param closed_test logical indicating whether p-values should be adjusted using a closed testing procedure. Default is \code{FALSE}, and in this case
#'	p-values will be adjusted by a single step procedure.
#'	With \eqn{k} hypotheses this involves the computation of \eqn{2^k} tests, which may require considerable computation time.
#' @param alternative_test character with possible values "tow.sided" (default) or "one-sided". Specifies whether hypothesis tests should be two-sided or one-sided. In the #'	latter case, \code{param_alternative} must be defined.
#' @param alternative_CI character with possible values "tow.sided" (default) or "one-sided". Specifies whether confidence intervals should be two-sided or one-sided.
#'	In the latter case, \code{param_alternative} must be defined.
#' @param haz_method character with possible values "local" or "muhaz". Specifies whether local hazard should be calculated under a local constant hazard assumption (default) #'	or using the function \code{\link[muhaz]{muhaz}} from the muhaz package. Only relevant when median or log(median) survival
#'	times are analysed.
#' @param rhs right-hand side vector of null hypotheses. Refers to log-scaled difference for ratios. Default is to consider for all null hypothesis a difference of 0.
#' @param perturb logical, indicating whether the perturbation based estiamte should be used instead of the asymptotic estimate to calculate the covariance matrix.
#'	Defaults to \code{FALSE}.
#' @param Kpert The number of perturbation samples to be used with the perturbation approach for covariance estimation.
#'
#' @return A list of class \code{nphparams} with elements:
#' \describe{
#'	\item{\code{est}}{Estimated differences (at log-scale in case of ratios).}
#'	\item{\code{V}}{Estimated covariance matrix of differences.}
#'	\item{\code{tab}}{A data frame with analysis results. Contains the parameter type (Parameter) and settings (Time_or_which_quantile), the estimated difference (Estimate),
#'		its standard error (SE), unadjusted confidence interval lower and upper bounds (lwr_unadjusted, upr_unadjusted), unadjusted p-values (p_unadj),
#'		mulitplicity adjusted confidence interval lower and upper bounds (lwr_adjusted, upr_adjusted), single-step multiplcity adjusted p-values (p_adj),
#'		closed-test adjusted p-values, if requested (p_adjusted_closed_test) and for comparison Bonferroni-Holm adjusted p-values (p_Holm).}
#'	\item{\code{param}}{The used parameter settings. If \code{param_par} was \code{NA} for "HR","avgHR" or "RMST", it is replaced by \code{minmaxt} here.}
#'	\item{\code{paramin}}{The parameter settings as provided to the function. The only difference to \code{param} is in \code{param_par}, as \code{NA} is not replaced here.}
#'	\item{\code{dat0}}{A data frame with information on all observed events in group 0. Contains time (t), number of events (ev), Nelson-Aalen estimate (NAsurv) 
#'		and Kaplan-Meier estimate (KMsurv) of survival, and the number at risk (atrisk).}
#'	\item{\code{dat1}}{A data frame with information on all observed events in group 1. Contains time (t), number of events (ev), Nelson-Aalen estimate (NAsurv) 
#'		and Kaplan-Meier estimate (KMsurv) of survival, and the number at risk (atrisk).}
#'	\item{\code{minmaxt}}{Minimum of the largest event times of the two groups.}
#'	\item{\code{est0}}{Estimated parameter values in group 0.}
#'	\item{\code{est1}}{Estimated parameter values in group 1.}
#'	\item{\code{V0}}{Estimated covariance matrix of parameter estimates in group 0.}
#'	\item{\code{V1}}{Estimated covariance matrix of parameter estimates in group 1.}
#' }
#' 
#' @author Robin Ristl, \email{robin.ristl@@meduniwien.ac.at}
#' @seealso \code{\link{print.nphparams}}, \code{\link{plot.nphparams}}
#' 
#' @examples
#' data(pembro)
#' set1<-nphparams(time=time, event=event, group=group,data=pembro,
#'	param_type=c("score","S"),
#'	param_par=c(3.5,2),
#'	param_alternative=c("less","greater"),
#'	closed_test=TRUE,alternative_test="one.sided")
#' print(set1)
#' plot(set1,trt_name="Pembrolizumab",ctr_name="Cetuximab")
#'
#' set2<-nphparams(time=time, event=event, group=group, data=pembro,  
#'	param_type=c("S","S","S","Q","RMST"),
#'	param_par=c(0.5,1,2,0.5,3.5))
#' print(set2)
#' plot(set2,showlines=TRUE,show_rmst_diff=TRUE)
#'
#' #Create a summary table for set2, showing parameter estimates for each group and the
#' #estimated differences between groups. Also show unadjusted and multiplicity adjusted
#' #confidence intervals using the multivariate normal method and, for comparison,
#' #Bonferroni adjusted confidence intervals:
#'
#' set2Bonf<-nphparams(time=time, event=event, group=group, data=pembro,  
#'	param_type=c("S","S","S","Q","RMST"),
#'	param_par=c(0.5,1,2,0.5,3.5),
#'	lvl=1-0.05/5)
#' KI_paste<-function(x,r) {
#'	x<-round(x,r)
#'	paste("[",x[,1],", ",x[,2],"]",sep="")
#' }
#' r<-3
#' tab<-data.frame(
#'	Parameter=paste(set2$tab[,1],set2$tab[,2]),
#'	Pembrolizumab=round(set2$est1,r),
#'	Cetuximab=round(set2$est0,r),
#'	Difference=round(set2$tab$Estimate,r),
#'	CI_undadj=KI_paste(set2$tab[,5:6],r),
#'	CI_adj=KI_paste(set2$tab[,8:9],r),
#'	CI_Bonf=KI_paste(set2Bonf$tab[,c(5:6)],r))
#' tab
#'
#' @export
nphparams<-function(time,event,group,data=parent.frame(),param_type,param_par=NA,param_alternative=NA,lvl=0.95,closed_test=FALSE,alternative_test="two.sided",alternative_CI="two.sided",haz_method="local",rhs=0,perturb=FALSE,Kpert=500) {
	call <- match.call()
	if (typeof(data) == "environment") {
		time<-time
		event<-event
		group<-group
	} else {
		#time
   		if(length(call$time) == 1) {
	            time.col <- which(colnames(data) == call$time)
            	if(length(time.col)>0) {
				time<-data[,time.col]
            	} else {
				time<-eval(call$time, envir = parent.frame())
			}
            }
		#event
		if(length(call$event) == 1) {
	            event.col <- which(colnames(data) == call$event)
            	if(length(event.col)>0) {
				event<-data[,event.col]
            	} else {
				event<-eval(call$event, envir = parent.frame())
			}
            }
		#group
		if(length(call$group) == 1) {
	            group.col <- which(colnames(data) == call$group)
            	if(length(group.col)>0) {
				group<-data[,group.col]
            	} else {
				group<-eval(call$group, envir = parent.frame())
			}
            }
	}

	if(length(time)!=length(event) | length(time)!=length(group)) stop("Arguments time, event and group must be vectors of equal legnth.")
	if(!all(group%in%c(0,1))) stop("Argument group must be a vector with entries 0 or 1.")
	if(all(group==0) | all(group==1)) stop("Only one group defined.")

	gr0_ind<-group==0
	gr1_ind<-group==1
	t0<-time[gr0_ind]
	t1<-time[gr1_ind]
	ev0<-event[gr0_ind]
	ev1<-event[gr1_ind]

	#Put together information on the requested paramters
	#if(is.null(param)) {
		param=data.frame(
      	  type=param_type,
	        par=param_par,
      	  alternative=param_alternative
		)
	#}
	allowed_parameters<-c("S","logS","cloglogS","Q","logQ","RMST","avgHR","score","HR")
	if(!all(param$type%in%allowed_parameters)) stop(paste("All requested parameters must be in c(",paste(allowed_parameters,collapse=" , "),").",sep=""))
	
	summary_parameters<-c("RMST","avgHR","score","HR")
	SorQ_parameters<-c("S","logS","cloglogS","Q","logQ")
	if(any(param$type%in%SorQ_parameters)) {
		param_need_par<-param[param$type%in%SorQ_parameters,]
		if(any(is.na(param_need_par$par))) stop(paste("Time point or quantile probability must be provided for parameters of types",paste(SorQ_parameters,collapse=", "),"."))
		
	}

	#ALTERNATIVE="two.sided"
	if(length(rhs)==1) rhs<-rep(rhs,dim(param)[1])
	paramin<-param
	#param$par[param$type=="HR"]<-NA
	minmaxt<-min(max(t0),max(t1))
	param$par[(param$type=="avgHR"|param$type=="RMST"|param$type=="HR"|param$type=="score") & (is.na(param$par)|param$par>minmaxt)]<-minmaxt
	
	est0<-est_fun(t=t0,ev=ev0)
  	est1<-est_fun(t=t1,ev=ev1)
	if(!perturb) {
		
	  	names(est0)<-paste(names(est0),"0",sep="")
  		names(est1)<-paste(names(est1),"1",sep="")
  		est01<-merge(est0,est1,by.x="t0",by.y="t1",all=TRUE)
  		#if(is.na(est01$atrisk0[1])) est01$atrisk0[1]<-max(est01$atrisk0,na.rm=TRUE)
  		#if(is.na(est01$atrisk1[1])) est01$atrisk1[1]<-max(est01$atrisk1,na.rm=TRUE)

  		est01$atrisk0<-fillo_at_risk(est01$atrisk0,ev=est01$ev0)
 		est01$atrisk1<-fillo_at_risk(est01$atrisk1,ev=est01$ev1)
 		est01$NAsurv0<-fillo(est01$NAsurv0)
 		est01$KMsurv0<-fillo(est01$KMsurv0)
 		est01$NAsurv1<-fillo(est01$NAsurv1)
 		est01$KMsurv1<-fillo(est01$KMsurv1)
 		est01$ev0[is.na(est01$ev0)]<-0
  		est01$ev1[is.na(est01$ev1)]<-0
		est01$ev_sqrt0<-sqrt(est01$ev0)
		est01$ev_sqrt1<-sqrt(est01$ev1)
		names(est01)[1]<-"t"
		head(est01)
		#est0<-est01[,1:6]
		#est1<-est01[,c(1,7:11)]
		names(est0)<-names(est1)<-c("t","ev","NAsurv","KMsurv","obs","atrisk")
		#est0$ev_sqrt<-sqrt(est0$ev)
		#est1$ev_sqrt<-sqrt(est1$ev)

	
		est0a<-est0[est0$ev>=1,]
		est1a<-est1[est1$ev>=1,]
		###est01a<-est01[est01$ev0==1 | est01$ev1==1,] #geordnet, nicht so gut

		#gruppe<-0
		for(gruppe in c(0,1)) {
			if(gruppe==0) est<-est0a
			if(gruppe==1) est<-est1a
			mat<-matrix(NA,nrow=dim(est)[1],ncol=dim(param)[1])
			val<-rep(NA,dim(param)[1])
			for(i in 1:dim(param)[1]) {
				temp<-getHval(i,est,param,est01,gruppe,t0,t1,ev0,ev1,haz_method)
				mat[,i]<-temp$H
				val[i]<-temp$val
			}
			#Variance with ties, assuming that the ties occured through rounding
			#x[1] = number at risk, x[2] = number of events
			vfun<-function(x)  ifelse(x[1]>0 & x[2]>0,sum(1/(x[1]-0:(x[2]-1))^2),0)

		

			if(gruppe==0) {
				mat0<-mat
				val0<-val
				#Diag0<-diag(est$ev/est$atrisk^2)
				Diag0<-diag(apply(cbind(est$atrisk,est$ev),1,vfun))

			}
			if(gruppe==1) {
				mat1<-mat
				val1<-val
				#Diag1<-diag(est$ev/est$atrisk^2)
				Diag1<-diag(apply(cbind(est$atrisk,est$ev),1,vfun))
			}
		}
		V0<-t(mat0)%*%Diag0%*%mat0
		V1<-t(mat1)%*%Diag1%*%mat1
		V<-V0+V1
		values<-val1-val0	
		obj<-list(est=values,V=V)
	} else {
		V0<-NA
		V1<-NA
		val0<-NA
		val1<-NA
		if(any(param$type%in%c("HR","score","logrank"))) stop("Perturbation approach is not available for Cox model HR or logrank test.")
		obj<-perturb_fun(t0=t0,t1=t1,ev0=ev0,ev1=ev1,param=param,Kpert=Kpert)
	}
	namen<-paste(paramin$type,paramin$par,sep="_")
	colnames(obj$V)<-rownames(obj$V)<-names(obj$est)<-namen

		
		if(alternative_test=="one.sided" | alternative_CI=="one.sided") {
			if(is.null(param$alternative)) stop("For one-sided inference, the direction of the alternative
				must be specified in the data frame handed to the argument param.")

			#change sign such that all alternatives point in same direction.
			change_sign<-ifelse(param$alternative=="less",-1,1)
			obj$est<-obj$est*change_sign
			diag_change<-diag(change_sign)
			obj$V<-diag_change%*%obj$V%*%diag_change
			if(!all(rhs==0)) rhs<-rhs*change_sign
			if(alternative_test=="one.sided") alternative_test<-"greater"
			if(alternative_CI=="one.sided") alternative_CI<-"greater"

		}
		class(obj)<-"myobj"
		G_KI<-glht(obj,alternative = alternative_CI,rhs=rhs,coef.=coef.myobj, vcov.=vcov.myobj)  #  "two.sided"
		#Not required if vcov.myobj and coef.myobj are in the namespace of the package,
		#may be required otherwise, e.h. if package is not loaded and the function is tested:
		#glht(obj,alternative = alternative_CI,rhs=rhs,coef.=coef.myobj, vcov.=vcov.myobj)
		if(alternative_test!=alternative_CI) {
			G_Test<-glht(obj,alternative = alternative_test,rhs=rhs,coef.=coef.myobj, vcov.=vcov.myobj)  #  "two.sided"
		} else {
			G_Test<-G_KI
		}
  		KI<-confint(G_KI,level=lvl)$confint
		p<-summary(G_Test)$test$pvalues

		#
		#Closed Test
		if(closed_test==TRUE) {
			clo_matrix<-clo_matrix_fun(length(obj$est))
			p_clo<-rep(NA,dim(clo_matrix)[1])
			for(i in 1:dim(clo_matrix)[1]) {
				obj_temp<-list(est=obj$est[clo_matrix[i,]],V=obj$V[clo_matrix[i,],clo_matrix[i,]])
				class(obj_temp)<-"myobj"
				G_temp<-glht(obj_temp,alternative = alternative_test,rhs=rhs[clo_matrix[i,]],coef.=coef.myobj, vcov.=vcov.myobj)
	  			p_clo[i]<-min(as.numeric(summary(G_temp)$test$pvalues))
			}
		
			#adjustierte p-Werte closed test:
			p_adj_CT<-rep(NA,dim(clo_matrix)[2])
			for(i in 1:length(p_adj_CT)) {
				p_adj_CT[i]<-max(p_clo[clo_matrix[,i]])
			}
		}

		
		#

		
		lookup<-data.frame(
			type=c("S","logS","cloglogS","cumhaz","Q","logQ","HR","avgHR","RMST","score"),
			new=c("Survival difference","Survival ratio","Cumulative-hazard ratio","Cumulative-hazard ratio",
			"Quantile difference","Quantile ratio",
			"Hazard ratio","Avg. hazard ratio","RMST difference","Score")
		)
		rownames(lookup)<-lookup$type
		#rownames(KI)<-gsub("log","ratio_",rownames(KI))
		KI<-as.data.frame(KI)
		KI$p_adjusted<-p
		KI$Parameter<-lookup[param$type,"new"]
		KI$Time_or_which_quantile<-param$par
		KI$SE<-sqrt(diag(obj$V))
		alpha<-1-lvl
		if(alternative_CI=="two.sided") {
			KI$lwr_unadj<-KI$Estimate-qnorm(1-alpha/2)*KI$SE
			KI$upr_unadj<-KI$Estimate+qnorm(1-alpha/2)*KI$SE
		} else {
			KI$lwr_unadj<-KI$Estimate-qnorm(1-alpha)*KI$SE
			KI$upr_unadj<-Inf
		}
		if(alternative_test=="two.sided") {
			KI$p_unadj<-2*(1-pnorm(abs((KI$Estimate-rhs)/KI$SE)))
		} else {
			KI$p_unadj<-1-pnorm((KI$Estimate-rhs)/KI$SE)
		}
		KI<-KI[,c(5,6,1,7:10,2:4)]
		names(KI)[names(KI)=="lwr"]<-"lwr_adjusted"
		names(KI)[names(KI)=="upr"]<-"upr_adjusted"

		if(closed_test) KI$p_adjusted_closed_test<-p_adj_CT
		KI$p_Holm<-p.adjust(KI$p_unadj,"holm")

		#reverse one sided CIs and estimates with alternative "less"
		if(alternative_CI=="greater" | alternative_test=="greater") {  #was changed from input "one.sided" to "greater" above
			#print("change")
			KI$Estimate<-KI$Estimate*change_sign
			set_change<-change_sign== -1
			KI_temp<-KI
			KI$upr_adjusted[set_change]<- -KI_temp$lwr_adjusted[set_change]
			KI$lwr_adjusted[set_change]<- -KI_temp$upr_adjusted[set_change]

			KI$upr_unadj[set_change]<- -KI_temp$lwr_unadj[set_change]
			KI$lwr_unadj[set_change]<- -KI_temp$upr_unadj[set_change]

		}

		index_exp<-grepl("lwr",colnames(KI)) | grepl("upr",colnames(KI)) | grepl("Estimate",colnames(KI))
		rows_exp<-grepl("HR",rownames(KI)) | grepl("log",rownames(KI)) | grepl("cumhaz",rownames(KI))

		KI[rows_exp,index_exp]<-exp(KI[rows_exp,index_exp])
		KI[rows_exp,"SE"]<-KI[rows_exp,"Estimate"]*KI[rows_exp,"SE"]

		rownames(KI)<-1:dim(KI)[1]

		#est01$obs0[is.na(est01$obs0)]<-0
		#est01$obs1[is.na(est01$obs1)]<-0
		est0$obs<-NULL
		est1$obs<-NULL
		ret<-list(est=obj$est,V=obj$V,tab=KI,param=param,paramin=paramin,dat0=est0,dat1=est1,minmaxt=minmaxt,est0=val0,est1=val1,V0=V0,V1=V1)
		class(ret)<-"nphparams"

		#if(plot) {
		#	gr<-rep(0:1,times=c(length(t0),length(t1)))
		#	km<-survfit(Surv(time=c(t0,t1),event=c(ev0,ev1))~gr)
		#	plot(km,col=1:2,ylab="Survival",xlab="Time",mark.time=TRUE)
		#	abline(v=param$par[grepl("S",param$type)],lty=2)
		#	abline(h=param$par[grepl("Q",param$type)],lty=2)
		#	abline(v=param$par[grepl("avgHR",param$type)],lty=3)
		#	abline(v=param$par[grepl("RMST",param$type)],lty=4)
		#	if(any(param$type=="HR")) abline(v=min(max(t0),max(t1)),lty=3)
		#}
	#}
	ret
}

plothelper<-function(est,col) {
	lines(c(0,est$t),c(1,est$NAsurv),type="s",col=col)
	diff_atrisk<-c(-diff(est$atrisk),0)
	
	zens<-est$ev==0 | est$ev < diff_atrisk
	points(est$t[zens],est$NAsurv[zens],pch=3,col=col)
}
atriskfun<-function(est,times) {
	fun_temp<-function(x) {
		if(any(est$t>=x)) max(est$atrisk[est$t>=x]) else 0
	}
	sapply(times,fun_temp)
}

#' @title Plot nphparams Objects
#'
#' @description Plots the estimated survival distributions, shows numbers at risk and indicates the requested parameters for quantifying differences between the survival curves.
#'
#' @param x an object of class \code{nphparams}.
#' @param xlim limits of the x-axis, must be a numeric vector of length two 
#' @param ylim limits of the y-axis, must be a numeric vector of length two
#' @param trt_name character, an optional name for group 1 to be shown with the number at risk table in the plot. Default is "Treatment". 
#' @param ctr_name character, an optional name for group 0 to be shown with the number at risk table in the plot. Default is "Control". 
#' @param xlab character, an optional label for the x-axis. Default is "Time". 
#' @param ylab character, an optional label for the x-axis. Default is "Survival". 
#' @param main character, an optional title of the plot. Default is "", showing  no title. 
#' @param col_ctr the color of the survival curve estimate of group 0. Default is 1 (black).
#' @param col_trt the color of the survival curve estimate of group 1. Default is 2 (red).
#' @param atrisktimes numeric vector of time-points for which the number at risk is displayed.
#' @param bold numeric, passed to linewidth and font settings. Default is 2, resulting in lines 
#'	of width 2 and boldfont. Use 1 for line-width 1 and standard font. 
#' @param showlines logical, indicating whether the time-points or the quantile-probabilites defined
#'	for the requested parametes should be shown in terms of vertical or horizontal lines. Default is \code{FALSE}.
#' @param show_rmst_diff logical, indicating whether the estiamted difference in restricted mean survival times should
#'	by visualized by a gray background area.
#' @param ... further arguments, not used
#'
#' @details
#' When setting \code{show_lines}, line type 2 (dashed) is used for survival probabilities and quantiles, line type 
#' 3 (dotted) is used for the score test, average hazard ratio and Cox model hazard ratio and
#' line type 5 (long dashed) is used for restricted mean survival time.
#'
#'
#' @author Robin Ristl, \email{robin.ristl@@meduniwien.ac.at}
#' 
#' @seealso \code{\link{nphparams}}, \code{\link{plot.nphparams}}
#'
#' @examples
#' data(pembro)
#' set1<-nphparams(time=time, event=event, group=group,data=pembro,
#'	param_type=c("score","S"),
#'	param_par=c(3.5,2),
#'	param_alternative=c("less","greater"),
#'	closed_test=TRUE,alternative_test="one.sided")
#' print(set1)
#' plot(set1,trt_name="Pembrolizumab",ctr_name="Cetuximab")
#'
#' set2<-nphparams(time=time, event=event, group=group, data=pembro,  
#'	param_type=c("S","S","S","Q","RMST"),
#'	param_par=c(0.5,1,2,0.5,3.5))
#' print(set2)
#' plot(set2,showlines=TRUE,show_rmst_diff=TRUE)
#'
#' @export
plot.nphparams<-function(x,xlim=NULL,ylim=c(0,1),trt_name="Treatment",ctr_name="Control",xlab="Time",ylab="Survival",main="",col_ctr=1,col_trt=2,atrisktimes=0:3,bold=2,showlines=FALSE,show_rmst_diff=FALSE,...) {
	
	if(is.null(xlim)) xlim<-c(0,max(x$dat0$t,x$dat1$t))
	COL<-c(col_ctr,col_trt)
	
	#par(lwd=2,font=2,font.lab=2,font.axis=2)

	old_par<-par(mar=c(5.1+5, 4.1+2, 4.1, 2.1),lwd=bold,font=bold,font.lab=bold,font.axis=bold)
	on.exit(par(old_par))
	plot.new()
	plot.window(xlim=xlim,ylim=ylim)
	axis(1,at=atrisktimes)
	axis(2)
	#title(xlab="Time",ylab="Survival",main=main)
	title(xlab=xlab,ylab=ylab,main=main)

	if(showlines) {
		abline(v=x$param$par[x$param$type%in%c("S","logS")],lty=2,col="grey")
		#abline(h=x$param$par[x$param$type%in%c("Q","logQ")],lty=2,col="grey")
		abline(h=1-x$param$par[x$param$type%in%c("Q","logQ")],lty=2,col="grey")
		abline(v=x$param$par[x$param$type%in%c("score","avgHR","HR")],lty=3,col="grey")
		abline(v=x$param$par[x$param$type%in%c("RMST")],lty=5,col="grey")
	}
	if(show_rmst_diff) {
		rmst_ind<-which(x$param$type=="RMST")
		if(length(rmst_ind)>0) {
			if(length(rmst_ind)>1) {
				rmst_ind<-rmst_ind[1]
				warning("RMST difference is included for more than one time-point. Only the first in the list is visualised in the plot.")
			}
			rmst<-x$tab$Estimate[rmst_ind]
			rmst_par<-x$param$par[rmst_ind]
			polygon(c(rmst_par-rmst,rmst_par,rmst_par,rmst_par-rmst),c(-0.2,-0.2,1.2,1.2),border=NA,col="lightgray")
		}
	}
	plothelper(x$dat0,col=COL[1])
	plothelper(x$dat1,col=COL[2])
	atrisk0<-atriskfun(x$dat0,atrisktimes)
	atrisk1<-atriskfun(x$dat1,atrisktimes)
	box()
	

	axis(1,line=4,at=atrisktimes,labels=atrisk1,lwd=0)
	axis(1,line=5,at=atrisktimes,labels=atrisk0,lwd=0)
	mtext(text="Number at risk",side=1,line=4,at=0,adj=0)
	par(old_par)
	old_par<-par(mar=c(5.1+5, 0, 4.1, 0))
	plot.window(new=TRUE,xlim=c(0,1),ylim=c(0,1))
	mtext(text=trt_name,side=1,line=5,at=-0.03,adj=0)		
	mtext(text=ctr_name,side=1,line=6,at=-0.03,adj=0)
	on.exit(par(old_par))
	par(old_par)
	old_par<-par(mar=c(0.1, 0, 0, 0))
	on.exit(par(old_par))
	plot.window(new=TRUE,xlim=c(0,1),ylim=c(0,1))
	legend("bottom",legend=c(trt_name,ctr_name),col=c(col_trt,col_ctr),lty=1,bty="n",horiz=TRUE)	
	par(lwd=1,font=1,font.lab=1,font.axis=1)
}


#' @title Print nphparams Objects
#'
#' @description Prints the results table of an nphparams object.
#' 
#' @param x an object of class \code{nphparams}.
#' @param ... further arguments, not used.
#' 
#' @details
#' Estiamtes corresponding to differences at a log scale are transformed by taking exp(), and are labelled as ratios. I.e. differnces in log urvival probabilites,
#' differences in log quantiles, cloglog survival differences (equivalent to the log cumulative hazard ratio), log average hazard ratios or Cox model log hazard ratioss are
#' transformed to survival probability ratios, quantile ratios, cumulative hazard ratios, average hazard ratios or Cox model hazard ratios, respectively. In the output,
#' the standard error at the backtransformed scale is calculated by the delta-method. Confidence interval bounds are calculated at the log-scale, though, and then
#' transformed by taking exp().
#'
#' @author Robin Ristl, \email{robin.ristl@@meduniwien.ac.at}
#' 
#' @seealso \code{\link{nphparams}}, \code{\link{plot.nphparams}}
#'
#' @examples
#' data(pembro)
#' set1<-nphparams(time=time, event=event, group=group,data=pembro,
#'	param_type=c("score","S"),
#'	param_par=c(3.5,2),
#'	param_alternative=c("less","greater"),
#'	closed_test=TRUE,alternative_test="one.sided")
#' print(set1)
#' plot(set1,trt_name="Pembrolizumab",ctr_name="Cetuximab")
#'
#' set2<-nphparams(time=time, event=event, group=group, data=pembro,  
#'	param_type=c("S","S","S","Q","RMST"),
#'	param_par=c(0.5,1,2,0.5,3.5))
#' print(set2)
#' plot(set2,showlines=TRUE,show_rmst_diff=TRUE)
#'
#' @export
print.nphparams<-function(x,...) {
	cat("Simultaneous inference for differences between two survival functions\n\n")
	print(format(x$tab,digits=4))
}

#NOTE: These need to be added to the namespace file!
vcov.myobj<-function(object,...) {
	#man muss die Namen richtig stellen, damit es immer mit glht funktioniert:
	V<-as.matrix(object$V)
	colnames(V)<-names(object$est)
	rownames(V)<-names(object$est)
	return(V)
}
coef.myobj<-function(object,...) object$est

vcov.nphparams<-function(object,...) {
	V<-as.matrix(object$V)
	colnames(V)<-names(object$est)
	rownames(V)<-names(object$est)
	return(V)
}


#Perturbation

#Variance with ties, assuming that the ties occured through rounding
#x[1] = number at risk, x[2] = number of events
vfun<-function(x)  ifelse(x[1]>0 & x[2]>0,sum(1/(x[1]-0:(x[2]-1))^2),0)

rmst_fun<-function(S,time,para) {
	pos<-sum(time<=para)
	set<-time<=para
	times<-time[set]
	S<-S[set]
	times<-c(0,times)
	S<-c(1,S)
	times<-c(times,para)
	dt<-diff(times)
	val<-sum(S*dt)
	val
}

perturb_fun<-function(t0,t1,ev0,ev1,param,Kpert=500) {
	km0<-survfit(Surv(time=t0,event=ev0)~1)
	km1<-survfit(Surv(time=t1,event=ev1)~1)
	#log_q<-log(q)

	time0<-km0$time
	time1<-km1$time
	dN_Y0<-km0$n.event/km0$n.risk
	dN_Y1<-km1$n.event/km1$n.risk
	NAest0_basis<-cumsum(dN_Y0)
	NAest1_basis<-cumsum(dN_Y1)
	#data.frame(km0$time,km0$n.event,km0$n.risk,km0$surv)
	sd_incr0<-sqrt(apply(cbind(km0$n.risk,km0$n.event),1,vfun))
	sd_incr1<-sqrt(apply(cbind(km1$n.risk,km1$n.event),1,vfun))
	
	Z0<-matrix(rnorm(length(km0$time)*Kpert,dN_Y0,sd_incr0),ncol=Kpert)
	Z1<-matrix(rnorm(length(km1$time)*Kpert,dN_Y1,sd_incr1),ncol=Kpert)

	#
	#Originaldaten als erste Spalte hinzufügen
	Z0<-cbind(dN_Y0,Z0)
	Z1<-cbind(dN_Y1,Z1)

	NA0<-apply(Z0,2,cumsum)
	NA1<-apply(Z1,2,cumsum)
	Pert0<-exp(-NA0)
	Pert1<-exp(-NA1)
	
	#matplot(time0,Pert0,type="l")
	val<-matrix(NA,Kpert+1,dim(param)[1])
	#i<-1
	#j<-2
	for(i in 1:dim(param)[1]) {
		for(j in 1:(Kpert+1)) {
			
			if(param$type[i]%in%c("S","logS","cloglogS")) {
				pos0<-sum(time0<=param$par[i])
				pos1<-sum(time1<=param$par[i])
				if(param$type[i]=="S") val[j,i]<-Pert1[pos1,j]-Pert0[pos0,j]
				if(param$type[i]=="logS") val[j,i]<-NA0[pos0,j]-NA1[pos1,j] #NA=-log S
				if(param$type[i]=="cloglogS") val[j,i]<-log(NA1[pos1,j])-log(NA0[pos0,j])
			}
			if(param$type[i]%in%c("Q","logQ")) {
				#pos<-sum(est$NAsurv>(1-param$par[i]))+1
				pos0<-sum(Pert0[,j]>param$par[i])+1
				pos1<-sum(Pert1[,j]>param$par[i])+1
				if(param$type[i]=="Q") val[j,i]<- time1[pos1]-time0[pos0]
				if(param$type[i]=="logQ") val[j,i]<- log(time1[pos1])-log(time0[pos0])
			}
			if(param$type[i]=="avgHR") {
				est0<-data.frame(t0=time0,NAsurv0=Pert0[,j],obs0=1)
				est1<-data.frame(t1=time1,NAsurv1=Pert1[,j],obs1=1)

				est01<-merge(est0,est1,by.x="t0",by.y="t1",all=TRUE)
		  		est01$NAsurv0<-fillo(est01$NAsurv0)
  				est01$NAsurv1<-fillo(est01$NAsurv1)
  
		 		#predictable
				weights01<-(c(1,est01$NAsurv0)*c(1,est01$NAsurv1))[-(1+dim(est01)[1])]
				weights0<-weights01[!is.na(est01$obs0)]
				weights1<-weights01[!is.na(est01$obs1)]
				L<-param$par[i]
			  	w_dN_Y0<-(weights0*Z0[,j])[time0<=L]
			  	w_dN_Y1<-(weights1*Z1[,j])[time1<=L]
				avgHR_est0<-sum(w_dN_Y0)
				avgHR_est1<-sum(w_dN_Y1) #kann bei kleiner Fallzhal in seltenen Faellen negative werden...
				val[j,i]<-log(avgHR_est1)-log(avgHR_est0)
			}
			if(param$type[i]=="RMST") val[j,i]<-rmst_fun(Pert1[,j],time1,param$par[i])-rmst_fun(Pert0[,j],time0,param$par[i])
		}
	}
	est<-val[1,]
	V<-cov(val[-1,],use="pairwise.complete.obs")
	list(est=est,V=V)
}

###Example Data
#' Reconstructed Data Set Based On Survival Curves In Burtess et al. 2019
#'
#' The data set was approximately reconstructed from the survival curves shown in Figure 2D of Burtness et al. 2019 (see references). It contains survival times and event #' indicator under treatment with pembrolizumab (group 1) versus cetuximab with chemotherapy (group 0).
#'
#' @docType data
#'
#' @usage data(pembro)
#'
#' @format A data frame.
#'
#' @keywords datasets
#'
#' @references
#' Burtness, B., Harrington, K. J., Greil, R., Soulières, D., Tahara, M., de Castro Jr, G., ... & Abel, E. (2019). Pembrolizumab alone or with chemotherapy
#' versus cetuximab with chemotherapy for recurrent or metastatic squamous cell carcinoma of the head and neck (KEYNOTE-048): a randomised, open-label,
#' phase 3 study. The Lancet, 394(10212), 1915-1928.
#'
#' @examples
#' data(pembro)
#' head(pembro)
"pembro"

