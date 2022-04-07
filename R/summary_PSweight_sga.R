#' Summarize a PSweight_sga object
#'
#' \code{summary.PSweight_sga} is used to summarize the results from \code{\link{PSweight_sga}}.
#' The output contains the average subgroup causal effects defined by specific contrasts, as well as their
#' standard error estimates.
#'
#' @param object a PSweight_sga object obtained from the \code{\link{PSweight_sga}} function.
#' @param contrast a vector or matrix specifying the causal contrast of interest. The average causal effects will be
#' defined by such contrats. For multiple treatments, the contrast parameters are explained in Li and Li (2019)
#' for estimating general causal effects. Default is all pairwise contrasts between any two treatment groups.
#' @param type a character specifying the target estimand. The most commonly seen additive estimand is specified
#' by \code{type = "DIF"}, abbreviated for weighted difference-in-means. This is the usual pairwise average treatment
#' effects as those defined in Li, Morgan, and Zaslavsky (2018) and Li and Li (2019). For binary (or count outcomes), we also
#' allow two ratio estimands: causal relative risk (\code{type = "RR"}) and causal odds ratio (\code{type = "OR"}).
#' Estimates for these two ratio estimands will be reported on the log scale (log relative risk and log
#' odds ratio) to improve the approximate for asymptotic normality. With binary outcomes, \code{"DIF"} is the same
#' as the average causal risk difference. Default is "DIF" if left empty.
#' @param het an indicator specifying whether to summarize the test of heterogeneity across subgroup levels. The default is FALSE.
#' @param ... further arguments passed to or from other methods.
#'
#' @details For the \code{contrast} argument, one specifies the contrast of interest and thus defines the target estimand
#' for comparing treatments. For example, if there are two treatment levels: A and B, the contrast A-B
#' (i.e., E[Y(A)] - E[Y(B)]) can be specified by \code{c(1,-1)}.
#'
#' For estimating the causal relative risk (\code{type = "RR"}), the contrast is specified at the log scale. For example,
#' the contrast A-B (specified by \code{c(1,-1)}) implies the estimation of log\{E[Y(A)]\} - log\{E[Y(B)]\}. For estimating the causal odds
#' ratio, the contrast is specified at the log odds scale. For example, the contrast A-B (specified by \code{c(1,-1)})
#' implies the estimation of log\{E[Y(A)]/E[1-Y(A)]\} - log\{E[Y(B)]/E[1-Y(B)]\}.
#'
#' The variance of the contrasts will be estimated by nonparametric bootstrap.
#'
#' The argument \code{type} takes one of three options: \code{"DIF"}, \code{"RR"}, or \code{"RR"}, with \code{"DIF"} as
#' the default option. Typically, \code{"RR"} is relavent for binary or count outcomes, and \code{"OR"} is relavent
#' only for binary outcomes. \code{"DIF"} applies to all types of outcomes.
#'
#' @return A list of following values:
#'
#' \describe{
#' \item{\code{ trtgrp}}{ a character indicating the treatment group, or target population under ATT weights.}
#'
#' \item{\code{ estimates}}{ a matrix of subgroup point estimates, standard errors and 95% confidence intervals
#' for contrasts of interest.}
#'
#'
#' \item{\code{ contrast}}{ a table listing the specified contrasts of interest.}
#'
#' \item{\code{ group}}{ a table of treatment group labels corresponding to the output point estimates, provided in results
#' obtained from \code{\link{PSweight_sga}}.}
#'
#' }
#'
#' @references
#' Yang, S., Lorenzi, E., Papadogeorgou, G., Wojdyla, D. M., Li, F., & Thomas, L. E. (2021).
#' Propensity score weighting for causal subgroup analysis. Statistics in medicine, 40(19), 4294-4309.
#'
#' Li, F., Morgan, K. L., Zaslavsky, A. M. (2018).
#' Balancing covariates via propensity score weighting.
#' Journal of the American Statistical Association, 113(521), 390-400.
#'
#'
#' @export
#'
#'
#' @importFrom  stats binomial coef cov formula glm lm model.matrix plogis poisson predict qnorm quantile sd
#' @importFrom  utils capture.output combn
#' @importFrom  graphics hist legend
#'
summary.PSweight_sga<-function(object,contrast=NULL,type='DIF',het=FALSE,...){

  muhat<-object$muhat
  muboot<-object$muboot
  muboot0<-muboot[[1]]
  muboot1<-muboot[[2]]
  trtgrp<-object$trtgrp
  group<-object$group
  sub_n<-object$sub_n

  sub_group<-rownames(muhat)



  #groups
  ngrp<-length(group)


  #transform contrast into matrix
  if(is.vector(contrast)){
    contrast<-t(contrast)
  }

  #error message for wrong contrast
  if(!is.null(contrast)){
    if(dim(contrast)[2]!=ngrp){
      cat('Contract length not equal to treatment groups, please check','\n')
      cat('\n')
      contrast<-NULL
    }
    cat('\n')
  }




  #if all contrast
  if(is.null(contrast)){
    ncst<-ngrp*(ngrp-1)/2
    contrasttmp<-matrix(0,ngrp-1,ngrp)
    contrasttmp[,1]<--1
    contrasttmp[1:(ngrp-1),2:ngrp]<-diag(1,ngrp-1,ngrp-1)
    contrast<-contrasttmp
    if(ngrp>2){
      for(i in 1:(ngrp-2)){
        ntmp<-nrow(contrasttmp)-1
        contrasttmp<-cbind(0,contrasttmp)[1:ntmp,1:ngrp]
        contrast<-rbind(contrast,contrasttmp)
      }
    }
  }

  rownames(contrast)<-paste('Contrast',1:nrow(contrast))
  colnames(contrast)<-group

  #define the bootstrap p value function
  pvalboot <- function(x,est) {
    esttmp<-abs(est)
    x<-abs(x-mean(x))
    return(mean(x>esttmp))
  }


  #calculate interval and test statistics for subgroup
  int_pval<-function(muboottmp,muhattmp,type=type,contrast=contrast){
    if(type=='DIF'){
      samp<-muboottmp%*%t(contrast)

      est.h<-c(contrast%*%muhattmp)
      se.h<-sqrt(diag(cov(samp)))
      lcl<-apply(samp,2,function(x) quantile(x,0.025,na.rm=TRUE))
      ucl<-apply(samp,2,function(x) quantile(x,0.975,na.rm=TRUE))
      p.value<-c()
      dimcon<-dim(samp)[2]
      for(j in 1:dimcon){
        esttmp<-c(est.h)[j]
        p.value<-c(p.value,pvalboot(samp[,j],esttmp))
      }

    }else if(type=='RR'){
      samp<-log(muboottmp)%*%t(contrast)
      tranmuhat<-log(muhattmp)
      tranest<-c(contrast%*%tranmuhat)
      est.h<-tranest
      se.h<-sqrt(diag(cov(samp)))
      lcl<-apply(samp,2,function(x) quantile(x,0.025))
      ucl<-apply(samp,2,function(x) quantile(x,0.975))
      p.value<-c()
      dimcon<-dim(samp)[2]
      for(j in 1:dimcon){
        esttmp<-c(est.h)[j]
        p.value<-c(p.value,pvalboot(samp[,j],esttmp))
      }

    }else if(type=='OR'){
      samp<-log(muboottmp/(1-muboottmp))%*%t(contrast)
      tranmuhat<-log(muhattmp/(1-muhattmp))
      tranest<-c(contrast%*%tranmuhat)
      est.h<-tranest
      se.h<-sqrt(diag(cov(samp)))
      lcl<-apply(samp,2,function(x) quantile(x,0.025))
      ucl<-apply(samp,2,function(x) quantile(x,0.975))
      p.value<-c()
      dimcon<-dim(samp)[2]
      for(j in 1:dimcon){
        esttmp<-c(est.h)[j]
        p.value<-c(p.value,pvalboot(samp[,j],esttmp))
      }

    }else{
      cat('type not found')
      cat('\n')
    }
    estimatestmp<-cbind(est.h,se.h,lcl,ucl,p.value)
    colnames(estimatestmp)<-c("Estimate","Std.Error","Lower.CL","Upper.CL","p.value")
    rownames(estimatestmp)<-rownames(contrast)
    return(estimatestmp)

  }

  estimates<-list()
  for(g in sub_group){
    muhattmp<-muhat[g,]
    muboottmp<-cbind(muboot0[,g],muboot1[,g])
    estimates[[g]]<-int_pval(muboottmp,muhattmp,type = type,contrast = contrast)
  }




  #test for hetrogenous effect
  het_eval<-NULL
  het_grp<-function(muboottmp,muhattmp){
    samp<-muboottmp%*%c(1,-1)
    est.h<-c(c(1,-1)%*%muhattmp)
    return(list(samp=samp,est.h=est.h))
  }

  if(het){
    het_eval<-list()

    #produce the point estimate and bootstrap sample for heterogeneous test
    het_est<-c()
    het_samp<-NULL
    for(g in sub_group[-1]){
      muhattmp<-muhat[g,]
      muboottmp<-cbind(muboot0[,g],muboot1[,g])
      het_tmp<-het_grp(muboottmp,muhattmp)
      het_est<-c(het_est,c(het_tmp$est.h))
      het_samp<-cbind(het_samp,het_tmp$samp)
    }


    count=1
    for(g1 in 1:(length(sub_n)-1)){ #remove the overall
      ngrp<-sub_n[g1+1]
      nametmp<-names(sub_n)[g1+1]
      idxtmp<-count:(count+ngrp-1)
      hetmutmp<-het_est[idxtmp]
      hetsamptmp<-het_samp[,idxtmp]
      subnametmp<-sub_group[1+idxtmp]

      count<-count+ngrp

      #specify the contrast for hetrogenous test
      ncst<-ngrp*(ngrp-1)/2
      contrasttmp<-matrix(0,ngrp-1,ngrp)
      contrasttmp[,1]<--1
      contrasttmp[1:(ngrp-1),2:ngrp]<-diag(1,ngrp-1,ngrp-1)
      contrast_het<-contrasttmp
      if(ngrp>2){
        for(i in 1:(ngrp-2)){
          ntmp<-nrow(contrasttmp)-1
          contrasttmp<-cbind(0,contrasttmp)[1:ntmp,1:ngrp]
          contrast_het<-rbind(contrast_het,contrasttmp)
        }
      }

      het_eval[[g1]]<-int_pval(hetsamptmp,hetmutmp,type = 'DIF',contrast = contrast_het)

      ones<-apply(contrast_het,1,function(x) which(x==1))
      nones<-apply(contrast_het,1,function(x) which(x==-1))
      rownames(het_eval[[g1]])<- paste(subnametmp[ones],'-',subnametmp[nones],": ")
    }
    names(het_eval)<-names(sub_n)[-1]
  }


  out<-list(estimates=estimates,contrast=contrast,group=group,trtgrp=trtgrp,het_eval=het_eval)
  class(out)<-'PSweightsum_sga'
  out
}
