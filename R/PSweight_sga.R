#' Estimate subgroup causal effects by propensity score weighting
#'
#' The function \code{PSweight_sga} is used to estimate the subgroup average potential outcomes corresponding to
#' each treatment group among the target population. The function currently implements
#' three types of weights: the inverse probability weights (target population is the combined population),
#' ATT weights (target population is the population receiving one treatment) and overlap weights (target
#' population is the overlap population at clinical equipoise).
#'
#' @param ps.formula an object of class \code{\link{formula}} (or one that can be coerced to that class):
#' a symbolic description of the propensity score model to be fitted. Additional details of model specification
#' are given under "Details". This argument is optional if \code{ps.estimate} is not \code{NULL}.
#' @param ps.estimate an optional matrix or data frame containing estimated (generalized) propensity scores for
#' each observation. Typically, this is an N by J matrix, where N is the number of observations and J is the
#' total number of treatment levels. Preferably, the column name of this matrix should match the name of treatment level,
#' if column name is missing or there is a mismatch, the column names would be assigned according to alphabatic order
#' of the treatment levels. A vector of propensity score estimates is also allowed in \code{ps.estimate}, in which
#' case a binary treatment is implied and the input is regarded as the propensity to receive the last category of
#' treatment by alphabatic order, unless otherwise stated by \code{trtgrp}.
#' @param subgroup a vector to specify name of subgroup variables by column index or column names
#' @param xname an optional character vector specifying the name of the covariates (confounders) in \code{data}. Only continuous and factors are accepted.
#' @param trtgrp an optional character defining the "treated" population for estimating the average treatment
#' effect among the treated (ATT). Only necessary if \code{weight = "ATT"}. This option can also be used to specify
#' the treatment (in a two-treatment setting) when a vector argument is supplied for \code{ps.estimate}.
#' Default value is the last group in the alphebatic order.
#' @param zname an optional character specifying the name of the treatment variable in \code{data}.
#' @param yname an optional character specifying the name of the outcome variable in \code{data}.
#' @param data an optional data frame containing the variables in the propensity score model.
#' @param R an optional integer indicating number of bootstrap replicates. Default is \code{R = 50}.
#' @param weight a character or vector of characters including the types of weights to be used.
#' \code{"ATE"} specifies the inverse probability weights for estimating the average treatment
#' effect among the combined population. \code{"ATT"} specifies the weights for estimating the
#' average treatment effect among the treated. \code{"ATO"} specifies the (generalized) overlap weights
#' for estimating the average treatment effect among the overlap population, or population at
#' clinical equipoise. Default is \code{"ATO"}.
#' @param method a character to specify the method for propensity model. When \code{ps.formula} is given, \code{"glm"} is the default; When \code{xname} is given, \code{"LASSO"} is the default.

#' @details A typical form for \code{ps.formula} is \code{treatment ~ terms} where \code{treatment} is the treatment
#' variable (identical to the variable name used to specify \code{zname}) and \code{terms} is a series of terms
#' which specifies a linear predictor for \code{treatment}. The \code{ps.formula} specifies generalized
#' linear models when \code{ps.estimate} is \code{NULL}. See \code{glm} for more details on generalized linear models.
#'
#' In Yang et al.(2021), the \code{term} is suggested to include all main effects and pairwise interactions between subgroup variables
#' and confounders. Then LASSO will be performed to select important interactions, and re-fit a logistic regression with the main
#' effects and LASSO selected interactions (pLASSO). If \code{xname}, \code{zname} and \code{subgroup} are provided, the function will automatically perform pLASSO, and the \code{ps.formula} is not required.
#'
#' When comparing two treatments, \code{ps.estimate} can either be a vector or a two-column matrix of estimated
#' propensity scores. If a vector is supplied, it is assumed to be the propensity scores to receive the treatment, and
#' the treatment group corresponds to the last group in the alphebatic order, unless otherwise specified by \code{trtgrp}.
#' In general, \code{ps.estimate} should have column names that indicate the level of the treatment variable,
#' which should match the levels given in \code{Z}.
#' If column names are empty or there is a mismatch, the column names will be created following
#' the alphebatic order of values in \code{Z}, and the rightmost coulmn of \code{ps.estimate} is assumed
#' to be the treatment group, when estimating ATT. \code{trtgrp} can also be used to specify the treatment
#' group for estimating ATT.
#'
#' The argument \code{zname} and/or \code{xname} is required when \code{ps.estimate} is not \code{NULL}.
#'
#' Current version of \code{PSweight_sga} allows for three types of propensity score weights used to estimate ATE, ATT and
#' ATO. These weights are members of larger class of balancing weights defined in Li, Morgan, and Zaslavsky (2018). The overlap weights can also be considered as
#' a data-driven continuous trimming strategy without specifying trimming rules, see Li, Thomas and Li (2019).
#' Additional details on balancing weights and generalized overlap weights for multiple treatment groups are provided in
#' Li and Li (2019).
#'
#' The variance will be calculated by nonparametric bootstrap, with \code{R} bootstrap
#' replications. The default of \code{R} is 50.
#'
#' @return PSweight_sga returns a \code{PSweight_sga} object containing a list of the following values:
#' estimated propensity scores, average subgroup potential outcomes corresponding to each treatment,
#' estimates in each bootstrap replicate, the label for each treatment group, and  the number of subgroup
#' levels defined by each subgrouping variable .
#' A summary of PSweight_sga can be obtained with \code{\link{summary.PSweight_sga}}.
#'
#' \describe{
#'
#'
#' \item{\code{ propensity}}{a data frame of estimated propensity scores.}
#'
#' \item{\code{ muhat}}{ average subgroup potential outcomes by treatment groups, with reference to specific target populations.}
#'
#' \item{\code{ muboot}}{list of point estimates in each bootstrap replicate.}
#'
#' \item{\code{ group}}{ a table of treatment group labels corresponding to the output point estimates \code{muhat}.}
#'
#' \item{\code{ trtgrp}}{a character indicating the treatment group.}
#'
#' \item{\code{ sub_n}}{a vector indicating the number of subgroup levels defined by each subgrouping variable .}
#'
#' \item{\code{ method}}{a character indicating the propensity score method used.}
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
#' Mao, H., Li, L., Greene, T. (2019). Propensity score weighting analysis and treatment effect discovery.
#' Statistical Methods in Medical Research, 28(8), 2439-2454.
#'
#' Li, F., Thomas, L. E., Li, F. (2019).
#' Addressing extreme propensity scores via the overlap weights. American Journal of Epidemiology, 188(1), 250-257.
#'
#'
#' @export
#'
#'
#' @import nnet
#' @import MASS
#' @import numDeriv
#' @importFrom  stats binomial coef cov formula glm lm model.matrix plogis poisson predict qnorm quantile sd as.formula var
#' @importFrom  utils capture.output combn
#' @importFrom  graphics hist legend
#'
PSweight_sga <- function(ps.formula=NULL,ps.estimate=NULL, subgroup=NULL, xname=NULL, trtgrp=NULL,zname=NULL,yname,data,R=50,weight="overlap",method="glm"){


  #####NO subgroup specified, abort the program########################################################################
    if(is.null(subgroup)){
      stop("No subgroup specified")
    }
  #####If xname is provided, automatically perform pLASSO ##############################################################
  if(is.null(xname)==F){
    if(  is.null(ps.formula)==F  ){
      warning("Both xname and ps.formula are deteced. The pLASSO will be performed and the ps.formula will be ignored.")
    }
    method <- "LASSO"
    exclude_var<-subgroup
    all_var<-xname
    include_var<-setdiff(all_var,exclude_var)
    ps.formula <- paste(zname,"~",paste0(c(paste0(include_var, collapse=" + "),paste0(subgroup, collapse=" + ")), collapse=" + "))

    for(i in 1:length(subgroup)){
      for (j in 1:length(include_var)){
        ps.formula=paste(ps.formula, paste(subgroup[i],":",include_var[j]), sep=" + ")
      }
    }
  }

    #####Extract zname, treatment label and relabel the treatment to 0/1##################################################
    if(is.null(ps.estimate)){
      ps.formula<-as.formula(ps.formula)
      zname<-all.vars(ps.formula)[1]

      #option to include all covariates
      if(all.vars(ps.formula)[2]=='.'){
        #exclude some variables
        exclude_var<-c(subgroup,zname,yname)
        all_var<-colnames(data)
        include_var<-setdiff(all_var,exclude_var)

        formulatmp<- paste(zname,'~')
        formulatmp<- paste(formulatmp,paste(include_var,collapse = "+"),'+')
        formulatmp<- paste(formulatmp,paste(subgroup,collapse = "+"))

        for(i in subgroup){
          #interaction terms
          inter_tmp<-paste0(include_var, ':', i)
          inter_tmp<-paste("+",paste(inter_tmp,collapse = " + "))
          formulatmp<-paste(formulatmp,inter_tmp)
        }
        ps.formula<-as.formula(formulatmp)
      }
    }

    #set ordered group to match group later
    facz<-as.factor(unlist(data[zname]))


    #set the treatment label
    oldlevel<-levels(facz)
    newlevel<-oldlevel

    z<-as.numeric(facz)-1
    n<-length(z) #total obs

    #set group for ATT
    if(weight=="treated"){
      if(is.null(trtgrp)){
        trtgrp<-oldlevel[2]
      }else{
        newlevel<-rev(unique(c(trtgrp,oldlevel)))
        facz<-factor(unlist(data[zname]),levels = newlevel)
        z<-as.numeric(facz)-1
      }
    }

    matchlevel<-match(oldlevel,newlevel) #match level for ATT

    #reassign value as numeric
    y<-unlist(data[yname])
    data[zname]<-z
    data_p<-data[,colnames(data)!=yname] #data without outcome



    ######## Prepare subgroup matrix for estimation#############################################################################
    submatrix<-as.matrix(data[,subgroup,drop=FALSE])
    submatrix<-cbind(1,submatrix) #include the overall effect
    colnames(submatrix)[1]<-"Overall"

    submatrix_level<-NULL
    level_name<-c()
    sub_n<-c()

    # matrix for subgroup
    for(r in 1:ncol(submatrix)){
      leveltmp <- sort(unique(submatrix[,r]))
      sub_n<-c(sub_n,length(leveltmp))
      for(g in leveltmp){
        submatrix_level<-cbind(submatrix_level,c(submatrix[,r]==g)*1)
        level_name<-c(level_name,paste0(colnames(submatrix)[r],"=",g))
      }
    }
    colnames(submatrix_level)<-level_name
    colnames(submatrix_level)[1]<-"Overall"


    ####### Estimate ps ########################################################################################################

    if(is.null(ps.estimate)){
      e.h <- do.call(PSmethod_sga,c(list(ps.formula = ps.formula, method=method, data=data_p, weight=weight)))
      if(method=="LASSO"){e.h <- e.h$e.h}
      e.h<-c(e.h[,2])
    }else{
      #the name for the propensity score
      if(length(ps.estimate)==n){
        e.h<-c(ps.estimate)
      }else if(!setequal(colnames(ps.estimate),newlevel)){
        ps.estimate<-as.matrix(ps.estimate)
        e.h<-c(ps.estimate[,2])
        warning("wrong column name set for ps.estimate, treatment set as: ",newlevel[1], " , ", newlevel[2])
      }else{
        ps.estimate<-ps.estimate[,match(newlevel,colnames(ps.estimate))]
        e.h<-c(ps.estimate[,2])
      }
    }

    ###### Genearate point estimate and evaluate uncertainty through bootstrap ##############################################
    if(weight=="entropy"){
      e.h<-pmax(e.h,1e-6)
      e.h<-pmin(e.h,1-1e-6)
    }

    #tilting function
    ftilt<-tiltbin(weight = weight)

    #calculate point
    muhat<-pt_est_sga(psest=e.h,y=y,z=z,submatrix_level=submatrix_level,ftilt=ftilt,newlevel=newlevel)


    #bootstrap
    mu1.boot<-NULL
    mu0.boot<-NULL
    for(i in 1:R){
      if(i %% 50==0){
        message("bootstrap ", i, " samples")
      }
      # estimate ps
      samp.b<-sample(n,n,replace = TRUE)
      data.b<-data_p[samp.b,]
      y.b<-y[samp.b]
      z.b<-z[samp.b]
      submatrix_level.b<-submatrix_level[samp.b,]
      if(is.null(ps.estimate)&method=='glm'){
        fit.b <- glm(formula = ps.formula, data=data.b, family = binomial(link = "logit"))
        e.b <- as.numeric(fit.b$fitted.values)
      }else{
        e.b<-e.h[samp.b]
      }


      #calculate point
      muhat.b<-pt_est_sga(psest=e.b,y=y.b,z=z.b,submatrix_level=submatrix_level.b,ftilt=ftilt,newlevel=newlevel)
      mu1.boot<-rbind(mu1.boot,c(muhat.b[,2]))
      mu0.boot<-rbind(mu0.boot,c(muhat.b[,1]))
    }
    rownames(mu1.boot)<-NULL
    rownames(mu0.boot)<-NULL

    muboot<-list(mu0.boot,mu1.boot)
    names(muboot)<-newlevel
    e.h<-cbind(1-e.h,e.h)
    colnames(e.h)<-newlevel
    names(sub_n)<-c("Overall",subgroup)

    #match back to original levels
    e.h<-e.h[,matchlevel]
    muhat<-muhat[,matchlevel]
    muboot<-muboot[matchlevel]

    out<-list(propensity=e.h, muhat=muhat, muboot=muboot, group=c(oldlevel), trtgrp=trtgrp, sub_n=sub_n, method=method)
    class(out)<-'PSweight_sga'
    out
  }






