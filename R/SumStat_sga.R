#' Calculate summary statistics for propensity score weighting in the subgroups
#'
#' \code{SumStat_sga} is used to generate distributional plots of the estimated propensity scores and balance
#' diagnostics after propensity score weighting.
#'
#' @param subgroup a vector to specify name of subgroup variables by column index or column names
#' @param xname an optional character vector specifying the name of the covariates (confounders) in \code{data}. Only continuous and factors are accepted.
#' @param ps.formula an object of class \code{\link{formula}} (or one that can be coerced to that class): a symbolic description of the propensity score model to be fitted. Additional details of model specification are given under "Details". This argument is optional if \code{ps.estimate} is not \code{NULL}.
#' @param ps.estimate an optional matrix or data frame containing estimated (generalized) propensity scores for each observation. Typically, this is an N by J matrix, where N is the number of observations and J is the total number of treatment levels. Preferably, the column names of this matrix should match the names of treatment level, if column names are missing or there is a mismatch, the column names would be assigned according to the alphabatic order of treatment levels. A vector of propensity score estimates is also allowed in \code{ps.estimate}, in which case a binary treatment is implied and the input is regarded as the propensity to receive the last category of treatment by alphabatic order, unless otherwise stated by \code{trtgrp}.
#' @param trtgrp an optional character defining the "treated" population for estimating the average treatment effect among the treated (ATT). Only necessary if \code{weight = "ATT"}. This option can also be used to specify the treatment (in a two-treatment setting) when a vector argument is supplied for \code{ps.estimate}. Default value is the last group in the alphebatic order.
#' @param zname an optional character specifying the name of the treatment variable in \code{data}.
#' @param yname an optional vector of characters including the names of outcome in \code{data}.
#' @param data an optional data frame containing the variables in the propensity score model. If not found in data, the variables are taken from \code{environment(formula)}.
#' @param weight a character or vector of characters including the types of weights to be used. \code{"ATE"} specifies the inverse probability weights for estimating the average treatment effect among the combined population. \code{"ATT"} specifies the weights for estimating the average treatment effect among the treated. \code{"ATO"} specifies the (generalized) overlap weights for estimating the average treatment effect among the overlap population, or population at clinical equipoise. Default is \code{"ATO"}.
#' @param method a character to specify the method for propensity model. When \code{ps.formula} is given, \code{"glm"} is the default; When \code{cov} is given, \code{"LASSO"} is the default.

#' @details A typical form for \code{ps.formula} is \code{treatment ~ terms} where \code{treatment} is the treatment
#' variable (identical to the variable name used to specify \code{zname}) and \code{terms} is a series of terms
#' which specifies a linear predictor for \code{treatment}. \code{ps.formula} specifies logistic or more flexible
#' models for estimating the propensity scores, when \code{ps.estimate} is \code{NULL}.
#'
#' When comparing two treatments, \code{ps.estimate} can either be a vector or a two-column matrix of estimated
#' propensity scores. If a vector is supplied, it is assumed to be the propensity scores to receive the treatment, and
#' the treatment group corresponds to the last group in the alphabetic order, unless otherwise specified by \code{trtgrp}.
#' In general, \code{ps.estimate} should have column names that indicate the level of the treatment variable,
#' which should match the levels given in \code{Z}.
#' If column names are empty or there is a mismatch, the column names will be created following
#' the alphabetic order of values in \code{Z}, and the rightmost column of \code{ps.estimate} is assumed
#' to be the treatment group, when estimating ATT. \code{trtgrp} can also be used to specify the treatment
#' group for estimating ATT.
#'
#' The argument \code{zname} and/or \code{yname} is required when \code{ps.estimate}
#' is not \code{NULL}.
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
#' @return SumStat_sga returns a \code{SumStat_sga} object including a list of the following value:
#' treatment group, propensity scores, propensity score weights, effective sample sizes,
#' and balance statistics.
#'
#' \describe{
#' \item{\code{ trtgrp}}{a character indicating the treatment group.}
#'
#' \item{\code{ propensity}}{a data frame of estimated propensity scores.}
#'
#' \item{\code{ ps.weight}}{a data frame of propensity score weights.}
#'
#' \item{\code{ ASD}}{ a table including absolute standardized mean differences in the overall sample and subgroups after weighting.}
#'
#' \item{\code{ ASD_bs}}{ a table including absolute standardized mean differences in the overall sample and subgroups before weighting.}
#'
#' \item{\code{ vif}}{ a vector indicating the approximated variance inflation in the overall sample and subgroups after weighting, see Yang et al. (2021)}
#'
#' \item{\code{ nsubg}}{ a vector indicating the subgroup sample sizes.}
#'
#' \item{\code{ ess}}{a table of effective sample sizes. This serves as a conservative measure to
#' characterize the variance inflation or precision loss due to weighting, see Yang et al. (2021).}
#'
#' \item{\code{ subgoup}}{ a vector indicating name of the specified subgroups.}
#'
#' \item{\code{ method}}{a character indicating the propensity score method used.}
#'
#' \item{\code{ nonzero_coef}}{a vector indicating the terms selected by LASSO. Only available when \code{method} is \code{LASSO}.}
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
#' Li, F., Thomas, L. E., Li, F. (2019).
#' Addressing extreme propensity scores via the overlap weights. American Journal of Epidemiology, 188(1), 250-257.
#'
#' @export
#'
#' @import nnet
#' @import PSweight
#' @import glmnet

#' @importFrom  stats binomial coef cov formula glm lm model.matrix plogis poisson predict qnorm quantile sd
#' @importFrom  utils capture.output combn
#' @importFrom  graphics hist legend

SumStat_sga <- function(subgroup=NULL,xname=NULL, ps.formula=NULL,ps.estimate=NULL,trtgrp=NULL, zname=NULL, yname=NULL, data=NULL,  weight= "overlap",method="glm"){



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
    allvarname<-all.vars(ps.formula)
    covname<-allvarname[-1]
  }

  #set ordered group
  facz<-as.factor(unlist(data[zname]))
  ncate<-nlevels(facz) #number of categories

  #set the treatment label
  dic<-levels(facz)

  znum<-as.numeric(facz) #numeric z
  z<-znum-1 #z to fit the model

  n<-length(z) #total obs

  #extract z
  data[,zname]<-z

  #set group for ATT
  if(weight=="treated"){
    if(is.null(trtgrp)){
      trtgrp<-dic[2]
    }
  }

  #pick out the treatment group index if treated is specified,
  #this could be useful when we have att or a single column of ps.estimate
  trt<-2
  if (!is.null(trtgrp)) trt<-which(dic==trtgrp)




  ###### Estimate propensity scores, could be external ps.estimate or estimated through formula ##############################################
  if(is.null(ps.estimate)){

    ####estimate ps######
    e.h<-do.call(PSmethod_sga,c(list(ps.formula = ps.formula, method=method, data=data, weight=weight)))
    if(method =='LASSO' ){
      nonzero_coef <- e.h$nonzero_coef
      e.h <-e.h$e.h
    }
    colnames(e.h)<-dic

  }else{
    #####imported ps######
    #the name for the propensity score
    if(length(ps.estimate)==n){
      if(trt==2){
        e.h<-c(ps.estimate)
      }else{
        e.h<-c(1-ps.estimate)
      }
    }else if(!setequal(colnames(ps.estimate),dic)){
      ps.estimate<-as.matrix(ps.estimate)
      e.h<-c(ps.estimate[,2])
      warning("wrong column name set for ps.estimate, treatment set as: ",dic[1], " , ", dic[2])
    }else{
      ps.estimate<-ps.estimate[,match(dic,colnames(ps.estimate))]
      e.h<-c(ps.estimate[,2])
    }
    e.h <- cbind(1-e.h,e.h)
    colnames(e.h)<-dic
  }

  #weight entropy needs extra clipping
  if(weight=="entropy"){
    e.h<-pmax(e.h,1e-6)
    e.h<-pmin(e.h,1-1e-6)
  }

  ######## Estimate weight ###################################################################################################

  #tilting function
  ftilt<-tiltbin(weight = weight)

  #evaluated tilting function
  tilt.h<-ftilt(c(e.h[,trt]))

  allwt<-(1/e.h)*tilt.h
  wt<-rep(0,n)
  for(i in 1:ncate){
    wt[znum==i]<-allwt[znum==i,i]
  }

  #also generate the baseline weight
  wt_bs<-rep(1,n)


  ######## Prepare design matrix for balance check#############################################################################

  if(is.null(xname)){
    # the design matrix for balance check
    covM<-as.data.frame(model.matrix(formula(ps.formula),data))
    covM<-as.matrix(covM)
    if (ncol(covM)>1){
      #drop intercept
      if(unique(covM[,1])==1){
        covM<-covM[,-1,drop=F]
      }
    }
  }else{
    ps.form_m <- paste(zname,"~",paste0(xname, collapse=" + "))
    covM<-as.data.frame(model.matrix(formula(ps.form_m),data))
    covM<-as.matrix(covM)
    if (ncol(covM)>1){
      #drop intercept
      if(unique(covM[,1])==1){
        covM<-covM[,-1,drop=F]
      }
    }
  }

  ncov <- dim(covM)[2]+1   #number of covariate, plus baseline


  ######## Prepare subgroup matrix for balance check#############################################################################
  submatrix<-as.matrix(data[,subgroup,drop=FALSE])

  # Access the overall and subgroup asd, all level all_sga_level, by expanding subgroup by their levels
  all_sga_level <- list()
  for (k in 1:ncol(submatrix)){
    all_sga_level[[k]]<- table(submatrix[,k])
  }
  # number of subgroup levels by each subgrouping variable

  # subgroup levels, how namy level each subgroup
  nlevel <- c(unlist(lapply(all_sga_level,length)))
  nsubg <- c(unlist(all_sga_level))


  nv <- length(nsubg)   # number of subgroups

  subg.name_raw <- names(nsubg) #raw subgroup levels


  subcovnames <- colnames(submatrix) #raw subgroup names

  subg.name <- paste0(rep(subcovnames,nlevel),'=',subg.name_raw)#combine variable name to levels

  names(all_sga_level) <- subcovnames

  names(nsubg)<-subg.name #assign new subgroup names to summary table



  ######## the asd function to calculate the balance##########################################################################
  abs_stand_diff <- function(x_j, z, w){
    # Inputs:
    # x_j: a vecor of the covariate
    # z: a vector of treatment indicator
    # w: a vector of weights
    if (anyNA(w)) {return (NA)} else
      x_j <- as.numeric(x_j)

    absdiff <- abs(sum(x_j*z*w)/sum(z*w) - sum(x_j*(1-z)*w)/sum((1-z)*w))
    tstat <- absdiff/sqrt((var(x_j[which(z==1)])+var(x_j[which(z==0)]))/2)
    return (tstat)
  }


  overall_asd <- apply(covM, 2, abs_stand_diff, z, wt) #weighted

  overall_asd_bs <- apply(covM, 2, abs_stand_diff, z, wt_bs) #baseline

  #Calculate balance per subgroup across covariates
  groups_asd <- c() #weighted
  groups_asd_bs <- c() #baseline
  names_col <- c()

  for(r in 1:ncol(submatrix)){
    level.name <- names(all_sga_level[[r]])
    for(g in 1:length(level.name)){
      find_g <- which(submatrix[,r]==level.name[g])
      g_asd <- apply(covM[find_g, ], 2, abs_stand_diff, z[find_g], wt[find_g])
      groups_asd <- cbind(groups_asd, g_asd)

      g_asd_bs<-apply(covM[find_g, ], 2, abs_stand_diff, z[find_g], wt_bs[find_g])
      groups_asd_bs<- cbind(groups_asd_bs, g_asd_bs)

    }
  }


  colnames(groups_asd) <-subg.name
  colnames(groups_asd_bs) <-subg.name


  ASD <- cbind(Overall=overall_asd,groups_asd) #weighted

  ASD_bs <- cbind(Overall=overall_asd_bs,groups_asd_bs) #baseline


  ############## calculate VIF, Effective sample size  ###########################################################################
  # weights approximation

  #Matrix for effective sample size
  eff.sample.size<-matrix(-1,nrow=ncate*(nv+1),ncol=length(weight)) #ncate: 2trt  nv:number of levels plus overall
  colnames(eff.sample.size)<-c(weight) #given a single weight, this is a long matrix

  rownames(eff.sample.size)<-c(t(outer(paste0(colnames(ASD),'_treatment_'),dic,paste0))) #subgroupname by treatment

  eff.sample.size_bs<-eff.sample.size #also calculate the baseline

  #function to calculate eff and vif
  eff_vif_cal<-function(wt=1,ztmp){
    eff_est<-c()
    nj<-0
    for (j in 1:ncate){
      eff_est<-c(eff_est,sum(wt*(ztmp==j))^2/sum((wt*(ztmp==j))^2))
      nj<-nj+1/sum(ztmp==j)
    }
    vif_est<-c(1/nj*sum(1/eff_est))
    return(c(eff_est,vif_est))
  }

  vif <- c()
  vif_bs<-c()

  #overall value
  overall_val<-eff_vif_cal(wt,znum) #use the numeric z (1/2)
  vif<-c(vif,overall_val[ncate+1])
  eff.sample.size[1:ncate,]<-overall_val[1:ncate]

  overall_val_bs<-eff_vif_cal(wt_bs,znum)
  vif_bs<-c(vif_bs,overall_val_bs[ncate+1])
  eff.sample.size_bs[1:ncate,]<-overall_val_bs[1:ncate]

  countertmp<-2 #temp counter to assign the index of vif and eff
  for(r in 1:ncol(submatrix)){
    level.name <- names(all_sga_level[[r]])
    for(g in 1:length(level.name)){
      find_g <- which(submatrix[,r]==level.name[g])

      overall_tmp<-eff_vif_cal(wt[find_g],znum[find_g])
      vif<-c(vif,overall_tmp[ncate+1])
      eff.sample.size[(1+countertmp):(2+countertmp),]<-overall_tmp[1:ncate]

      overall_tmp_bs<-eff_vif_cal(wt_bs[find_g],znum[find_g])
      vif_bs<-c(vif_bs,overall_tmp_bs[ncate+1])
      eff.sample.size_bs[(1+countertmp):(2+countertmp),]<-overall_tmp_bs[1:ncate]

      countertmp<-countertmp+2

    }
  }

  names(vif)<-colnames(ASD)
  names(vif_bs)<-colnames(ASD)

  nsubg <- c(length(z), nsubg)   # subgroup sample size including overall

  names(nsubg)[1]<-'Overall'
  if(method=="LASSO" & is.null(ps.estimate)){
    output<-list(trtgrp=trtgrp, propensity=e.h, ps.weights= weight, ASD=ASD, ASD_bs=ASD_bs, vif=vif,vif_bs=vif_bs, nsubg=nsubg, ess=eff.sample.size,ess_bs=eff.sample.size_bs,subgroup=subgroup, method=method, nonzero_coef=nonzero_coef)

  }else{

  output<-list(trtgrp=trtgrp, propensity=e.h, ps.weights= weight, ASD=ASD, ASD_bs=ASD_bs, vif=vif,vif_bs=vif_bs, nsubg=nsubg, ess=eff.sample.size,ess_bs=eff.sample.size_bs,subgroup=subgroup, method=method)
  }
  class(output)<-"SumStat_sga"
  output

}


