#' Fitting propensity scores with different models
#'
#' The function \code{PSmethod_sga} is an internal function to estimate the propensity scores given a specified model through formula.
#' It is built into function \code{Sumstat_sga}, and \code{PSweight_sga}.
#'
#' @param ps.formula an object of class \code{\link{formula}} (or one that can be coerced to that class):
#' a symbolic description of the propensity score model to be fitted. Additional details of model specification
#' are given under "Details". This argument is optional if \code{ps.estimate} is not \code{NULL}.
#' @param method a character to specify the method for propensity model. When \code{ps.formula} is given, \code{"glm"} is the default; When \code{xname} is given, \code{"LASSO"} is the default.
#' @param weight a character or vector of characters including the types of weights to be used. \code{"IPW"} specifies the inverse probability weights for estimating the average treatment effect among the combined population (ATE). \code{"treated"} specifies the weights for estimating the average treatment effect among the treated (ATT). \code{"overlap"} specifies the (generalized) overlap weights for estimating the average treatment effect among the overlap population (ATO), or population at clinical equipoise. \code{"matching"} specifies the matching weights for estimating the average treatment effect among the matched population (ATM). \code{"entropy"} specifies the entropy weights for the average treatment effect of entropy weighted population (ATEN). Default is \code{"overlap"}.
#' @param data a data frame containing the variables in the propensity score model.
#'
#' @details A typical form for \code{ps.formula} is \code{treatment ~ terms} where \code{treatment} is the treatment
#' variable (identical to the variable name used to specify \code{zname}) and \code{terms} is a series of terms
#' which specifies a linear predictor for \code{treatment}. \code{ps.formula} by default specifies generalized
#' linear models given the default argument \code{method = "glm"}.  It fits the logistic regression. The argument
#' \code{method} allows user to choose model other than glm to fit the propensity score models. In Yang et al.(2021),
#' the \code{term} is suggested to include all main effects and pairwise subgroup-confounder interactions, combined with
#' the \code{method}=\code{"LASSO"} to select important interactions. We have included \code{"LASSO"} and \code{"gbm"},
#' through the \code{method} argument. Note that the current code does not handle multiple treatment groups.
#
#'
#' @return
#'
#' \describe{
#'
#' \item{\code{ e.h}}{a data frame of estimated propensity scores.}
#'
#' \item{\code{ nonzero_coef}}{the LASSO selected interactions when \code{method = "LASSO"}.}
#' }
#'
#' @references
#' Yang, S., Lorenzi, E., Papadogeorgou, G., Wojdyla, D. M., Li, F., & Thomas, L. E. (2021).
#' Propensity score weighting for causal subgroup analysis. Statistics in medicine, 40(19), 4294-4309.
#'
#' @export
#'
#' @import glmnet
#'




PSmethod_sga <-function(ps.formula=ps.formula, method="glm", weight='overlap', data=data){

  ps.formula<-as.formula(ps.formula)
  zname<-all.vars(ps.formula)[1]
  facz<-as.factor(data[,zname])
  #creat a dictionary for the original and recoded values in Z
  dic<-levels(facz)


  if (method=="glm"){
    ############## logistic #############################################################
    #change z to 0/1
    dataps<-data
    dataps[,zname]<- as.numeric(facz)-1

    fitglm <- glm(formula = ps.formula, data=dataps,family = binomial(link = "logit"))
    e.h <- fitglm$fitted.values
    e.h <- cbind(1-e.h,e.h)
  }else if(method =='LASSO'){
    ############## pLASSO #############################################################
    #change z to 0/1
    dataps<-data
    z<- as.numeric(facz)-1
    dataps[,zname]<-z

    fullmatrix <- model.matrix(ps.formula, data)

    # only penalize interactions
    # create an idx for penalty.factor
    idx <- rep(0, dim(fullmatrix)[2])
    idx[grep('\\:',colnames(fullmatrix))] <- 1

    fitLASSO <- cv.glmnet(y=factor(z), x=fullmatrix, penalty.factor=idx, family="binomial", maxit=50000)
    # nonzero_coef <- rownames(coef(fitLASSO, s='lambda.min'))[which(coef(fitLASSO, s='lambda.min')!=0)][-1]
    tmp_coeffs <- coef(fitLASSO, s = "lambda.min")
    nonzero_coef <- tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1]

    if(length(nonzero_coef)==0){
      warning("no coefficient ")
      fitpLASSO <- glm(factor(z)~1, family = binomial(link = "logit"))
    }else{
      fitpLASSO <- glm(factor(z)~., data=data.frame(fullmatrix[,nonzero_coef]), family = binomial(link = "logit"))
    }
    e.h <- fitpLASSO$fitted.values
    e.h <- cbind(1-e.h,e.h)
  }

  colnames(e.h)<-dic
  if(method =='LASSO'){ out=list(e.h=e.h, nonzero_coef=nonzero_coef)}else{
    out=e.h
  }
  return(out=out)
}





















