#' Plot the distribution of propensity scores and subgroup balance statistics via Connect-S plot
#'
#' Summarize the SumStat_sga object, plot the subgroup balance statistics under weighting versus no weighting.
#' In a Connect-S plot (Yang et al. 2021), each row represents a subgroup and each column represents a confounder;
#' each dot is shade-coded according to the value of the ASMD corresponding to the specific subgroup and confounder.
#' Dark gray and black dots in the Connect-S plot flag meaningful covariate imbalance. The last two columns show the
#' subgroup sample size and estimated variance inflation due to weighting.
#'
#' @param x a \code{SumStat_sga} object obtained with \code{\link{SumStat}} function.
#' @param varlist an optional vector specifying the variables to be plotted on the x-axis by index or variable names.
#'  The default plots all main effect variables specified by \code{ps.formula} in the \code{SumStat_sga} function,
#'  excluding subgrouping variables.
#' @param base an indicator to specify whether to print the unadjusted ASMDs before weighting. The default is FALSE.
#' @param plotsub an indicator to specify whether to treat subgrouping variables as confounders and print them on the x-axis.
#'  The default is FALSE.
#' @param ... further arguments passed to or from other methods.
#'
#' @details For the Connect-S plot, the shade is categorized based on the common ASMD threshold of 0.1 and 0.2
#' following Austin and Stuart (2015), with darker shade implying more severe imbalance.
#'
#' @return Plot of the indicated type.
#'
#' @references
#' Yang, S., Lorenzi, E., Papadogeorgou, G., Wojdyla, D. M., Li, F., & Thomas, L. E. (2021).
#' Propensity score weighting for causal subgroup analysis. Statistics in medicine, 40(19), 4294-4309.
#'
#' Austin, P.C. and Stuart, E.A. (2015). Moving towards best practice when using inverse probability of treatment weighting (IPTW) using the propensity score to estimate causal treatment effects in observational studies.
#' Statistics in Medicine, 34(28), 3661-3679.
#'
#' @export
#'
#'
#' @import ggplot2
#' @importFrom  stats binomial coef cov formula glm lm model.matrix plogis poisson predict qnorm quantile sd
#' @importFrom  utils capture.output combn
#' @importFrom  graphics hist legend par
#' @importFrom  grDevices rgb
#' @importFrom rlang .data

plot.SumStat_sga <- function(x,varlist=NULL,base=FALSE,plotsub=FALSE,...){
  VarName  <- Subgroups <- lcolor <- c()
  vif<-x$vif
  ASD<-x$ASD
  subgroup<-x$subgroup

  if(base){
    vif<-x$vif_bs
    ASD<-x$ASD_bs
  }
  nsubg<-x$nsubg
  VIFW <- round(vif,2)


  if(!is.null(varlist)){
    ASD<-ASD[varlist,]
  }else{
    ASD<-ASD[!grepl('\\:',rownames(ASD)),]
    #print not subgroup as confounder
    if(!plotsub){
      namelist<-c()
      for (i in subgroup){
        namelist<-c(namelist,which(i==rownames(ASD)))
      }
      namelist<-unique(namelist)
      ASD<-ASD[-namelist,]
    }
  }


  # connect-S plot
  nv<-length(colnames(ASD))
  ncov<-length(rownames(ASD))
  mydata <- data.frame(ASMD = c(ASD),
                       Subgroups = rep(colnames(ASD), each=nrow(ASD)),
                       VarName = rep(rownames(ASD), nv))

  mydata$VarName <- factor(mydata$VarName, levels=unique(mydata$VarName))
  mydata$Subgroups <- factor(mydata$Subgroups, levels=unique(mydata$Subgroups))

  # toplot<- plyr::mutate(mydata, lcolor = ifelse(ASMD <=0.1 |is.na(ASMD), "0",
  #                                         ifelse((ASMD <=0.15), "1",
  #                                                ifelse((.data$ASMD <=0.20), "2","3"))))

  mydata$lcolor <- ifelse(mydata$ASMD <=0.1 |is.na(mydata$ASMD), "0",
                          ifelse((mydata$ASMD <=0.15), "1",
                                 ifelse((mydata$ASMD <=0.20), "2","3")))
  toplot <- mydata
  toplot$lcolor <- factor(toplot$lcolor, levels=c("0","1","2","3"))

  group.colors <- c("0" = "White", "1" = "grey85", "2" ="grey50","3"="black")

  g1 <- ggplot(toplot, aes(x =VarName , y =Subgroups , fill=lcolor)) +
    geom_dotplot(binaxis = "y",binwidth = 0.75, stackdir = "center")+labs(x="Covariate", y="Subgroup"  )+
    scale_fill_manual(values=group.colors,name="ASMD",  labels=c("<0.1", "0.1-0.15", "0.15-0.20", ">0.20"), drop=FALSE)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    theme(axis.text.x = element_text( angle=-45,vjust = 0, hjust=0.3),axis.text = element_text( size = 14 ),
          axis.title = element_text( size = 20, face = "bold" ),legend.position="bottom",legend.title = element_text(size=18),
          legend.text=element_text(size=18), legend.key.size = unit(2,"line"))+
    annotate("text", x = rep(ncov+1,nv), y = seq(1, nv,1), label = format(nsubg,zero.print = T),size = 5)+
    annotate("text", rep(ncov+2,nv), y = seq(1, nv,1), label = format(VIFW,zero.print = T),size = 5)+
    annotate("text", x = c(ncov+1,ncov+2), y = rep(nv+1,2), label = c("Size", "VIF"),size = 5)+
    expand_limits(x= c(0, length(levels(toplot$VarName)) + 3),y= c(0, length(levels(toplot$Subgroups)) + 1.5) )

  splot <- g1+ theme(panel.border = element_rect(linetype = "solid", fill = NA))

  print(splot)
}


