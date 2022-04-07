# PSweight.sga: Propensity Score Weighting for Causal Subgroup Analysis

Performs propensity score weighting for causal subgroup analysis of observational studies and randomized trials. Enables the estimation and inference of subgroup average causal effects with binary treatments using overlap weights (S-ATO), inverse probability of treatment weights (S-ATE), and average subgroup treatment effect among the treated weights (S-ATT). When the covariates and subgrouping variables are provided, this package allows to automatically perform the Post-LASSO to select important covariate-subgroup interactions and generate the Connect-S plot, introduced in Yang et al. (2021) <https://doi.org/10.1002/sim.9029 >

In the design stage, the SumStat_sga function is used to generate the estimated propensity scores and balance diagnostics after propensity score weighting. The print and plot functions are available to tabulate and plot the Connect-S plot introduced in Yang et al. (2021) to visualize weighted balance statistics in subgroups.

In the analysis stage, through the PSweight_sga function, the average potential outcomes for each subgroup and each treatment group is estimated using weighting, and the summary function generates point estimates, standard errors and confidence intervals for the desired subgroup causal contrasts of interest. 

For binary outcomes, both the additive and ratio estimands (causal relative risk and odds ratio) are considered. To allow for additional flexibility in specifying the propensity score and outcome models, the function can also work with user-supplied propensity score estimates. The variance is estimated by nonparametric bootstrap that ignores the variability in estimating these nuisances.

# How to install
The latest development version can be installed directly from Github using devtools:
``` Ruby
if (!require("devtools")) install.packages("devtools")
devtools::install_github(repo="siyunyang/PSweight.sga")
library(PSweight.sga)
```
# Link to Paper

You can access the paper with examples on real-world dateset at https://onlinelibrary.wiley.com/doi/full/10.1002/sim.9029

# Example
The following example uses a simulated observational data with two treatment arms to illustrate the unity of PSweight.sga The simulated dataset includes 3000 rows, with each row represents information recorded from each individual. There are 22 variables (columns). The treatment is the variable Treatment. The outcome of interest is variable Y. X1-X20 are pre-treatment covariates among which X1-X8 are binary, and X9-X20 are continuous.

``` Ruby
# library(devtools)
# install_github(repo="siyunyang/PSweight.sga")
library(PSweight.sga)


data(psdata_sga)
# pre-specify the confounder names: X1-X20
xname <- paste0('X',1:20)
# pre-specify subgroups of interest by column names
subgroup<-paste0('X',c(1,2,3, 19,20))


##### Design
# OW post-LASSO PS model
set.seed(123)
p <- SumStat_sga( subgroup=subgroup,xname=xname, zname="Treatment",data=psdata_sga, method='LASSO',weight="overlap")
summary(p)
p$nonzero_coef

# main effect logistic regression PS model
ps.form_m <- paste("Treatment~",paste0("X",1:20,collapse = "+"))
p <- SumStat_sga(ps.formula = ps.form_m, subgroup = subgroup, data=psdata_sga, method="glm",weight="overlap")
summary(p)

# user provided ps
e.h <- runif(3000,0,1)
p <- SumStat_sga(subgroup=subgroup,xname=xname,zname = "Treatment",ps.estimate = e.h, data=psdata_sga, weight="IPW")

# balance diagnostic
plot(p)
# the unweighted connect-s: 
plot(p, base = T)


##### Analysis
p1<-PSweight_sga(ps.formula=ps.form_m,subgroup=subgroup,yname="Y",data=psdata_sga,R=50,weight="overlap")
# can specify arbitrary contrast statement for subgroup causal effects
summary(p1,contrast = rbind(c(1,-1)))
# test of HTE across subgroup levels
summary(p1,het = T)

# to extract the summary statistics
# subgroup average potential outcomes by treatment arms
p1$muhat
# subgroup causal effects
tmp <-do.call(rbind, summary(p1)$estimates)
rownames(tmp) <- rownames(p1$muhat)

# OW-pLASSO
p2 <- PSweight_sga(xname=xname, subgroup=subgroup, yname="Y",zname="Treatment",data=psdata_sga,R=50, method='LASSO',weight="overlap")
summary(p2)

# binary outcome
psdata_sga$Y2 <- ifelse(psdata_sga$Y < 2, 0 ,1)

p <- PSweight_sga(ps.formula = ps.form_m, subgroup = subgroup,yname = "Y2", data=psdata_sga, weight="overlap")

summary(p)
summary(p,type='RR',contrast = rbind(c(1,-1),c(0.5,-1)))
summary(p,type='OR')
end
```
