#' @encoding UTF-8
#' @title logistic regressions for pvInd objetcs
#' @name pvLogistic
#' @description pvLogistic implements a detection strategy based on the use of the logistic regressions. One major advantage is that it makes it possible to easily account for covariates without requiring using Mantel Haenzel approximates. This function will compute odds ratios, pvalues and for each adverse event count pairs. Since it deals with individuals data (instead of cell counts), this function is much more computationnaly intensive. This amount of time can be greatly decreased by the setup of parallel package which makes it possible to use multiple cores of the computers. By contrast with pvPen, this function is not particularly RAM demanding so it can be run on standard desktop machines.
#'  This function heavily relies on the highly efficient \code{glmnet} package. \code{pvPen} is computationnaly highly demanding. It is also strongly advised to use a linux or a mac computer in order to use several cores. Note also that this regression approach should be used with ae and drugs both having a reasonable number of reports. 
#' @param object an object of class PvInd
#' @param cov a character vector indicating which covariate have to be used in the analyses. 
#' Such variables are first use in a logistic regression model and the residual is then used as the outcome 
#' in the penalized regression
#' @param aeLab The label of the adverse event to be regressed. By default, the \code{pvPen} will regress all ae, one at a time and this can be very time consuming and RAM demanding. 
#' @param detectCriter Can be either BIC or AIC. These criteria are used to select the final model.
#' @param posConst If TRUE (default), the regression coefficients are constrained to be non non negative, this to ensure that the final model only contains drugs "increasing" the risk of a given ae.

NA