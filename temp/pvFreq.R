#' @encoding UTF-8
#' @title Fisher's exact test, Reporting Odds Ratio & False Discovery Rate
#' @name pvOR
#' @param object An object of class pvInd or pvCont
#' @param method Either 
#' @param stat Statistic used to order the drug-events pairs: 
#' \itemize{
#' \item pvalue: P-value 
#' \item Lb95: Lower bound of the 95\% two sided confidence interval of log(ror).
#' }
#' @param detectCriter Decision rule for the signal generation based on
#' \itemize{
#' \item FDR: a prespecified value (\code{detectCriter}) for the false discovery rate (Default value) 
#' \item nbSig: a prespecified number (\code{detectCriter}) of signals
#' \item assocMeasure: a prespecified value (\code{detectCriter}) for the Ranking statistic (\code{assocMeasure})
#' }
#' @param  criterThres Threshold for \code{detectCriter}. Ex 0.05 for FDR.
#' @param nMin Minimum number of spontaneous reports for a drug-event pair to be potentially considered as a signal. By default, \code{nMin}=1.
#' @param strat Character vector containing the name of the covariates to be used for a stratified analysis
#' @param allRes Logical value indicating whether all drug-event combination results must be provided.
#' @description Reporting Odds Ratio extended to the multiple comparison framework.
#' @author Youness Ergaibi & Ismaïl Ahmed
#' @return allSig Data.frame summarizing the results for all drug-event combinations with at least \code{nMin} spontaneous reports ordered by \code{assocMeasure}. This output is provided if \code{allRes=TRUE}. Operating characteristics are estimated according to Ahmed et al (Biometrics 2010).
#' @return sig Same as \code{allSig} but restricted to the list of generated signals.
#' @return nSig Number of generated signals.
#' @return call Arguments entered when calling \code{ror}.
#' @keywords ror 


setGeneric(
  name="ror",
  def=function(object, ...){standardGeneric("pvFreq")}
)
