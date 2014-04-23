#' @encoding UTF-8
#' @title ADRminer: Adverse Drug event Reporting systems miner
#' @author Ismaïl Ahmed, Antoine Poncet and Youness Ergaibi (developer)
#' @docType package
#' @name ADRminer
#' @description \code{ADRminer} is dedicated to the automated generation of drug safety signals from spontaneous reporting databases. While currently under development, this package will likely replace the PhViD package. Great efforts were made to re-think the structure of the objects manipulated. In particular, ADRminer handles individual spontaneous reports and thus makes it  possible to conduct stratified analyses according to several covariates. From a more technical point of view, objects manipulated (\code{\link{pvInd}} and \code{\link{pvCont}}) are of class S4 and rely on the sparse matrix representation (Matrix package) in order to save as much as random access memory as possible. 
#' @details For the frequentist methods, the package requires the LBE procedure that is stored in the Bioconductor website \url{http://bioconductor.org/}. LBE can be installed by entering \code{source("http://bioconductor.org/biocLite.R")} \code{biocLite("LBE")} in the R console.
#' @import LBE Matrix glmnet methods
#' @references
#' Ahmed I, Thiessard F, Miremont-Salamé G, Bégaud B, Tubert-Bitter P. Pharmacovigilance data mining with methods based on false discovery rates: a comparative simulation study. Clin. Pharmacol. Ther. 2010 Oct;88(4):492-498. 
#' 
#' Ahmed I, Dalmasso C, Haramburu F, Thiessard F, \enc{Broët}{Broet} P, Tubert-Bitter P. False discovery rate estimation for frequentist pharmacovigilance signal detection methods. Biometrics. 2010 Mar;66(1):301-309. 
#' 
#' Ahmed I, Haramburu F, Fourrier-Réglat A, Thiessard F, Kreft-Jais C, Miremont-Salamé G, Bégaud B, Tubert-Bitter P. Bayesian pharmacovigilance signal detection methods revisited in a multiple comparison setting. Stat Med. 2009 Jun 15;28(13):1774-1792.
#'  
#' Bate A, Lindquist M, Edwards IR, Olsson S, Orre R, Lansner A, De Freitas RM. A Bayesian Neural Network Method for Adverse Drug Reaction Signal Generation European Journal of Clinical Pharmacology, 1998, 54, 315-321.  
#' 
#' Dalmasso C, Broet P, Moreau T (2005), A simple procedure for estimating the false discovery rate, Bioinformatics, Bioinformatics, 21: 660 - 668. 
#' 
#' DuMouchel W. Bayesian Data Mining in Large Frequency Tables, with an Application to the FDA Spontaneous Reporting System. The American Statistician. 1999, 53. 177.

.onAttach=function(libname, pkgname){ 
  packageStartupMessage("Loaded ADRminer ", as.character(packageDescription("ADRminer")[["Version"]]),"\n")
}