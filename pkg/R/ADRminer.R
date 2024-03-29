#' ADRminer: Adverse Drug event Reporting system miner
#' 
#' @docType package
#' @name ADRminer
# @import LBE
#' @import Matrix 
#' @import glmnet
#' @import methods
#' @import data.table
#' @description \code{ADRminer} is dedicated to the automated generation of drug safety signals from spontaneous reporting databases. While currently under development, this package will likely replace the PhViD package. Great efforts were made to re-think the structure of the objects manipulated. In particular, ADRminer handles individual spontaneous reports which makes it possible to conduct stratified analyses according to several covariates. From a more technical point of view, objects manipulated (\code{\link{pvInd}} and \code{\link{pvCont}}) are of class S4 and rely on the sparse matrix representation (Matrix package) in order to save as much as random access memory as possible. This package also contains an lasso based detection strategy recently which relies on the very efficient glmnet package and exploits parallel computing.
# @details For the frequentist methods, the package requires the LBE procedure that is stored in the Bioconductor website \url{http://bioconductor.org/}. LBE can be installed by entering \code{source("http://bioconductor.org/biocLite.R")} \code{biocLite("LBE")} in the R console.
# @import LBE Matrix glmnet methods
NULL