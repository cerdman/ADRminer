#' @encoding UTF-8
#' @title pvIndResize
#' @description \code{pvIndResize} makes it possible to select sets of drugs and/or aes. It also allows to remove drug(s) and/or ae(s) whose margin are less
#' than a given number.
#' @param pvIndObj an object of class "pvInd" (see \code{\link{pvInd-class}} and \code{\link{pvInd}})
#' @param drugId can be a vector of characters containing drug labels to be selected or a vector of column indices or the drug matrice (object$drug)
#' @param aeId can be a vector of characters containing ae labels to be selected or a vector of column indices or the ae matrice (object$ae)
#' @param aeMarginMin Minimum number of reports in which an ae must be involved in order to be kept
#' @param drugMarginMin Minimum number of reports in which a drug must be involved in order to be kept
#' @return an object of class pvInd
#' @author Ismaïl Ahmed
#' @export

pvIndResize <- function(pvIndObj, drugId=NULL, aeId=NULL, aeMarginMin=1, drugMarginMin=1){
  pvInd <- pvIndObj
  if (!inherits(pvInd, "pvInd")) stop("object must be of class pvInd")  
  if (length(drugId) > ncol(pvInd$drug)) stop("The length of drugId exceeds the number of drugs in pvIndObj ")
  if (length(aeId) > ncol(pvInd$ae)) stop("The length of aeId exceeds the number of aes in pvIndObj ")
  if (is.numeric(drugId) && (max(drugId)>ncol(pvInd$drug))) stop("Problem with drugId indices (at least one is too large)")
  if (is.numeric(aeId) && (max(aeId)>ncol(pvInd$ae))) stop("Problem with aeId indices (at least one is too large)")  
  if (aeMarginMin > max(pvInd@aeMargin)) stop("The value specified for aeMarginMin exceeds the number of reports of the most frequent ae ")
  if (drugMarginMin > max(pvInd@drugMargin)) stop("The value specified for drugMarginMin exceed the number of reports of the most frequent drug")
  if (!is.null(drugId)){
    if (is.character(drugId)) {
      idxDrug <- match(drugId, colnames(pvInd$drug))
      if (all(is.na(idxDrug))) stop("None of the DrugIds match with any of the drugs in the pvIndObj")
      if (sum(is.na(idxDrug))>0) warning("Some of the DrugIds do not match with any of the drugs in the pvIndObj")
      drugId <- drugId[!is.na(idxDrug)]
    }
    pvInd$drug <- pvInd$drug[,drugId]
    pvInd$drugMargin <- pvInd$drugMargin[drugId]
  }
  if (!is.null(aeId)){
    if (is.character(aeId)) {
      idxAe <- match(aeId, colnames(pvInd$ae))
      if (all(is.na(idxAe))) stop("None of the aeIds match with any of the drugs in the pvIndObj")
      if (sum(is.na(idxAe))>0) warning("Some of the aeIds do not match with any of the drugs in the pvIndObj")
      aeId <- aeId[!is.na(idxAe)]
    }
    pvInd$ae <- pvInd$ae[,aeId]
    pvInd$aeMargin <- pvInd$aeMargin[aeId]
  }    
  
  if (aeMarginMin>1){
    idxAe <- which(pvInd@aeMargin < aeMarginMin)
    if (length(idxAe)>0){
      pvInd@ae <- pvInd@ae[,-idxAe]
      pvInd@aeMargin <- pvInd@aeMargin[-idxAe]
    }
  }
  if (aeMarginMin>1){
    idxD <- which(pvInd@drugMargin < drugMarginMin)
    if (length(idxD)>0){
      pvInd@drug <- pvInd@drug[,-idxD]
      pvInd@drugMargin <- pvInd@drugMargin[-idxD]
    }
  }
  return(pvInd)
}