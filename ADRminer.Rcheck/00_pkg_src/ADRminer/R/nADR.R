#' @encoding UTF-8
#' @title Number of adverse drug reaction associated with adverse event drug pair
#' @description \code{nADR} allows the user to get the number of spontaneous reports associated with an ae - drug pair. It works on both pvInd and pvCont object.
#' @name nADR
#' @param object an object of class pvInd or pvCont
#' @param drugId the label of the drug of interest (character)
#' @param aeId the label of the adverse event of interest (character)
#' @param ... not of any use for now.
#' @author Ismaïl Ahmed
#' @return a pvInd object
#' @keywords pvInd
#' @docType methods
#' @export
#' @aliases nADR.pvInd nADR.pvCont
#' @usage
#' \method{nADR}{pvCont}(object, drugId, aeId, \dots)
#' \method{nADR}{pvInd}(object, drugId, aeId, \dots)

# nADR definition ----------------------------------------------------------
nADR <- function (object, drugId, aeId, ...) UseMethod("nADR")

# nADR pvInd----------------------------------------------------------------
nADR.pvInd <- function(object, drugId, aeId, ...){
  if (!is.character(drugId)) stop("drugId must be a character")
  if (length(drugId)>1) stop("For now, nADR only handle one drug at a time")
  if (!is.character(aeId)) stop("aeId must be a character")
  if (length(aeId)>1) stop("For now, nADR only handle one ae at a time")
  
  idxDrug <- match(drugId, colnames(object@drug))
  if (is.na(idxDrug)) stop("drugId doesn't match any of the drug label in object")
  idxAe <- match(aeId, colnames(object@ae))
  if (is.na(idxAe)) stop("aeId doesn't match any of the ae label in object")
  
  nADR <- as.numeric(t(object@ae[,idxAe]) %*% object@drug[,idxDrug])
  nADR    
}
  
# nADR pvCont----------------------------------------------------------------

nADR.pvCont <- function(object, drugId, aeId, ...){
  if (!is.character(drugId)) stop("drugId must be a character")
  if (length(drugId)>1) stop("For now, nADR only handle one drug at a time")
  if (!is.character(aeId)) stop("aeId must be a character")
  if (length(aeId)>1) stop("For now, nADR only handle one ae at a time")
  
  
  presDrug <- drugId %in% object@drugLab
  if (!presDrug) stop("drugId doesn't match any of the drug label in object")
  presAe <- aeId %in% object@aeLab
  if (!presAe) stop("aeId doesn't match any of the ae label in object")
  
  idxDrug <- object@drugLab == drugId
  idxAe <- object@aeLab == aeId
  idxADR <- idxDrug & idxAe
  if(!any(idxADR)) return(0)
  
  nADR <- as.numeric(object@n[idxADR,1])
  nADR
    
}
