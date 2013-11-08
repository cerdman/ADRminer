#' @encoding UTF-8
#' @title pvIndResize
#' @description \code{pvIndResize} makes it possible to remove drug(s) and/or ae(s) whose margin are less
#' than a given number.
#' @param pvIndObj an object of class "pvInd" (see \code{\link{pvInd-class}} and \code{\link{pvInd}})
#' @param aeMarginMin Minimum number of reports in which an ae must be involved in order to be kept
#' @param dMarginMin Minimum number of reports in which a drug must be involved in order to be kept
#' @return an object of class pvInd
#' @author Ismaïl Ahmed
#' @export

pvIndResize <- function(pvIndObj, aeMarginMin=1, dMarginMin=1){
  pvInd <- pvIndObj
  if (!inherits(pvInd, "pvInd")) stop("object must be of class pvInd")
  if (aeMarginMin > max(pvInd@aeMargin)) stop("The value specified for aeMarginMin exceeds the number of reports of the most frequent ae ")
  if (dMarginMin > max(pvInd@dMargin)) stop("The value specified for dMarginMin exceed the number of reports of the most frequent drug")
  if (aeMarginMin>1){
    idxAe <- which(pvInd@aeMargin < aeMarginMin)
    pvInd@ae <- pvInd@ae[,-idxAe]
    pvInd@aeMargin <- pvInd@aeMargin[-idxAe]
  }
  if (aeMarginMin>1){
    idxD <- which(pvInd@dMargin < dMarginMin)
    pvInd@drug <- pvInd@drug[,-idxD]
    pvInd@dMargin <- pvInd@dMargin[-idxD]
  }
  return(pvInd)
}