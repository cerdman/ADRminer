#' \code{setMargin} function
#' @param PvInd a pvInd object
#' @param dMarginLim minimum number of spontaneous reports for a drug to be included in the PvInd object
#' @param aeMarginLim minimum number of spontaneous reports for an adverse evebt to be included in the PvInd object
#' @description setMargin creates a pvInd object removing drugs and adverse events involved in less than dMarginLim and aeMargineLim.
#' @author Youness Ergaibi & Ismaïl Ahmed
#' @return A pvInd object
#' 
setMargin <- function(PvInd, dMarginLim=1, aeMarginLim=1)
{
  dMargin <- apply(PvInd@drug,2, sum)
  dIdx <- which(dMargin >= dMarginLim)
  PvInd@drug <- PvInd@drug[,dIdx]
  dMargin<-apply(PvInd@drug,2, sum)
  aeMargin<-apply(PvInd@ae,2, sum)
  aeIdx <- which(aeMargin >= aeMarginLim)
  PvInd@ae <- PvInd@ae[,aeIdx]
  aeMargin<-apply(PvInd@ae,2, sum)
  PvInd<-new(Class="PvInd", drug=PvInd@drug, ae=PvInd@ae, dMargin=dMargin, aeMargin=aeMargin)
  return(PvInd)  
}
