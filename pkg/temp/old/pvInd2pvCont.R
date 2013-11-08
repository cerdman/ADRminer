#' @title \code{pvInd2pvCont} function
#' @param x an object of class PvInd
#' @description The function pvInd2pvCont converts pvInd into pvCont.
#' @author Youness Ergaibi & Ismaïl Ahmed
#' @return  A pvCont object. 
#' @keywords pvCont

pvInd2pvCont <- function(x) {
  if (!is(x, "PvInd")) stop("x must be an object of class PvInd")
  data<-t(x@drug)%*%x@ae
  if(nrow(data)==1)stop("not enough observations")
  coord <- which(data!= 0, arr.ind=TRUE)
  coord <- coord[order(coord[,1]),] 
  labeldrug<- rownames(data)[coord[,1]] # on conserve les libellés des médicaments qui restent
  labeldrug<-as.factor(labeldrug)
  labelae <- colnames(data)[coord[,2]] # on conserve les libellés des effets indésirables qui restent
  labelae<-as.factor(labelae)
  n11 <- data[coord] 
  obj<- pvCont(labeldrug,labelae,n11)
  dMargin<-obj@dMargin
  aeMargin<-obj@aeMargin
  expN<-obj@expN
  
  PvCont<-new(Class="PvCont", drugLab=labeldrug, aeLab=labelae,n=n11, dMargin=dMargin, aeMargin=aeMargin, expN=expN)
  return(PvCont)
}

