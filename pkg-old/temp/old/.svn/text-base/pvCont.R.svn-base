#' @name pvCont
#' @title pvCont
#' @aliases is.PvCont
#' @usage 
#' pvCont(drugLab, aeLab, n)
#' 
#' is.PvCont(x)
#' @param drugLab is a vector of drug labels
#' @param aeLab is a vector of adverse event labels
#' @param n is the number of spontaneous reports of the corresponding ae-d pair
#' @description \code{pvCont} converts a data.frame into a pvCont object that can be used in the signal detection method functions.
#' @author Youness Ergaibi & Ismaïl Ahmed
#' @return dMargin a Margin of drugs
#' @return aeMargin a Margin of adverse event 
#' @return expN an expected number 
#' @keywords pvCont
#' @rdname pvCont

# pvCont method definition ---------------------------------------------------- 
pvCont <- function(drugLab, aeLab, n)
{
  #tester si les 3 vecteurs ont meme longeur
  if(length(drugLab)!=length(aeLab))
    stop("les trois vecteurs n\'ont pas meme longeur")
  if(length(aeLab)!=length(n))
    stop("les trois vecteurs n\'ont pas meme longeur")
  #mettre chaque vecteur dans son prototype
  drugLab<-as.factor(drugLab)
  aeLab<-as.factor(aeLab)
  n<-as.numeric(n)
  data<-data.frame(drugLab,aeLab,n)
  names(data) <-c("drugLab","aeLab","n")
  ts<-xtabs(n~.,data=data, sparse = T)
  n1._mat <- apply(ts,1,sum) # on recalcule les marges des lignes...
  n.1_mat <- apply(ts,2,sum) # ...et des colonnes
  coord <- which(ts!=0, arr.ind=TRUE) 
  coord <-coord[order(coord[,1]),]
  n <- ts[coord] 
  N <- sum(n) # le nombre total de notifications
  dMargin <- n1._mat[coord[,1]] # on affecte à chaque notification sa marge "ligne"...
  aeMargin <- n.1_mat[coord[,2]] # ... et sa marge "colonne"
  expN<-dMargin*aeMargin/N
  PvCont<-new(Class="PvCont",drugLab=drugLab, aeLab=aeLab, n=n, dMargin= dMargin, aeMargin=aeMargin, expN=expN)
  return(PvCont)           
}
