#' @name pvInd2pvContStrat
#' @title pvInd2pvContStrat
#' @param PvInd a pvInd object
#' @param strat Character vector containing the name of the covariates to be used for a stratified analysis
#' @description pvInd2pvContStrat converts a pvInd object into a list of pvCont objects for each combination of \code{strat}
#' @author Youness Ergaibi & Ismaïl Ahmed
#' @return A list of pvCont objects.  
#' @keywords pvCont

pvInd2pvContStrat <- function(PvInd, strat)
{
  cont <- t(PvInd@drug)%*%PvInd@ae
  coord <- which(cont!= 0, arr.ind=TRUE)
  coord <- coord[order(coord[,1]),] 
  
  idx <- match(strat, c(colnames(PvInd@cov)))
  
  if (sum(is.na(idx)>0)) stop("at least one of the stratification covariates does not exist in the object")
  varStrat <- factor(apply( PvInd@cov[,idx, drop=F], 1, paste, collapse=" "))
  L <- nlevels(varStrat)
  #print(table(varStrat))
  PvContList <- vector("list",L)
  aeSplit <- vector("list", length=nlevels(varStrat))
  drugSplit <- vector("list", length=nlevels(varStrat))
  for (i in 1: nlevels(varStrat)){
    idx <- varStrat == levels(varStrat)[i]
    if(length(idx !=0)){
      drugSplit[[i]] <- PvInd@drug[idx, ,drop=F]
      aeSplit[[i]] <- PvInd@ae[idx, ,drop=F]
    }
  }
  
  labeldrug <- vector("list", L)
  labelae <- vector("list", L)
  n11 <- vector("list", L)
  for(l in 1:L)
  {
    cont <- t(drugSplit[[l]])%*%aeSplit[[l]]
    labeldrug[[l]] <- as.factor(rownames(cont)[coord[,1]])
    labelae[[l]] <- as.factor(colnames(cont)[coord[,2]])
    n11[[l]] <- cont[coord] 
    #obj <- pvCont(labeldrug[[l]],labelae[[l]], n11[[l]])
    n1._mat <- apply(cont,1,sum) # on recalcule les marges des lignes...
    n.1_mat <- apply(cont,2,sum)
    dMargin <- n1._mat[coord[,1]] # on affecte à chaque notification sa marge "ligne"...
    aeMargin <- n.1_mat[coord[,2]] # ... et sa marge "colonne"
    N <- sum(n11[[l]]) # le nombre total de notifications
    expN <- dMargin*aeMargin/N
    PvContList[[l]] <- new("PvCont",drugLab=labeldrug[[l]], aeLab=labelae[[l]],n=n11[[l]], dMargin=dMargin,  aeMargin= aeMargin,expN=expN)
    names(PvContList)[[l]] <- levels(varStrat)[l]
  }
  return(PvContList)                    
}