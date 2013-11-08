#' @title pvInd
#' @name pvInd
#' @aliases is.PvInd
#' @usage 
#' pvInd(drugMat, aeMat, covMat=NULL, dMarginLim=1, aeMarginLim=1, typeCompact=TRUE)
#' 
#' is.PvInd(x)
#' @param drugMat is a matrix of drugs
#' @param aeMat is a matrix of  the adverse events
#' @param covMat is a matrix containing a sex, age and other variables
#' @param dMarginLim is a limit margin of drugs
#' @param aeMarginLim is a limit margin of the adverse events
#' @param typeCompact ... 
#' @description Function that converts some variables into an pvIndobject that can be used in the signal detection method functions.
#' @author Youness Ergaibi & Ismaïl Ahmed
#' @return  return PvInd object
#' @keywords pvInd


# pvInd method definition ---------------------------------------------------- 
pvInd <- function(drugMat, aeMat, covMat=NULL, dMarginLim=1, aeMarginLim=1, typeCompact=TRUE){
  if(nrow(drugMat)==0 | nrow(aeMat)==0 | nrow(covMat)==0) { stop("drugMat and ae must contain at least one row")}
  ########### transfor drugMat  et aeMat en sparse Mat #################
  if(typeCompact==TRUE){
    
    #on transforme la matrice drugMat en matrice sparse
    colnD <- levels(factor(drugMat[,2]))
    rownD <- levels(factor(drugMat[,1]))
    idxR <- as.numeric(factor(drugMat[,1]))
    idxC <- as.numeric(factor(drugMat[,2]))
    
    sD <- sparseMatrix(idxR,idxC, x=1)
    sD[sD>1] <- 1
    rownames(sD) <- rownD
    colnames(sD) <- colnD
    #on transforme la matrice aeMat en matrice sparse
    colnD <- levels(factor(aeMat[,2]))
    rownD <- levels(factor(aeMat[,1]))
    idxR <- as.numeric(factor(aeMat[,1]))
    idxC <- as.numeric(factor(aeMat[,2]))
    sAe <- sparseMatrix(idxR,idxC, x=1)
    sAe[sAe>1] <- 1
    rownames(sAe) <- rownD
    colnames(sAe) <- colnD
    
  }
  else {
    #on transforme la matrice drugMat en matrice sparse
    rownD<-rownames(drugMat)
    colnD<-colnames(drugMat)
    sD<-Matrix(drugMat,sparse=T)
    sD[sD>1] <- 1
    rownames(sD)<-rownD
    colnames(sD)<-colnD
    
    #on transforme la matrice aeMat en matrice sparse
    rownD<-rownames(aeMat)
    colnD<-colnames(aeMat)
    sAe<-Matrix(aeMat,sparse=T)
    sAe[sAe>1] <- 1
    rownames(sAe)<-rownD 
    colnames(sAe)<-colnD
  }
  sex <- NULL
  age <- NULL
  date <- NULL
  
  if (!is.null(covMat)){
    idx1 <- match(row.names(sAe), covMat[,1])
    if(sum(is.na(idx1))) stop("labels between Ae matrix and covariate matrix missmatch")
    covMatSort <- covMat[idx1,]
    
    id <- which(colnames(covMatSort) =="sex")
    sex <- if(length(id)==1) as.factor(covMatSort[,id])
    id<-which(colnames(covMatSort)=="age")
    age <- if(length(id)==1) as.numeric(covMatSort[,id])
    #id<-which(colnames(covMatSort)=="date")
    #date <- if(length(id)==1) covMatSort[,id]
  }
  nSp <- nrow(sD)
  dMargin <- drop(Matrix(1, nrow=1, ncol=nSp) %*% sD)
  aeMargin <- drop(Matrix(1, nrow=1, ncol=nSp) %*% sAe)
  
  PvInd<-new(Class="PvInd", drug=sD, ae=sAe,age=age, sex=sex, date=NULL, dMargin=dMargin, aeMargin=aeMargin, cov=covMatSort)  
  return(PvInd)
}

