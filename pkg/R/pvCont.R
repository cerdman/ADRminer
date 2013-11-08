#' @encoding UTF-8
#' @name pvCont
#' @title pvCont
#' @param object a pvInd object or a data.frame. The latter should have 3 columns, i.e. the drug labels, the ae labels and the corresponding number of spontaneous reports. 
#' @param strat Character vector containing the name of the covariates to be used for a stratified analysis.  \code{strat} is only meaningful when used along with a pvInd object.
#' @description pvCont converts an pvInd object into pvCont object. 
#' @author Ismaïl Ahmed
#' @return an object of class \code{\link{pvCont}}.
#' @exportMethod pvCont
#' @keywords pvCont
#' @docType methods
#' @aliases pvCont,pvInd,character-method
#' @aliases pvCont,pvInd,missing-method
#' @aliases names,pvCont
#' @aliases shows,pvCont
#' @aliases pvCont,data.frame,missing-method
#' @aliases pvCont,matrix,missing-method
setGeneric(
  name = "pvCont",
  def = function(object, strat){standardGeneric("pvCont")}
)

setMethod(
  f = "pvCont",
  signature = c(object="pvInd", strat="missing"),
  definition = function(object, strat){
    cont <- t(object@drug)%*%object@ae
    coord <- which(cont!= 0, arr.ind=TRUE)
    coord <- coord[order(coord[,1]),] 
    dLab <- rownames(cont)[coord[,1]]
    dLab <- as.factor(dLab)
    aeLab <- colnames(cont)[coord[,2]] 
    aeLab <- as.factor(aeLab)
    n11 <- matrix(cont[coord], ncol=1) 
    
    n1. <- as.numeric(cont %*% matrix(1, nrow=ncol(cont)))
    n.1 <- as.numeric(matrix(1, ncol=nrow(cont)) %*% cont)      
    dMargin <- matrix(n1.[coord[,1]], ncol=1)
    aeMargin <- matrix(n.1[coord[,2]], ncol=1)
    N <- sum(n11) 
    expN <- dMargin*aeMargin/N
    pvCont<-new(Class="pvCont", dLab=dLab, aeLab=aeLab, n=n11, dMargin=dMargin, aeMargin=aeMargin, expN=expN, N=N, strat=NULL, coord=coord)
    return(pvCont)
  }  
)

setMethod(
  f = "pvCont",
  signature = c(object="pvInd", strat="character"),
  definition = function(object, strat){
    cont <- t(object@drug)%*%object@ae
    coord <- which(cont!= 0, arr.ind=TRUE)
    coord <- coord[order(coord[,1]),] 
    dLab <- rownames(cont)[coord[,1]]
    dLab <- as.factor(dLab)
    aeLab <- colnames(cont)[coord[,2]] 
    aeLab <- as.factor(aeLab)    
    idx <- match(strat, c(colnames(object@cov)))
    
    if (sum(is.na(idx)>0)) stop("at least one of the stratification covariates does not exist in the object")
    varStrat <- factor(apply( object@cov[,idx, drop=F], 1, paste, collapse=" "))
    L <- nlevels(varStrat)
    print(table(varStrat))

    aeSplit <- vector("list", length=nlevels(varStrat))
    drugSplit <- vector("list", length=nlevels(varStrat))
    for (i in 1: nlevels(varStrat)){
      idx <- varStrat == levels(varStrat)[i]
      if(length(idx !=0)){
        drugSplit[[i]] <- object@drug[idx, ,drop=F]
        aeSplit[[i]] <- object@ae[idx, ,drop=F]
      }
    }
    n11 <- matrix(0, nrow=nrow(coord),  ncol=L)
    expN <- matrix(0, nrow=nrow(coord),  ncol=L)
    dMargin <- matrix(0, nrow=nrow(coord),  ncol=L)
    aeMargin <- matrix(0, nrow=nrow(coord),  ncol=L)
    N <- vector("numeric", length=L)
    for(l in 1:L)
    {
      cont <- t(drugSplit[[l]])%*%aeSplit[[l]]
      n11[,l] <- cont[coord] 
      n1._mat <- apply(cont,1,sum) 
      n.1_mat <- apply(cont,2,sum)
      dMargin[,l] <- n1._mat[coord[,1]] 
      aeMargin[,l] <- n.1_mat[coord[,2]] # ... et sa marge "colonne"
      N[l] <- sum(n11[,l]) 
      expN[,l] <- dMargin[,l]*aeMargin[,l]/N[l]
    }
    names(N) <- levels(varStrat)
    new("pvCont", dLab=dLab, aeLab=aeLab, n=n11, dMargin=dMargin,  aeMargin= aeMargin, expN=expN, N=N, strat=levels(varStrat))                
  }
)

setMethod(
  f = "pvCont",
  signature = c(object="data.frame", strat="missing"),
  definition = function(object, strat){
    if (ncol(object)!=3) stop("object must contain 3 columns : the drug labels, the ae labels and the corresponding adr counts")
    
    dLab <- factor(object[,1])
    aeLab <- factor(object[,2])
    n <- as.numeric(object[, 3])
    data <- data.frame(dLab, aeLab, n)
    
    cont <- xtabs(n~.,data=data, sparse = T)
    
    n1. <- as.numeric(cont %*% matrix(1, nrow=ncol(cont)))
    n.1 <- as.numeric(matrix(1, ncol=nrow(cont)) %*% cont ) 
    
    coord <- which(cont!=0, arr.ind=TRUE) 
    coord <-coord[order(coord[,1]),]
    n <- matrix(cont[coord], ncol=1) 
    N <- sum(n) # le nombre total de notifications
    dMargin <- matrix(n1.[coord[,1]], ncol=1)  # on affecte à chaque notification sa marge "ligne"...
    aeMargin <- matrix(n.1[coord[,2]], ncol=1) # ... et sa marge "colonne"
    expN<-dMargin*aeMargin/N
    pvCont<-new(Class="pvCont",dLab=dLab, aeLab=aeLab, n=n, dMargin= dMargin, aeMargin=aeMargin, expN=expN, N=N, strat=NULL)
    return(pvCont)
  }  
)

setMethod(
  f = "pvCont",
  signature = c(object="matrix", strat="missing"),
  definition = function(object, strat){
    object <- as.data.frame(object)
    pvCont(object)
  }
)
    