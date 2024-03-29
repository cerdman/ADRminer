#' @encoding UTF-8
#' @name pvCont
#' @title pvCont constructor
#' @param object a pvInd object or a data.frame. The latter should have 3 columns, i.e. the drug labels, the ae labels and the corresponding number of spontaneous reports. 
#' @param strat Character vector containing the name of the covariates to be used for a stratified analysis.  \code{strat} is only meaningful when used along with a pvInd object.
#' @description pvCont is to be used to convert raw aggregated spontaneaous reporting data into a pvCont object (\code{pvCont-class}). This constructor can also be used to convert pvInd object into pvCont object. 
#' @author Ismaïl Ahmed
#' @return an object of class \code{\link{pvCont}}.
#' @include pvCont-class.R
#' @include pvInd-class.R
#' @exportMethod pvCont
#' @keywords pvCont
#' @docType methods
#' @aliases pvCont,pvInd,character-method pvCont,pvInd,missing-method pvCont,data.frame,missing-method pvCont,matrix,missing-method
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
    drugMargin <- matrix(n1.[coord[,1]], ncol=1)
    aeMargin <- matrix(n.1[coord[,2]], ncol=1)
    N <- sum(n11) 
    expN <- drugMargin*aeMargin/N
    pvCont<-new(Class="pvCont", drugLab=dLab, aeLab=aeLab, n=n11, drugMargin=drugMargin, aeMargin=aeMargin, expN=expN, N=N, strat=character()) #, coord=coord)
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
    for (i in 1: length(strat)){
      if (!is.factor(object$cov[[strat[i]]])) {
        warning(strat[i] ," is not defined as factor - Attempt to convert this covariate into a factor", immediate. = TRUE)
        object$cov[[strat[i]]] <- factor(object$cov[[strat[i]]])
        if (nlevels(object$cov[[strat[i]]]) >= 10) stop("Too many levels for the factor")  
      }
    }
    varStrat <- factor(apply( object@cov[,idx, drop=F], 1, paste, collapse=" "))
    L <- nlevels(varStrat)
    if (L >= 10) warning("The number of strates is ", L, ". Be careful as it may create numerical and statistical instabilities", immediate. = TRUE)
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
    drugMargin <- matrix(0, nrow=nrow(coord),  ncol=L)
    aeMargin <- matrix(0, nrow=nrow(coord),  ncol=L)
    N <- vector("numeric", length=L)
    for(l in 1:L)
    {
      cont <- t(drugSplit[[l]])%*%aeSplit[[l]]
      n11[,l] <- cont[coord] 
      n1._mat <- apply(cont,1,sum) 
      n.1_mat <- apply(cont,2,sum)
      drugMargin[,l] <- n1._mat[coord[,1]] 
      aeMargin[,l] <- n.1_mat[coord[,2]] 
      N[l] <- sum(n11[,l]) 
      expN[,l] <- drugMargin[,l]*aeMargin[,l]/N[l]
    }
    names(N) <- levels(varStrat)
    new("pvCont", drugLab=dLab, aeLab=aeLab, n=n11, drugMargin=drugMargin,  aeMargin= aeMargin, expN=expN, N=N, strat=levels(varStrat))                
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
    drugMargin <- matrix(n1.[coord[,1]], ncol=1)  # on affecte à chaque notification sa marge "ligne"...
    aeMargin <- matrix(n.1[coord[,2]], ncol=1) # ... et sa marge "colonne"
    expN<-drugMargin*aeMargin/N
    pvCont<-new(Class="pvCont",drugLab=dLab, aeLab=aeLab, n=n, drugMargin= drugMargin, aeMargin=aeMargin, expN=expN, N=N, strat=character())
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


