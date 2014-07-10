#' @title pvInd
#' @encoding UTF-8
#' @name pvInd
#' @description This function (which takes the same name that its corresponding formal class \code{\link{pvInd-class}}) should be used to convert a data.frame containing individual spontaneous reports (and optionnaly supplementary individual information) into an pvInd object.
#' @param adr a data.frame with three columns:
#' \itemize{
#' \item spontaneous report identifier
#' \item drug label
#' \item adverse event label
#' }
#' @param drug data.frame or matrix with two columns
#' \itemize{
#' \item spontaneous report identifier
#' \item drug label
#' }
#' @param ae data.frame or matrix with two columns
#' \itemize{
#' \item spontaneous report identifier
#' \item adverse event label
#' }
#' @param cov data.frame which contains individual covariates that can be used for stratified analysis. It is assumed that the first column of the data.frame corresponds to the spontaneous report identifier.
#' @param ... not of any use
#' @description \code{pvInd} is used to convert raw data  (\code{adr} and ) \code{cov}) into an pvInd object that can be used in the signal detection method functions.
#' @author Ismaïl Ahmed
#' @return a pvInd object
#' @keywords pvInd
#' @export
#' @docType methods
# @usage
# pvInd(drug, ae, adr, cov, ...)
# pvInd(drug, ae, ...)
# pvInd(drug, ae, cov, ...)
# pvInd(adr, cov, ...)
# pvInd(adr, ...)
#' @seealso \code{\link{pvInd-class}}
#' @aliases 
#' pvInd,data.frame,data.frame,missing,missing-method 
#' pvInd,data.frame,data.frame,missing,data.frame-method
#' pvInd,missing,missing,data.frame,missing-method
#' pvInd,missing,missing,data.frame,data.frame-method
#' 
#' #' @aliases pvCont,pvInd,character-method pvCont,pvInd,missing-method pvCont,data.frame,missing-method pvCont,matrix,missing-method
#' 
# @param drugMarginLim Minimum number of spontaneous reports for a drug to be included (default is 1)
# @param aeMarginLim Minimum number of spontaneous reports for an adverse event to be included (default is 1) 


# pvInd method definition ---------------------------------------------------- 

# generic definition ------------------------------------------------------
#' @export
setGeneric(
  name="pvInd",
  def=function(drug, ae, adr, cov, ...){standardGeneric("pvInd")}
)


# drug = "data.frame", ae = "data.frame", cov = "missing", adr = "missing" -----
setMethod(
  f="pvInd",
  signature = c(drug = "data.frame", ae = "data.frame", cov = "missing", adr = "missing"),
  definition = function(drug, ae){
    
    if(ncol(drug)!=2) { stop("drug must contain two columns")}
    if(ncol(ae)!=2) { stop("ae must contain two columns")}
    
    drug[,1] <- factor(drug[,1])
    drug[,2] <- factor(drug[,2])
    ae[,1] <- factor(ae[,1])
    ae[,2] <- factor(ae[,2])
    
    # drug Matrix
    colnD <- levels(drug[,2])
    rownD <- levels(drug[,1])
    idxR <- as.numeric(drug[,1])
    idxC <- as.numeric(drug[,2])
    sD <- sparseMatrix(idxR, idxC, x=1)
    idx <- which(sD>1, arr.ind = T)
    sD[idx] <- 1
    rownames(sD) <- rownD
    colnames(sD) <- colnD
    
    # ae Matrix
    colnAe <- levels(ae[,2])
    rownAe <- levels(ae[,1])
    idxR <- as.numeric(ae[,1])
    idxC <- as.numeric(ae[,2])
    sAe <- sparseMatrix(idxR,idxC,x=1)
    idx <- which(sAe>1, arr.ind = T)    
    sAe[sAe>1] <- 1
    rownames(sAe) <- rownAe
    colnames(sAe) <- colnAe
    
    idx <- match(rownD, rownAe)
    if(sum(is.na(idx))) warning("Some observations are missing or mismatch between the drug and ae datasets")
    sD <- sD[!is.na(idx), ]
    sAe <- sAe[idx, ]
    
    nSp <- nrow(sD)
    drugMargin <- as.numeric(Matrix(1, nrow=1, ncol=nSp) %*% sD)
    aeMargin <- as.numeric(Matrix(1, nrow=1, ncol=nSp) %*% sAe)
    
    pvInd<-new(Class="pvInd", drug=sD, ae=sAe, drugMargin=drugMargin, aeMargin=aeMargin)
    return(pvInd)
  }
)

# drug = "data.frame", ae = "data.frame", cov = "data.frame", adr = "missing"" ----
setMethod(
  f="pvInd",
  signature = c(drug = "data.frame", ae = "data.frame", cov = "data.frame", adr = "missing"),
  definition = function(drug, ae, cov){
    
    pvInd <- pvInd(drug, ae)
    if (ncol(cov)>0){
      idx <- match(row.names(pvInd@ae), as.character(cov[,1]))
      if(sum(is.na(idx))) warning("Some observations are missing or mismatch between adr and covariates datasets")
      covSort <- cov[idx,]
      row.names(covSort) <- covSort[,1]
      covSort <- covSort[,-1]      
    }else{
      covSort <- data.frame()
    }
    pvInd$cov <- covSort
    return(pvInd)
  }      
)


# adr = "data.frame", cov = "missing", drug = "missing", ae = "missing") ----
#' @export
setMethod(
  f="pvInd",
  signature = c(adr = "data.frame", cov = "missing", drug = "missing", ae = "missing"),
  definition = function(adr){
    
    if(nrow(adr)==0) { stop("adr must contain at least one row")}
    if(ncol(adr)!=3) { stop("adr must contain three columns")}
    
    adr[,1] <- factor(adr[,1])
    adr[,2] <- factor(adr[,2])
    adr[,3] <- factor(adr[,3])
    
    # drug Matrix
    colnD <- levels(adr[,2])
    rownD <- levels(adr[,1])
    idxR <- as.numeric(adr[,1])
    idxC <- as.numeric(adr[,2])
    sD <- sparseMatrix(idxR,idxC,x=1)
    sD[sD>1] <- 1
    rownames(sD) <- rownD
    colnames(sD) <- colnD
    
    # ae Matrix
    colnAe <- levels(adr[,3])
    idxC <- as.numeric(adr[,3])
    sAe <- sparseMatrix(idxR,idxC,x=1)
    sAe[sAe>1] <- 1
    rownames(sAe) <- rownD
    colnames(sAe) <- colnAe
    
    nSp <- nrow(sD)
    drugMargin <- as.numeric(Matrix(1, nrow=1, ncol=nSp) %*% sD)
    aeMargin <- as.numeric(Matrix(1, nrow=1, ncol=nSp) %*% sAe)
    
    #     if (ncol(adr)>3){
    #       sel <- !duplicated(adr[,1])
    #       cov <- adr[sel, 4:ncol(adr), drop=F]
    #       rownames(cov) <- as.character(adr[sel,1])
    #     }
    cov <-data.frame()
    pvInd<-new(Class="pvInd", drug=sD, ae=sAe, drugMargin=drugMargin, aeMargin=aeMargin, cov=cov)
    return(pvInd)
  }
)

# adr = "data.frame", cov="data.frame", drug = "missing", ae = "missing") ----
setMethod(
  f="pvInd",
  signature = c(adr = "data.frame", cov="data.frame", drug = "missing", ae = "missing"),
  definition = function(adr, cov){
    pvInd <- pvInd(adr)
    
    if (ncol(cov)>0){
      idx <- match(row.names(pvInd@ae), as.character(cov[,1]))
      if(sum(is.na(idx))) warning("Some observations are missing or mismatch between adr and covariates datasets")
      covSort <- cov[idx,]
      row.names(covSort) <- covSort[,1]
      covSort <- covSort[,-1]      
    }else{
      covSort <- data.frame()
    }
    pvInd$cov <- covSort
    return(pvInd)
  }
)


