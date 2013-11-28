#' @encoding UTF-8
#' @title Cross-validation for glmnet and pvInd object
#' @description This function is mainly a wrapper of the cv.glmnet function (glmnet package) for pvInd object
#' @param object an object of class pvInd
#' @param aeId The label of the adverse event to be regressed. By default, the \code{cv.glmnet} will regress all ae, one at a time and this can be very time consuming and RAM demanding. This parameter can be either the name of the ae(s), or the index of the ae column
#' @param covId a character vector indicating which covariate have to be used in the analyses. 
#' @param posConst If TRUE (default), the regression coefficients are constrained to be non non negative, this to ensure that the final model only contains drugs "increasing" the risk of a given ae.
#' @param nDrugMax Maximum number of drugs to be included in the model. In addition to be computationnaly intensive, odels with to many drugs are likely to be very unstable.
#' @param parallel Whether parallel computations should be used to speed up the calculation. Be careful as it will be more RAM demanding. Only available for linux or Mac Os
#' @param ... Further arguments to be passed to the cv.glmnet function (not in use for now)
#' @export

cv.glmnet.pvInd <- function(object, aeId = "all", covId = NULL, posConst = TRUE, nDrugMax = 20, parallel = FALSE, ...){
  if(!inherits(object, "pvInd")) stop("x must be a pvInd object.")
  lower.limit <- ifelse(posConst, 0, -Inf)
  pvPenCall = match.call(expand.dots = TRUE)
  
  x <- object@drug
  if (aeId[1] == "all") { ## match.arg ?
    ae <- object@ae
  }else{
    ae <- object@ae[,aeId, drop = F]
  }
  nAe <- ncol(ae)
  nD <- ncol(x)
  nObs <- nrow(ae)
  
  # Transform vector of covariates into a dummy matrix for glmnet -----------
  if(!is.null(covId)){
    if (!is.character(covId)) stop("covId must be a vector of character")
    idxCov <- match(covId, colnames(object@cov))
    if(sum(is.na(idxCov))) stop("Covariates listed in covId must match these in object (see names(object))")
    covGlmnet <- c()  
    for (i in 1:length(cov)){  
      if (is.factor(object@cov[[covId[i]]])) {
        temp <- model.matrix(~object@covId[[cov[i]]]-1)
        colnames(temp) <- paste(covId[i], "_", levels(object@cov[[covId[i]]]), sep = "")
        covGlmnet <- cbind(covGlmnet, temp[,-1])
      } else  {
        covGlmnet <- cbind(covGlmnet, object@cov[,covId[i], drop = F])
      }
    }    
  }
  
  
  if(!is.null(covId)){
    penalty.factor <- c(rep(0, ncol(covGlmnet)), rep(1, ncol(x)))
    lower.limits <- c(rep(-Inf, ncol(covGlmnet)), rep(lower.limit, ncol(x)))
  }else{
    lower.limits <- lower.limit
  }
  
  
  #jerrGlmnet <- vector("numeric", length = nAe)
  resGlmnet <- vector("list", length = nAe) 
  
  for(i in 1:nAe){
    cat("adverse event: ", colnames(ae)[i], "\n")
    y <- as.numeric(ae[,i])
    
    if(!is.null(covId)){
      resGlmnet[[i]] <- cv.glmnet(cBind(covGlmnet,x), y, family="binomial", standardize = F, penalty.factor = penalty.factor, dfmax = nDrugMax, lower.limits = lower.limits, parallel=parallel, ...)
    }else{
      resGlmnet[[i]] <- cv.glmnet(x, y, family = "binomial", standardize = F, dfmax = nDrugMax, lower.limits = lower.limits, parallel=parallel, ...)
    }
    #jerrGlmnet[i] <- resGlmnet[[i]]$jerr    
  }  
  if (nAe == 1) {
    res <- resGlmnet[[1]]
  }else{
    res <- resGlmnet
    names(resGlmnet) <- colnames(ae)
  }
  return(res)
}