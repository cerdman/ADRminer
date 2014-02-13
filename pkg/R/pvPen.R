#' @encoding UTF-8
#' @title Penalized regression detection strategy
#' @name pvPen
#' @aliases pvPen.pvInd
#' @description pvPen implements a detection strategy based on the use of the Lasso regression along with BIC
#'  type detection criteria. The methodology follows two main steps: 
#'  \itemize{
#'  \item Lasso regression (with a positive constraint on the regression coefficients) for a grid of constraint 
#'  parameters
#'  \item selection on the drug associated with the ae of interest according to the BIC. 
#'  }
#'  This function heavily relies on the highly efficient \code{glmnet} package. \code{pvPen} is computationnaly highly demanding. It is also strongly advised to use a linux or a mac computer in order to use several cores. Note also that this regression approach should be used with ae and drugs both having a reasonable number of reports. 
#' @param object an object of class pvInd
#' @param aeId The label of the adverse event to be regressed. By default, the \code{pvPen} will regress all ae, one at a time and this can be very time consuming and RAM demanding. This parameter can be either the name of the ae(s), or the index of the ae column
#' @param covId a character vector indicating which covariate have to be used in the analyses.
#' @param criter Can be either BIC or AIC. These criteria are used to select the final model.
#' @param posConst If TRUE (default), the regression coefficients are constrained to be non non negative, this to ensure that the final model only contains drugs "increasing" the risk of a given ae.
#' @param nDrugMax Maximum number of drugs to be included in the model. In addition to be computationnaly intensive, odels with to many drugs are likely to be very unstable.
#' @param parallel Whether parallel computations should be used to speed up the calculation. Be careful as it will be more RAM demanding. Only available for linux or Mac Os
#' @param nCores If parallel is true, allows to specify the number of cores to be used for parallel computing.
#' @param ... Further arguments to be passed  (None for the moment)
#' @export
#' @usage
#' \method{pvPen}{pvInd}(object, aeId = "all", covId = NULL, criter = c("BIC", "AIC", "eBIC"), posConst = TRUE, nDrugMax = 20, parallel = require(parallel), nCores = NULL, \dots)



# pvPen definition --------------------------------------------------------
pvPen <- function (object, ...) UseMethod("pvPen")

#' @export
# pvPen pvInd -------------------------------------------------------------
pvPen.pvInd <- function(object, aeId = "all", covId = NULL, criter = c("BIC", "AIC", "eBIC"), posConst = TRUE, nDrugMax = 20, parallel = require(parallel), nCores = NULL, ...){
  
  if(parallel && !require(parallel)) stop("parallel package requested but not installed")
  if(parallel && is.null(nCores)) nCores <- parallel::detectCores()
  if (!inherits(object, "pvInd")) stop("object must be of class PvInd")
  
  criter <- match.arg(criter)
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
  
  dev <- vector("list", length = nAe)
  bic <- vector("list", length = nAe)
  aic <- vector("list", length = nAe)
  ebic <- vector("list", length = nAe)
  nParam <- vector("list", length = nAe)
  resGlmBIC <- vector("list", length = nAe)
  resGlmAIC <- vector("list", length = nAe)
  resGlmEBIC <- vector("list", length = nAe)
  jerrGlmnet <- vector("numeric", length = nAe)
  resGlmnet <- vector("list", length = nAe)  
  
  # Glmnet adjustment -------------------------------------------------------
  
  for(i in 1:nAe){
    cat("adverse event: ", colnames(ae)[i], "\n")
    y <- as.numeric(ae[,i])
    
    if(!is.null(covId)){
      resGlmnet[[i]] <- glmnet(cBind(covGlmnet,x), y, family="binomial", standardize = F, penalty.factor = penalty.factor, dfmax = nDrugMax, lower.limits = lower.limits, ...)
    }else{
      resGlmnet[[i]] <- glmnet(x, y, family = "binomial", standardize = F, dfmax = nDrugMax, lower.limits = lower.limits, ...)
    }
    jerrGlmnet[i] <- resGlmnet[[i]]$jerr
    
    if ((jerrGlmnet[i] == 0) && (max(resGlmnet[[i]]$df>0)))  { ## check convergence of glmnet
      
      idxDf <- !duplicated(resGlmnet[[i]]$df)
      betaCoef <- resGlmnet[[i]]$beta[,idxDf] != 0
      betaCoef <- betaCoef[,apply(betaCoef,2,sum) >= 1] #remove s0
      
      nGlm <- ncol(betaCoef)
      res <- vector("list", length = nGlm)
      dev[[i]] <- vector("numeric", length = nGlm)
      bic[[i]] <- vector("numeric", length = nGlm)
      aic[[i]] <- vector("numeric", length = nGlm)
      ebic[[i]] <- vector("numeric", length = nGlm)
      nParam[[i]] <- vector("numeric", length = nGlm)
      
      if (parallel == FALSE){
        for (j in 1:nGlm){
          betaSel <- betaCoef[,j] == 1
          if(!is.null(covId)){
            xGlm <- as.matrix(cBind(covGlmnet,x)[,betaSel])
          }else{
            xGlm <- as.matrix(x[,betaSel])
          }
          colnames(xGlm) <- row.names(betaCoef)[betaSel]
          xGlm <- as.data.frame(xGlm)          
          res[[j]] <- glm(y~., data = xGlm, family = "binomial")
        } 
      }else{
        xGlm <- vector("list", length = nGlm)
        for (j in 1:nGlm) {
          betaSel <- betaCoef[,j] == 1
          xGlm[[j]] <- x[,betaSel, drop = F]
          colnames(xGlm[[j]]) <- row.names(betaCoef)[betaSel]
        }
        
        if(!is.null(covId)){
          res <- mclapply(xGlm, .glmPar, y = y, cov = object@cov[, covId], mc.cores = nCores)
        }else{
          res <- mclapply(xGlm, .glmPar, y = y, cov = NULL, mc.cores = nCores)
        }
      } ## end if PARALLEL
      for (j in 1:nGlm) { ## Best Glm model
        dev[[i]][j] <- res[[j]]$deviance
        #bic[[i]][j] <- dev[[i]][j] + (ncol(xGlm[[j]])+1)*log(nObs)
        bic[[i]][j] <- BIC(res[[j]])
        aic[[i]][j] <- AIC(res[[j]])      
        nParam[[i]][j] <- ncol(xGlm[[j]])+1
        ebic[[i]][j] <- bic[[i]][j] + 2 * lchoose(nD,nParam[[i]][j])
      }
#       idxMin <- switch(criter, 
#                        BIC =  which.min(bic[[i]]),
#                        AIC =  which.min(aic[[i]]),
#                        eBIC = which.min(ebic[[i]])
#       )
#       
#       if (idxMin == nGlm) warning(paste("The best", criter, "is obtained with about nDrugMax variables."))
#       resGlm[[i]] <- res[[idxMin]]
      idxMinBIC <-  which.min(bic[[i]])
      idxMinAIC <-  which.min(aic[[i]])
      idxMinEBIC <- which.min(ebic[[i]])
      if (idxMinBIC == nGlm) warning("The best BIC is obtained with about nDrugMax variables.")
      if (idxMinAIC == nGlm) warning("The best AIC is obtained with about nDrugMax variables.")
      resGlmBIC[[i]] <- res[[idxMinBIC]]
      resGlmEBIC[[i]] <- res[[idxMinEBIC]]
      resGlmAIC[[i]] <- res[[idxMinAIC]]
      
      
    } ## en if on check glmnet convergence
  } ## end loop on AE
  
  resFinal <- vector("list", length = nAe)
  
  for (i in 1: nAe){
    resFinal[[i]] <- vector("list")
    resFinal[[i]]$BIC <- bic[[i]]
    resFinal[[i]]$AIC <- aic[[i]]
    resFinal[[i]]$eBIC <- ebic[[i]]
    resFinal[[i]]$jerrGlmnet <- jerrGlmnet[i]
    resFinal[[i]]$nParam <- nParam[[i]]
    #resFinal[[i]]$glmnet <- resGlmnet[[i]]  
    resFinal[[i]]$bestGlmBIC <- coefficients(summary(resGlmBIC[[i]]))
    resFinal[[i]]$bestGlmEBIC <- coefficients(summary(resGlmEBIC[[i]]))
    resFinal[[i]]$bestGlmAIC <- coefficients(summary(resGlmAIC[[i]]))
  }
  if (nAe == 1) {
    resFinal <- unlist(resFinal, recursive = F)
  }else{
    names(resFinal) <- colnames(ae)
  }
  
  resFinal  
}

.glmPar <- function(x, y, cov = NULL, family = "binomial"){
  x <- as.matrix(x)
  if (!is.null(cov)) x <- cbind(x, cov)
  x <- as.data.frame(x)
  res <- glm(y~., data = x, family =  family)
  res
  ## pas de warnings sur les poscontraints
}
