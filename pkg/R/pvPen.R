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
#' @param object an object of class PvInd
#' @param covId a character vector indicating which covariate have to be used in the analyses. As current implementation of the lasso (glmnet) does not 
#' Such variables are first use in a logistic regression model and the residual is then used as the outcome 
#' in the penalized regression
#' @param aeId The label of the adverse event to be regressed. By default, the \code{pvPen} will regress all ae, one at a time and this can be very time consuming and RAM demanding. This parameter can be either the name of the ae(s), or the index of the ae column
#' @param detectCriter Can be either BIC or AIC. These criteria are used to select the final model.
#' @param posConst If TRUE (default), the regression coefficients are constrained to be non non negative, this to ensure that the final model only contains drugs "increasing" the risk of a given ae.
#' @param nDrugMax Maximum number of drugs to be included in the model. In addition to be computationnaly intensive, odels with to many drugs are likely to be very unstable.
#' @param parallel Whether parallel computations should be used to speed up the calculation. Be careful as it will be more RAM demanding. Only available for linux or Mac Os
#' @param nCores If parallel is true, allows to specify the number of cores to be used for parallel computing.
#' @param ... Further arguments to be passed to other functions (None for the moment)
#' @export
#' @usage
#' \method{pvPen}{pvInd}(object, aeId="all", covId=NULL, detectCriter=c("BIC", "AIC"), posConst=TRUE, nDrugMax=20, parallel=require(parallel), nCores=NULL, \dots)
#' 
# multiple regression.
# biclustering


# pvPen definition --------------------------------------------------------
pvPen <- function (object, ...) UseMethod("pvPen")

#' @export
# pvPen pvInd -------------------------------------------------------------
pvPen.pvInd <- function(object, aeId="all", covId=NULL, detectCriter=c("BIC", "AIC"), posConst=TRUE, nDrugMax=20, parallel=require(parallel), nCores=NULL, ...){
  
  if(parallel && !require(parallel)) stop("parallel package requested but not installed")
  if(parallel && is.null(nCores)) nCores <- parallel::detectCores()
  if (!inherits(object, "pvInd")) stop("object must be of class PvInd")
  
  detectCriter <- match.arg(detectCriter)
  lower.limit <- ifelse(posConst, 0, -Inf)
  pvPenCall = match.call(expand.dots = TRUE)
  #which = match(c("type.measure", "nfolds", "foldid", "grouped", 
#                  "keep"), names(glmnet.call), F)
  #print(pvPenCall)
  #print(names(pvPenCall))
  x <- object@drug
  if (aeId[1]=="all") { ## match.arg ?
    ae <- object@ae
  }else{
    ae <- object@ae[,aeId, drop=F]
  }
  nAe <- ncol(ae)
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
        colnames(temp) <- paste(covId[i], "_", levels(object@cov[[covId[i]]]), sep="")
        covGlmnet <- cbind(covGlmnet, temp[,-1])
      } else  {
        covGlmnet <- cbind(covGlmnet, object@cov[,covId[i], drop=F])
      }
    }    
  }
  
  resGlmnet <- vector("list", length=nAe)  
  if(!is.null(covId)){
    penalty.factor <- c(rep(0, ncol(covGlmnet)), rep(1, ncol(x)))
    lower.limits <- c(rep(-Inf, ncol(covGlmnet)), rep(lower.limit, ncol(x)))
  }else{
    lower.limits <- lower.limit
  }
  
  dev <- vector("list", length=nAe)
  bic <- vector("list", length=nAe)
  nParam <- vector("list", length=nAe)
  resGlm <- vector("list", length=nAe)
  jerrGlmnet <- vector("numeric", length=nAe)
  
  # Glmnet adjustment -------------------------------------------------------
  
  for(i in 1:nAe){
    cat("adverse event: ", colnames(ae)[i], "\n")
    y <- as.numeric(ae[,i])
    
    if(!is.null(covId)){
      resGlmnet[[i]] <- glmnet(cBind(covGlmnet,x), y, family="binomial", standardize=F, penalty.factor=penalty.factor, dfmax=nDrugMax, lower.limits=lower.limits, ...)
    }else{
      resGlmnet[[i]] <- glmnet(x, y, family="binomial", standardize=F, dfmax=nDrugMax, lower.limits=lower.limits, ...)
    }
    jerrGlmnet[i] <- resGlmnet[[i]]$jerr
    #print(resGlmnet[[i]]$jerr)
    #print(max(resGlmnet[[i]]$df))
    if ((jerrGlmnet[i] == 0) && (max(resGlmnet[[i]]$df>0)))  { ## check convergence of glmnet
      
      idxDf <- !duplicated(resGlmnet[[i]]$df)
      betaCoef <- resGlmnet[[i]]$beta[,idxDf]!=0
      betaCoef <- betaCoef[,apply(betaCoef,2,sum)>=1] #remove s0
      
      res <- vector("list", length=ncol(betaCoef))
      dev[[i]] <- vector("numeric", length=ncol(betaCoef))
      bic[[i]] <- vector("numeric", length=ncol(betaCoef))
      nParam[[i]] <- vector("numeric", length=ncol(betaCoef))
      
      if (parallel==FALSE){
        for (j in 1:ncol(betaCoef)){
          bSel <- betaCoef[,j]==1
          if(!is.null(covId)){
            xGlm <- as.matrix(cBind(covGlmnet,x)[,bSel])
            #glm.lower.limits <- c(rep(-Inf,ncol(covGlmnet)), rep(lower.limit, (ncol(xGlm)-ncol(covGlmnet))))
          }else{
            xGlm <- as.matrix(x[,bSel])
            #glm.lower.limits <- lower.limit
          }
          colnames(xGlm) <- row.names(betaCoef)[bSel]
          xGlm <- as.data.frame(xGlm)
          #print(dim(xGlm))
          #res[[j]] <- glmnet(xGlm, y, family="binomial", standardize=F, lower.limits=glm.lower.limits, lambda=0)
          #res[[j]] <- glmnet(xGlm, y, family="binomial", standardize=F, lower.limits=glm.lower.limits, lambda=0)
          
          res[[j]] <- glm(y~., data=xGlm, family="binomial")
          #print(res[[j]]$coefficients)
          if (posConst & (sum(res$coefficients[-1]<0)>0)) warning("negative coef in the regression step")
          dev[[i]][j] <- res[[j]]$deviance
          bic[[i]][j] <- dev[[i]][j] + (ncol(xGlm)+1)*log(nObs)
          nParam[[i]][j] <- ncol(xGlm)+1
        }
      }else{
        xGlm <- vector("list", length=ncol(betaCoef))
        for (j in 1:ncol(betaCoef)) {
          bSel <- betaCoef[,j]==1
          xGlm[[j]] <- x[,bSel, drop=F]
          colnames(xGlm[[j]]) <- row.names(betaCoef)[bSel]
        }

        if(!is.null(covId)){
          res <- mclapply(xGlm, .glmPar, y=y, cov=object@cov[, covId], mc.cores=nCores)
        }else{
          res <- mclapply(xGlm, .glmPar, y=y, cov=NULL, mc.cores=nCores)
        }

        for (j in 1:ncol(betaCoef)) { 
          dev[[i]][j] <- res[[j]]$deviance
          bic[[i]][j] <- dev[[i]][j] + (ncol(xGlm[[j]])+1)*log(nObs)
          nParam[[i]][j] <- ncol(xGlm[[j]])+1
        }
      }
      idxMin <- which.min(bic[[i]])
      if (idxMin==length(bic[[i]])) warning("The best BIC is obtained with about nDrugMax variables")
      resGlm[[i]] <- res[[idxMin]]
    }
  }

  
  resFinal <- vector("list")
  resFinal$bic <- bic
  resFinal$dev <- dev
  resFinal$jerrGlmnet <- jerrGlmnet
  resFinal$nParam <- nParam
  resFinal$glmnet <- resGlmnet  
  resFinal$bestGlm <- resGlm
  resFinal  
}

.glmPar <- function(x, y, cov=NULL, family="binomial"){
  x <- as.matrix(x)
  if (!is.null(cov)) x <- cbind(x, cov)
  x <- as.data.frame(x)
  res <- glm(y~., data=x, family =  family)
  res
  ## pas de warnings sur les poscontraints
}
# 
# 
# .LogisConst<-function(x,y,wt=rep(1,length(y)),intercept=TRUE,lower=rep(0,dim(as.data.frame(x))[2])) {
#   ## this function was largely adapted from Venables, W. N. and Ripley, B. D. (2002) Modern Applied Statistics with S. New York: Springer.
#   cat("function fmin ok")
#   fmin<-function(betas,X,y,w){
#     p<-plogis(X%*%betas)
#     -sum(2*w*ifelse(y,log(p),log(1-p)))
#   }
#   gmin<-function(betas,X,y,w){
#     eta<-X%*%betas
#     p<-plogis(eta)
#     gam<--2*(w*dlogis(eta)*ifelse(y,1/p,-1/(1-p)))
#     t(gam) %*%X
#   }
#   
#   init<-rep(0,dim(as.data.frame(x))[2])
#   if(is.null(dim(x))) dim(x)<-c(length(x),1)
#   dn<-colnames(x)
#   print(dn)
#   if(!length(dn)) dn<-paste("Var", 1:dim(as.data.frame(x))[2],sep="")
#   p<-dim(as.data.frame(x))[2] + intercept
#   if(intercept) {x<-cbind(1,x);cat(p);dn<-c("(Intercept)",dn);bas<-c(-Inf,lower)}
#   if(is.factor(y)) y<-unclass(y)!=1
#   #fit<-nlminb(rep(0,dim(as.data.frame(x))[2]),fmin,gmin,lower=bas,X=x,y=y,w=wt)
#   fit<-optim(par=rep(0,ncol(as.data.frame(x))),fn=fmin,gr=gmin,lower=bas,method="L-BFGS-B",X=x,y=y,w=wt)
#   names(fit$par)<-dn
#   cat("\nCoefficients:\n");print(fit$par)
#   #cat("\nResidual Deviance:",format(fit$objective),"\n")
#   cat("\nConvergence message:", fit$convergence,"\n")
#   invisible(fit)
# }
# 
# .paraStabSel <- function(idSel, x, y, alpha, family, lambda) {
#   resGlmnet <- glmnet(x[idSel,], y[idSel], alpha = alpha, family = family, lambda = lambda, lower.limits=0)
#   resGlmnet$beta>0
# }
# 
# .paraStabSel2 <- function(idSel, x, y, alpha, family, lambda) {
#   resGlmnet <- glmnet(x[idSel,], y[idSel], alpha = alpha, family = family, lambda = lambda, lower.limits=0)
#   resGlmnet$beta>0
# }
# 
# # if (detectCriter=="stabSel"){
# #   idSel <- as.data.frame(matrix(sample(1:nObs, nBoot*sizeBoot, replace=T), ncol=nBoot, nrow=sizeBoot))
# #   print(head(idSel[,1:10]))
# #   if (parallel){
# #     resBeta <- mclapply(idSel, .paraStabSel2, x, y, 1, family ="gaussian")
# #   }else{
# #     resBeta <- lapply(idSel, .paraStabSel2, x, y, 1, family ="gaussian")
# #   }
# #   minCol <- min(unlist(lapply(resBeta, ncol)))
# #   print(maxCol) 
# #   betaSel <- Reduce('+', resBeta)
# # }
