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
#' @param cov a character vector indicating which covariate have to be used in the analyses. As current implementation of the lasso (glmnet) does not 
#' Such variables are first use in a logistic regression model and the residual is then used as the outcome 
#' in the penalized regression
#' @param aeLab The label of the adverse event to be regressed. By default, the \code{pvPen} will regress all ae, one at a time and this can be very time consuming and RAM demanding. 
#' @param detectCriter Can be either BIC or AIC. These criteria are used to select the final model.
#' @param posConst If TRUE (default), the regression coefficients are constrained to be non non negative, this to ensure that the final model only contains drugs "increasing" the risk of a given ae.
#' @param nDrugMax Maximum number of drugs to be included in the model. In addition to be computationnaly intensive, odels with to many drugs are likely to be very unstable.
#' @param parallel Whether parallel should be used in order to speed up the calculation. Be careful as it will be more RAM demanding. Only available for linux or Mac Os
#' @param nCore If parallel is true, allows to specify the number of cores to be used for parallel computing.
#' @param ... Further arguments to be passed to other functions (None for the moment)
#' @export
#' @usage
#' \method{pvPen}{pvInd}(object, aeLab="all", cov=NULL, detectCriter=c("BIC", "AIC"), posConst=TRUE, nDrugMax=20, parallel=FALSE, nCore=NULL, \dots)
#' 
# multiple regression.
# biclustering


# pvPen definition --------------------------------------------------------
pvPen <- function (object, ...) UseMethod("pvPen")

#' @export
# pvPen pvInd -------------------------------------------------------------
pvPen.pvInd <- function(object, aeLab="all", cov=NULL, detectCriter=c("BIC", "AIC"), posConst=TRUE, nDrugMax=20, parallel=FALSE, nCore=NULL, ...){
  if (!inherits(object, "pvInd")) stop("object must be of class PvInd")
  #if (!is.null(cov) & !is.null(strat)) stop("strat and cov cannot be both null")
  if(!is.null(cov)){
    if (!is.character(cov)){stop("object must be a vector of character")}
    covMat <- c()  
    for (i in 1:length(cov)){  
      if (is.factor(object@cov[[cov[i]]])) {
        temp <- model.matrix(~object@cov[[cov[i]]]-1)
        colnames(temp) <- paste(cov[i], "_", levels(object@cov[[cov[i]]]), sep="")
        covMat <- cbind(covMat, temp[,-1])
      } else  {
        covMat <- cbind(covMat, object@cov[,cov[i], drop=F])
      }
    }
    
  }

  #print(head(covMat))
  
  detectCriter <- match.arg(detectCriter)
  lower.limit <- ifelse(posConst, 0, -Inf)

  x <- object@drug
  if (aeLab[1]=="all") {
    ae <- object@ae
  }else{
    ae <- object@ae[,aeLab, drop=F]
  }
  nAe <- ncol(ae)
  nObs <- nrow(ae)
  
  resGlmnet <- vector("list", length=nAe)  
  penalty <- c(rep(0, ncol(covMat)), rep(1, ncol(x)))
  lowerLim <- c(rep(-Inf, ncol(covMat)), rep(0, ncol(x)))
  
  dev <- vector("list", length=nAe)
  bic <- vector("list", length=nAe)
  nParam <- vector("list", length=nAe)
  for(i in 1:nAe){
    cat(i, ",", colnames(object@ae)[i], "\n")
    y <- as.numeric(ae[,i])
    resLogCov <- glm(y~., data= object@cov[,cov], family="binomial")
    print(resLogCov$aic)
    print(logLik(resLogCov))
    print(deviance(resLogCov))
    resGlmnet[[i]] <- glmnet(cBind(covMat,x), y, family="binomial", standardize=F, penalty.factor=penalty, dfmax=nDrugMax, lower.limits=lowerLim, nlambda=200)
    print(resGlmnet[[i]])
    print(resGlmnet[[i]]$nulldev)
    print((1-resGlmnet[[i]]$dev.ratio)*resGlmnet[[i]]$nulldev)
    idxDf <- !duplicated(resGlmnet[[i]]$df)
    betaCoef <- resGlmnet[[i]]$beta[,idxDf]!=0
    print(dim(betaCoef))
    print(apply(betaCoef,2,sum))
    
    res <- vector("list", length=ncol(betaCoef))
    dev[[i]] <- vector("numeric", length=ncol(betaCoef))
    bic[[i]] <- vector("numeric", length=ncol(betaCoef))
    nParam[[i]] <- vector("numeric", length=ncol(betaCoef))

    for (j in 1:ncol(betaCoef)){
      
      xGlm <- as.matrix(cBind(covMat,x)[,betaCoef[,j]==1])
      #print(dim(xGlm))
      lowerLogisConst <- c(rep(-Inf,ncol(covMat)), rep(0, (ncol(xGlm)-ncol(covMat))))
      #res[[j]] <- .LogisConst(xGlm,y,wt=rep(1,length(y)),intercept=TRUE,lower=lowerLogisConst)
      res[[j]] <- glmnet(xGlm,y,family="binomial", standardize=F, penalty.factor=0, lower.limits=lowerLim, lambda=0)
      #print(res[[j]])
      #print(res[[j]]$nulldev)
      dev[[i]][j] <- (1-res[[j]]$dev.ratio)*res[[j]]$nulldev
      bic[[i]][j] <- dev[[i]][j] + (ncol(xGlm)+1)*log(nObs)
      nParam[[i]][j] <- ncol(xGlm)+1
    }
    #plot(apply(betaCoef, 2, sum), bic)
  }
  resFinal <- vector("list")
  resFinal$bic <- bic
  resFinal$dev <- dev
  resFinal$nParam <- nParam
  resFinal$glmnet <- resGlmnet  
  resFinal
}


.LogisConst<-function(x,y,wt=rep(1,length(y)),intercept=TRUE,lower=rep(0,dim(as.data.frame(x))[2])) {
  ## this function was largely adapted from Venables, W. N. and Ripley, B. D. (2002) Modern Applied Statistics with S. New York: Springer.
  cat("function fmin ok")
  fmin<-function(betas,X,y,w){
    p<-plogis(X%*%betas)
    -sum(2*w*ifelse(y,log(p),log(1-p)))
  }
  gmin<-function(betas,X,y,w){
    eta<-X%*%betas
    p<-plogis(eta)
    gam<--2*(w*dlogis(eta)*ifelse(y,1/p,-1/(1-p)))
    t(gam) %*%X
  }
  
  init<-rep(0,dim(as.data.frame(x))[2])
  if(is.null(dim(x))) dim(x)<-c(length(x),1)
  dn<-colnames(x)
  print(dn)
  if(!length(dn)) dn<-paste("Var", 1:dim(as.data.frame(x))[2],sep="")
  p<-dim(as.data.frame(x))[2] + intercept
  if(intercept) {x<-cbind(1,x);cat(p);dn<-c("(Intercept)",dn);bas<-c(-Inf,lower)}
  if(is.factor(y)) y<-unclass(y)!=1
  #fit<-nlminb(rep(0,dim(as.data.frame(x))[2]),fmin,gmin,lower=bas,X=x,y=y,w=wt)
  fit<-optim(par=rep(0,ncol(as.data.frame(x))),fn=fmin,gr=gmin,lower=bas,method="L-BFGS-B",X=x,y=y,w=wt)
  names(fit$par)<-dn
  cat("\nCoefficients:\n");print(fit$par)
  #cat("\nResidual Deviance:",format(fit$objective),"\n")
  cat("\nConvergence message:", fit$convergence,"\n")
  invisible(fit)
}

.paraStabSel <- function(idSel, x, y, alpha, family, lambda) {
  resGlmnet <- glmnet(x[idSel,], y[idSel], alpha = alpha, family = family, lambda = lambda, lower.limits=0)
  resGlmnet$beta>0
}

.paraStabSel2 <- function(idSel, x, y, alpha, family, lambda) {
  resGlmnet <- glmnet(x[idSel,], y[idSel], alpha = alpha, family = family, lambda = lambda, lower.limits=0)
  resGlmnet$beta>0
}

# if (detectCriter=="stabSel"){
#   idSel <- as.data.frame(matrix(sample(1:nObs, nBoot*sizeBoot, replace=T), ncol=nBoot, nrow=sizeBoot))
#   print(head(idSel[,1:10]))
#   if (parallel){
#     resBeta <- mclapply(idSel, .paraStabSel2, x, y, 1, family ="gaussian")
#   }else{
#     resBeta <- lapply(idSel, .paraStabSel2, x, y, 1, family ="gaussian")
#   }
#   minCol <- min(unlist(lapply(resBeta, ncol)))
#   print(maxCol) 
#   betaSel <- Reduce('+', resBeta)
# }
