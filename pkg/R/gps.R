#' @encoding UTF-8
#' @title Gamma Poisson Shrinkage
#' @param object An object of class pvInd or pvCont
#' @param rr0 The tested risk. By default, rr0=1.
#' @param assocMeasure Statistic used to order the drug-events pairs: 
#' \itemize{
#' \item postH0: Posterior probability of the null hypothesis
#' \item lb05: 5\% quantile of the posterior distribution of lambda
#' \item postE: Posterior expectation of log(lambda,2)
#' }
#' @param detectCriter Decision rule for the signal generation based on:
#' FDR: (Default value) 
#' nSig: number of signals 
#' assocMeasure: ranking statistic. See \code{assocMeasure}
#' @param criterThres Threshold used for \code{detectCriter}. Ex 0.05 for FDR.
#' @param nMin Minimum number of spontaneous reports for a drug-event pair to be potentially considered as a signal. By default, \code{nMin}=1.
#' @param truncThres Used for the hyper parameter calculations. The marginal likelihood involves mixtures of two (possibly truncated) negative binomial. 
#' Default value is set to 0 as the information regarding drug-event pairs not associated with any adr is not used. By setting \code{truncThres} to NULL, the whole contingency table 
#' is rebuild and all data (thus including the 0s) are used to estimate the hyperparameters.
#' in order to account for that drug-event pairs associated with no spontaneous reports are usually not represented in the data. In order not to use any truncation, set \code{truncThres} to NULL
#' @param hyperparamInit Initialization vector for the prior mixture parameters (alpha1, beta1, alpha2, beta2, w). By default, hyperparamInit = c(alpha1 = 0.2, beta1 = 0.06, alpha2 = 1.4, beta2 = 1.8, w = 0.1), i.e. the prior parameters provided in DuMouchel (1999).
#' @param hyperparam Hyper parameter values for the 2 gamma mixture model. By default, hyperparam = NULL which means that the hyperparameters are calculated by maximising the marginal likelihood. This parameter is meant to avoid this maximization step
#' @param strat Character vector containing the name of the covariates to be used for a stratified analysis.
#' @param allRes Logical value indicating whether all drug-event combination results must be provided. 
#' @param ... Further arguments to be passed to other functions (None for the moment)
#' @description Gamma Poisson Shrinkage model proposed by DuMouchel (1999) extended to the multiple comparison framework.
#' @return allSig Data.frame summarizing the results for all drug-event combinations with at least \code{nMin} spontaneous reports ordered by \code{assocMeasure}. This output is provided if \code{allRes=TRUE}. Operating characteristics are estimated according to Ahmed et al (Stat Med 2009). 
#' @return sig Same as \code{allSig} but restricted to the list of generated signals.
#' @return nSig Number of generated signals.
#' @return call specified argumetnts
#' @author Youness Ergaibi & Ismaïl Ahmed
#' @keywords gps
#' @export
#' @docType methods
#' @aliases gps.pvInd
#' @aliases gps.pvCont
#' @rdname gps
#' @usage
#' \method{gps}{pvCont}(object, rr0=1, assocMeasure=c("postH0","lb05","postE"), detectCriter=c("FDR","nSig","assocMeasure"), criterThres = 0.05, nMin=1, truncThres = 0, hyperparamInit = c(alpha1= 0.2, beta1= 0.06, alpha2=1.4, beta2=1.8, w=0.1), hyperparam = NULL, allRes=F, \dots)
#' 
#' \method{gps}{pvInd}(object, rr0=1, assocMeasure=c("postH0","lb05","postE"), detectCriter=c("FDR","nSig","assocMeasure"), criterThres = 0.05, strat=NULL, nMin=1, truncThres = NULL, hyperparamInit = c(alpha1= 0.2, beta1= 0.06, alpha2=1.4, beta2=1.8, w=0.1), hyperparam = NULL, allRes=F, \dots)

# gps definition ----------------------------------------------------------
gps <- function (object, ...) UseMethod("gps")

# gps pvCont --------------------------------------------------------------
#' @export
gps.pvCont <- function(object, rr0=1, assocMeasure=c("postH0","lb05","postE"), detectCriter=c("FDR","nSig","assocMeasure"), criterThres = 0.05, nMin=1, truncThres = 0, hyperparamInit = c(alpha1= 0.2, beta1= 0.06, alpha2=1.4, beta2=1.8, w=0.1), hyperparam = NULL, allRes=F, ...){
  
  if(!inherits(object, "pvCont")) stop("x must be a pvCont object.")
  assocMeasure <- match.arg(assocMeasure) # keep only the first argument
  detectCriter <- match.arg(detectCriter)
  n <- object@n
  drugMargin <- object@drugMargin
  aeMargin <- object@aeMargin
  expN <- object@expN
  dLab <- object@drugLab
  aeLab <- object@aeLab
  N <- object@N
  # only keep the pairs associated with at least nMin sr
  if(nMin > 1) {
    expN <- expN[n >= nMin]
    drugMargin <- drugMargin[n >= nMin]
    aeMargin <- aeMargin[n >= nMin]
    n <-  n[n >= nMin]
  }
  nPair <- length(n)
  postH0 <- vector(length=nPair)
  if (is.null(hyperparam)){         
    likTrunc  <- function(p, n, expN, trunc){
      a1 <- p[1]
      b1 <- p[2]/(p[2]+expN)
      a2 <- p[3]
      b2 <- p[4]/(p[4]+expN)
      w <- p[5]
      sum(
        -log(
          w * (dnbinom(n, a1, b1) / pnbinom(trunc, a1, b1, lower.tail=F)) + (1-w) * (dnbinom(n, a2, b2)/ pnbinom(trunc, a2, b2, lower.tail=F))
        )
      )
    }
    lik  <- function(p, n, expN){
      a1 <- p[1]
      b1 <- p[2]/(p[2]+expN)
      a2 <- p[3]
      b2 <- p[4]/(p[4]+expN)
      w <- p[5]
      sum(
        -log(
          w * dnbinom(n, a1, b1) + (1-w) * dnbinom(n, a2, b2)
        )
      )
    }
    
    if (is.null(truncThres)){
      cont <- xtabs(n ~ dLab + aeLab, sparse=T)
      n1. <- apply(cont, 1, sum)
      n.1 <- apply(cont, 2, sum)
      n1.c <- rep(n1., times = length(n.1))
      n.1c <- rep(n.1, each = length(n1.))
      expNc <- n1.c * n.1c/N
      nc <- as.vector(cont)
      
      pOut  <- suppressWarnings(nlm(lik, p=hyperparamInit, n=nc, expN=expNc, iterlim=1000))
    }else{
      pOut  <- suppressWarnings(nlm(likTrunc, p=hyperparamInit, n=n[n>truncThres], expN=expN[n>truncThres], trunc=truncThres, iterlim=1000))
    }      
    hp <- pOut$estimate
    convergence <- pOut$code
  }
  
  Q <- hp[5]*dnbinom(n,size=hp[1],prob=hp[2]/(hp[2]+expN)) / 
    (hp[5]*dnbinom(n,size=hp[1],prob=hp[2]/(hp[2]+expN)) + (1-hp[5])* 
       dnbinom(n,size=hp[3],prob=hp[4]/(hp[4]+expN)))
  
  postH0 <- Q*pgamma(rr0,hp[1]+n,hp[2]+expN) +(1-Q)*pgamma(rr0,hp[3]+n,hp[4]+expN)     
  # Posterior Expectation of Lambda
  postE <- log(2)^(-1)*(Q*(digamma(hp[1]+n)-log(hp[2]+expN))+(1-Q)*(digamma(hp[3]+n)-log(hp[4]+expN)))
  
  # lb05
  lb05 <- .quantGps(0.05, Q, hp[1]+n, hp[2]+ expN, hp[3]+n, hp[4]+expN)   
  
  # Assignment based on the method (pp/postE/lb05)
  switch(assocMeasure, 
         postH0 = {assocM <- postH0; decrease <- F},
         lb05 = {assocM <- lb05; decrease <- T},
         postE = {assocM <- postE; decrease <- T}
  )
  
  FDR <- (cumsum(postH0[order(assocM, decreasing = decrease)]) / (1:length(postH0)))
  #FNR <- rev(cumsum((1-postH0)[order(1-assocM,decreasing=decrease)])) / (nPair - 1:length(postH0))
  #Se <- cumsum((1-postH0)[order(assocM,decreasing=decrease)]) / (sum(1-postH0))
  #Sp <- rev(cumsum(postH0[order(1-assocM,decreasing=decrease)])) / (nPair - sum(1-postH0))
  
  # Number of signals according to the decision rule 
  if (detectCriter =="FDR") nSig <- sum(FDR <= criterThres)
  if (detectCriter == "nSig") nSig <- min(criterThres,nPair)
  if (detectCriter == "assocMeasure") {
    if (assocMeasure =="postH0") nSig <- sum(assocM <= criterThres, na.rm = TRUE)
    if (assocMeasure =="lb05" | assocMeasure=="postE") nSig <- sum(assocM >= criterThres, na.rm = TRUE)
  }
  
  # output ------------------------------------------------------------------
  allSig <- data.frame( 
    dLab[order(assocM, decreasing=decrease)],
    aeLab[order(assocM, decreasing=decrease)],
    n[order(assocM, decreasing=decrease)],
    expN[order(assocM, decreasing=decrease)],
    assocM[order(assocM, decreasing=decrease)],
    (n/expN)[order(assocM, decreasing=decrease)],
    drugMargin[order(assocM, decreasing=decrease)],
    aeMargin[order(assocM, decreasing=decrease)],
    FDR,
    postH0[order(assocM, decreasing=decrease)] )
  
  if (assocMeasure=="postH0") colnames(allSig) <- c("drug","event","n","expected", "postH0","rrr","drugMargin","aeMargin","FDR","postH0")
  if (assocMeasure=="lb05") colnames(allSig) <- c("drug","event","n","expected", "lb05","rrr","drugMargin","aeMargin","FDR","postH0")
  if (assocMeasure=="postE") colnames(allSig) <- c("drug","event","n","expected", "postE","rrr","drugMargin","aeMargin","FDR","postH0")
  
  if (convergence > 2) warning("Prior parameter optimization algorithm did not not converge properly (see ?nlm()). Maybe try another set of values for hyperparamInit")
  res <- vector(mode="list")
  res$sig <- allSig[1:nSig,]
  if (allRes) res$allSig <- allSig  
  res$hyperparam <- hp
  res$convergence <- convergence    
  res$nSig <- nSig
  res$call <- match.call()
  res
}

# gps pvInd ---------------------------------------------------------------
#' @export
gps.pvInd <- function(object, rr0=1, assocMeasure=c("postH0","lb05","postE"), detectCriter=c("FDR","nSig","assocMeasure"), criterThres = 0.05, strat=NULL, nMin=1, truncThres = NULL, hyperparamInit = c(alpha1= 0.2, beta1= 0.06, alpha2=1.4, beta2=1.8, w=0.1), hyperparam = NULL, allRes=F, ...){
  if(!inherits(object, "pvInd")) stop("x must be a pvInd object.")
  assocMeasure <- match.arg(assocMeasure) #keep only the first argument
  detectCriter <- match.arg(detectCriter)
  objectPvCont <- pvCont(object)
  if (!is.null(strat)) {
    temp <- pvCont(object, strat)
    print(head(temp@expN))
    expN <-apply(temp@expN, 1, sum)
    objectPvCont@expN <- matrix(expN, ncol=1)
  }  
  #return(gps(objectPvCont, rr0=rr0, assocMeasure=assocMeasure, detectCriter=detectCriter, criterThres=criterThres, nMin=nMin, truncThres=truncThres, hyperparamInit=hyperparamInit, hyperparam=hyperparam, allRes=allRes))
  NextMethod(generic="gps", object=objectPvCont, rr0=rr0, assocMeasure=assocMeasure, detectCriter=detectCriter, criterThres=criterThres, nMin=nMin, truncThres=truncThres, hyperparamInit=hyperparamInit, hyperparam=hyperparam, allRes=allRes)
}



.quantGps <- function(Seuil, Q, a1, b1, a2, b2) {
  m <- rep(-100000,length(Q))
  M <- rep(100000,length(Q))
  x <- rep(1,length(Q))
  .cdfGps <- function(p,Seuil,Q,a1,b1,a2,b2) {
    Q*pgamma(p,shape=a1,rate=b1)+(1-Q)*pgamma(p,shape=a2,rate=b2)-Seuil
  }
  cdfGpsRes  <- .cdfGps(x,Seuil,Q, a1, b1, a2, b2)
  while (max(round(cdfGpsRes*1e4))!=0)  {
    S <- sign(cdfGpsRes)
    xnew <- (1+S)/2*((x+m)/2)+(1-S)/2*((M+x)/2)
    M <- (1+S)/2*x+(1-S)/2*M
    m <- (1+S)/2*m+(1-S)/2*x
    x <- xnew
    cdfGpsRes <- .cdfGps(x, Seuil, Q, a1, b1, a2, b2)
  }
  x
}
