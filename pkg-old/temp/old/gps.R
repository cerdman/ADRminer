#' @encoding UTF-8
#' @name gps
#' @title Gamma Poisson Shrinkage
#' @param object An object of class pvInd or pvCont
#' @param rr0 The tested risk. By default, rr0=1.
#' @param assocMeasure Statistic used to order the drug-events pairs: 
#' \itemize{
#' \item postH0: Posterior probability of the null hypothesis
#' \item LB: 5\% quantile of the posterior distribution of lambda
#' \item postE: Posterior expectation of log(lambda,2)
#' }
#' @param detectCriter Decision rule for the signal generation based on:
#' FDR: (Default value) 
#' nSig: number of signals 
#' assocMeasure: ranking statistic. See \code{assocMeasure}
#' @param criterThres Threshold used for \code{detectCriter}. Ex 0.05 for FDR.
#' @param nMin Minimum number of spontaneous reports for a drug-event pair to be potentially considered as a signal. By default, \code{nMin}=1.
#' @param truncThres Used for the hyper parameter calculations. The marginal likelihood involved mixtures of two (possibly truncated) negative binomial. The default value is NULL 
# in order to account for that drug-event pairs associated with no spontaneous reports are usually not represented in the data. In order not to use any truncation, set \code{truncThres} to NULL
#' @param  priorInit initialization vector for the prior mixture parameters (alpha1, beta1, alpha2, beta2, w). By default, priorInit = c(alpha1 = 0.2, beta1 = 0.06, alpha2 = 1.4, beta2 = 1.8, w = 0.1), i.e. the prior parameters provided in DuMouchel (1999).
#' @param priorParam Hyper parameter values for the 2 gamma mixture model. By default, priorParam = NULL which means that the hyperparameters are calculated by maximising the marginal likelihood. This parameter is meant to avoid this maximization step
#' @param strat Character vector containing the name of the covariates to be used for a stratified analysis.
#' @param allRes Logical value indicating whether all drug-event combination results must be provided. 
#' @description Gamma Poisson Shrinkage model proposed by DuMouchel (1999) extended to the multiple comparison framework.
#' @author Youness Ergaibi & Ismaïl Ahmed
#' @return allSig Data.frame summarizing the results for all drug-event combinations with at least \code{nMin} spontaneous reports ordered by \code{assocMeasure}. This output is provided if \code{allRes=TRUE}. Operating characteristics are estimated according to Ahmed et al (Stat Med 2009). 
#' @return sig Same as \code{allSig} but restricted to the list of generated signals.
#' @return nSig Number of generated signals.
#' @return call Arguments entered when calling \code{gps}.
#' @keywords gps
#' 

gps <- function(object, rr0 = 1, assocMeasure=c("postH0","LB","postE"), detectCriter=c("FDR","nSig","assocMeasure"), criterThres = 0.05, nMin=1, truncThres = NULL, priorInit = c(alpha1= 0.2, beta1= 0.06, alpha2=1.4, beta2=1.8, w=0.1), priorParam = NULL, strat=NULL, allRes=F) {
  if(!is.PvCont(object) & !is.PvInd(object)) stop("object must be of PvInd or PvCont class")
  if(is.PvInd(object)) objectPvCont<-pvInd2pvCont(object)
  n<-objectPvCont@n
  dMargin<-objectPvCont@dMargin
  aeMargin<-objectPvCont@aeMargin
  expN<-objectPvCont@expN
  drugLab<-objectPvCont@drugLab
  aeLab<-objectPvCont@aeLab
  
  # Algorithm allowing to keep only the combinations associated with at least nMin sr
  if(nMin > 1) {
    expN <- expN[n >= nMin]
    dMargin <- dMargin[n >= nMin]
    aeMargin <- aeMargin[n >= nMin]
    n<- n[n >= nMin]
  }       
  if (!is.null(strat)){
    if(is.PvCont(object))stop("Stratified analysis cannot be performed with aggregated data (PvCont) ")
    temp <- pvInd2pvContStrat(object,strat)
    L <- length(temp)
    expN <- 0
    # print(L)
    for(l in 1:L){
      expN <- expN + temp[[l]]@expN
    }
  }
  nPair <- length(n)
  postH0 <- vector(length=nPair)
  
  if (is.null(priorParam)){         
    likTrunc <-function(p, n, expN, truncThres){
      sum(
        -log(
          p[5] * ( dnbinom(n, size=p[1], prob=p[2]/(p[2]+expN)) / pnbinom(trunk, size=p[1], prob=p[2]/(p[2]+expN)) )
          + (1-p[5]) * ( dnbinom(n, size=p[3], prob=p[4]/(p[4]+expN))/ pnbinom(trunk, size=p[1], prob=p[2]/(p[2]+expN)) )
        )
      , na.rm=T)
    }
    lik <-function(p, n, expN){
      sum(
        -log(
        p[5] * dnbinom(n, size=p[1], prob=p[2]/(p[2]+expN))
        + (1-p[5]) * dnbinom(n, size=p[3], prob=p[4]/(p[4]+expN))
        )
        )
    }
    
    if (is.null(truncThres)){
      pOut <-nlm(lik, p=priorInit, n=n, expN=expN, iterlim=1000)
    }else{
      if (truncThres>=0){
        pOut <-nlm(likTrunc, p=priorInit, n=n[n>=truncThres], expN=expN[n>=truncThres], trunc,iterlim=1000)
      }else{stop("truncThres must be non negative or null")}
    }
    priorParam <- pOut$estimate
    convergence <- pOut$code
  }
  
  ## fast determination of priorParam. discretization
  #expNQ <- quantile(expN, probs=(seq(0,1,0.01))) 
  #expNQM <- (expNQ[-101]+expNQ[-1])/2
  #expNTable <- table(cut(expN , expNQ))
  # Posterior probability of the null hypothesis
  Q <- priorParam[5]*dnbinom(n,size=priorParam[1],prob=priorParam[2]/(priorParam[2]+expN)) / 
    (priorParam[5]*dnbinom(n,size=priorParam[1],prob=priorParam[2]/(priorParam[2]+expN)) + (1-priorParam[5])* 
       dnbinom(n,size=priorParam[3],prob=priorParam[4]/(priorParam[4]+expN)))
  
  postH0 <- Q*pgamma(rr0,priorParam[1]+n,priorParam[2]+expN) +(1-Q)*pgamma(rr0,priorParam[3]+n,priorParam[4]+expN) # proba a posteriori de H0
  
  # Posterior Expectation of Lambda
  postE <- log(2)^(-1)*(Q*(digamma(priorParam[1]+n)-log(priorParam[2]+expN))+(1-Q)*(digamma(priorParam[3]+n)-log(priorParam[4]+expN)))
  
  # Algorithm allowing the calculation of the CI Lower Bound (at Seuil%)    
  # Algorithm Emmanuel Roux
  QuantileDuMouchel<-function(Seuil,Q,a1,b1,a2,b2) {
    m<-rep(-100000,length(Q))
    M<-rep(100000,length(Q))
    x<-rep(1,length(Q))
    Cout<-FCoutQuantileDuMouchel(x,Seuil,Q,a1,b1,a2,b2)
    while (max(round(Cout*1e4))!=0)  {
      S<-sign(Cout)
      xnew<-(1+S)/2*((x+m)/2)+(1-S)/2*((M+x)/2)
      M<-(1+S)/2*x+(1-S)/2*M
      m<-(1+S)/2*m+(1-S)/2*x
      x<-xnew
      Cout<-FCoutQuantileDuMouchel(x,Seuil,Q,a1,b1,a2,b2)
    }
    x
  }
  FCoutQuantileDuMouchel<-function(p,Seuil,Q,a1,b1,a2,b2) {
    Q*pgamma(p,shape=a1,rate=b1)+(1-Q)*pgamma(p,shape=a2,rate=b2)-Seuil
  }
  
  # Calculation of LB
  LB <- QuantileDuMouchel(0.05,Q,priorParam[1]+n,priorParam[2]+expN,priorParam[3]+n,priorParam[4]+expN)   
  
  # Assignment based on the method (pp/postE/LB)
  if (assocMeasure=="postH0") assocM <- postH0
  if (assocMeasure=="LB") assocM <- LB
  if (assocMeasure=="postE") assocM <- postE
  
  # FDR, FNR, Se and Sp based on the method (pp/postE/LB)
  if (assocMeasure=="postH0") {
    FDR <- (cumsum(postH0[order(assocM)]) / (1:length(postH0)))
    FNR <- rev(cumsum((1-postH0)[order(1-assocM)])) / (nPair - 1:length(postH0))
    Se <- cumsum((1-postH0)[order(assocM)]) / (sum(1-postH0))
    Sp <- rev(cumsum(postH0[order(1-assocM)])) / (nPair - sum(1-postH0))
  } 
  
  if (assocMeasure=="LB" | assocMeasure=="postE") {
    FDR <- (cumsum(postH0[order(assocM,decreasing=TRUE)]) / (1:length(postH0)))
    FNR <- rev(cumsum((1-postH0)[order(1-assocM,decreasing=TRUE)])) / (nPair - 1:length(postH0))
    Se <- cumsum((1-postH0)[order(assocM,decreasing=TRUE)]) / (sum(1-postH0))
    Sp <- rev(cumsum(postH0[order(1-assocM,decreasing=TRUE)])) / (nPair - sum(1-postH0))
  }
  
  # Number of signals according to the decision rule 
  if (detectCriter =="FDR") nSig<- sum(FDR <= criterThres)
  if (detectCriter == "nSig") nSig <- min(criterThres,nPair)
  if (detectCriter == "assocMeasure") {
    if (assocMeasure=="postH0") nSig <- sum(assocM <= criterThres, na.rm = TRUE)
    if (assocMeasure=="LB" | assocMeasure=="postE") nSig <- sum(assocM >= criterThres, na.rm = TRUE)
  }
  
  # Output ------------------------------------------------------------------
  
  if (assocMeasure=="postH0") {
    allSig <- data.frame( drugLab[order(assocM)],
                          aeLab[order(assocM)],
                          n[order(assocM)],
                          expN[order(assocM)],
                          assocM[order(assocM)],
                          (n/expN)[order(assocM)],
                          dMargin[order(assocM)],
                          aeMargin[order(assocM)],
                          FDR, FNR, Se, Sp )
    colnames(allSig) <- c("drug","event","n","expected","postH0",
                          "rrr","drugMargin","aeMargin","FDR","FNR","Se","Sp")
  }   
  
  if (assocMeasure=="LB" | assocMeasure=="postE") {
    allSig <- data.frame( drugLab[order(assocM, decreasing=TRUE)],
                          aeLab[order(assocM, decreasing=TRUE)],
                          n[order(assocM, decreasing=TRUE)],
                          expN[order(assocM,decreasing=TRUE)],
                          assocM[order(assocM,decreasing=TRUE)],
                          (n/expN)[order(assocM,decreasing=TRUE)],
                          dMargin[order(assocM,decreasing=TRUE)],
                          aeMargin[order(assocM,decreasing=TRUE)],
                          FDR, FNR, Se, Sp,
                          postH0[order(assocM,decreasing=TRUE)] )
    if (assocMeasure=="LB") colnames(allSig) <- c("drug","event","n","expected",                                                        "Q_0.05(lambda)","rrr","drugMargin","aeMargin","FDR","FNR","Se","Sp","postH0")
    if (assocMeasure=="postE") colnames(allSig) <- c("drug","event","n","expected", "postE(Lambda)","rrr","drugMargin","aeMargin","FDR","FNR","Se","Sp","postH0")
  }
  res <- vector(mode="list")
  res$sig <- allSig[1:nSig,]
  if (allRes) res$allSig <- allSig  
  res$priorParam <- priorParam
  res$convergence <- convergence
  res$call <- match.call()      
  res$nSig <- nSig
  res
}
