#' \code{bcpnn} function
#' @title Bayesian confidence propagation neural network
#' @param object is an object(pvInd or pvCont)
#' @param rr0 is the Value of the tested risk. By default, RR0=1.
#' @param assocMeasure is a Statistic used for ranking the couples: 
#' post.H0:Posterior probability of the null hypothesis
#' LB:2.5% quantile of the posterior distribution of IC.
#' @param detectCriter is a Decision rule for the signal generation based on:
#' FDR: (Default value) 
#' nbSig:Number of signals 
#' assocMeasure:Ranking statistic. See \code{assocMeasure}
#' @param  criterThres is a Threshold for \code{detectCriter}. Ex 0.05 for FDR.
#' @param nMin is a Minimum number of notifications for a couple to be potentially considered as a signal. By default, \code{nMin=1}.
#'  @param NB.MC  If MC=TRUE, NB.MC indicates the number of Monte Carlo simulations to be done
#'@param  MC If MC=TRUE, the statistic of interest (see RANKSTAT) is calculated by Monte Carlo simulations which can be very long. If MC=FALSE, IC is approximated by a normal distribution (which can be very crude for small counts).
#' @param allRes is a logical value
#' @description Bayesian confidence propagation neural network (Bate et al. 1998, Noren et al. 2006) extended to the multiple comparison framework.
#' @author Youness Ergaibi & Ismaïl Ahmed
#' @return allSig Data.frame summarizing the results of all couples with at least \code{nMin} notifications ordered by \code{assocMeasure}. It contains notably the labels, the cell counts, the expected counts n1. * n.1 / N
#' @return sig Same Data.frame as \code{allSig} but restricted to the list of generated signals.
#' @return nbSig Number of generated signals.
#' @return call Parameters entered in the function.
#' @keywords bcpnn

`bcpnn` <-
  function(object, rr0 = 1, assocMeasure=c("post.H0","LB"), detectCriter=c("FDR","nbSig","assocMeasure"), criterThres = 0.05, nMin=1, MC=FALSE, NB.MC=10000,allRes=TRUE) {
    
    if(!is.PvCont(object) & !is.PvInd(object)) stop("object must be of PvInd or PvCont class")
    if(is.PvInd(object)) object<-pvInd2pvCont(object)
    n<-object@n
    dMargin<-object@dMargin
    aeMargin<-object@aeMargin
    expN<-object@expN
    drugLab<-object@drugLab
    aeLab<-object@aeLab
    N <- sum(n)
    n10 <- dMargin - n
    n01 <- aeMargin - n
    n00 <- N - (n+n10+n01)
    if( nMin>1) {
      if(nMin > max(n)) stop ("nMin is too large")
      expN <- expN[n >= nMin]
      dMargin <- dMargin[n >= nMin]
      aeMargin <- aeMargin[n >= nMin]
      n10 <- n10[n >= nMin]
      n01 <- n01[n >= nMin]
      n00 <- n00[n >= nMin]
      
      n<- n[n >= nMin]
    }
    Nb.Cell <- length(n)
    
   
    if (MC == FALSE) {
      post.H0 <- matrix(nrow=Nb.Cell,ncol=length(rr0))
      p1  <- 1 +  dMargin 
      p2  <- 1 + N -  dMargin 
      q1  <- 1 +  aeMargin 
      q2  <- 1 + N -  aeMargin 
      r1  <- 1 + n
      r2b <- N - n -1 + (2+N)^2/(q1*p1)
      EICb <- log(2)^(-1) * (digamma(r1) - digamma(r1+r2b) - (digamma(p1) - digamma(p1+p2) + digamma(q1) - digamma(q1+q2)))
      VICb <- log(2)^(-2) * (trigamma(r1) - trigamma(r1+r2b) + (trigamma(p1) - trigamma(p1+p2) + trigamma(q1) - trigamma(q1+q2)))
      post.H0 <- pnorm(log(rr0),EICb,sqrt(VICb))
      # Calculation of the Lower Bound
      LB <- qnorm(0.025,EICb,sqrt(VICb))
    }
    
    if (MC == TRUE) { # Advanced option MC
      require(MCMCpack)
      dMargin <- n + n10
      aeMargin <- n + n01
      Nb_Obs <- length(n)
      
      ## Nouvelles priors
      q1. <- (dMargin +.5)/(N +1)
      q.1 <- (aeMargin +.5)/(N +1)
      q.0 <- (N - aeMargin +.5)/(N +1)
      q0. <- (N - dMargin +.5)/(N +1)
      
      a.. <- .5/(q1.*q.1) ## le .5 devrait pouvoir être changé
      
      a11 <- q1.*q.1* a..
      a10 <- q1.*q.0* a..
      a01 <- q0.*q.1* a..
      a00 <- q0.*q.0* a..
      
      g11 <- a11 + n
      g10 <- a10 + n10
      g01 <- a01 + n01
      g00 <- a00 + n00
      g1. <- g11 + g10
      g.1 <- g11 + g01
      
      post.H0 <- vector(length=length(n))
      LB <- vector(length=length(n))
      quantile <- vector("numeric",length=length(n))
      for (m in 1 : length(n)){
        p <- rdirichlet(NB.MC,c(g11[m],g10[m],g01[m],g00[m]))
        p11 <- p[,1]
        p1. <- p11 + p[,2]
        p.1 <- p11 + p[,3]	
        IC_monte <- log(p11/(p1.* p.1))
        temp <- IC_monte < log(rr0)
        post.H0[m] <- sum(temp)/NB.MC
        LB[m] <- sort(IC_monte)[round(NB.MC * 0.025)]
      }
      rm(p11,p1.,p.1,IC_monte,temp)
      gc()
    }
    if (assocMeasure=="post.H0") assocM <- post.H0
    if (assocMeasure=="LB") assocM <- LB
    
   
    if (assocMeasure=="post.H0") {
      FDR <- (cumsum(post.H0[order(assocM)]) / (1:length(post.H0)))
      FNR <- rev(cumsum((1-post.H0)[order(1-assocM)])) / (Nb.Cell - 1:length(post.H0))
      Se <- cumsum((1-post.H0)[order(assocM)]) / (sum(1-post.H0))
      Sp <- rev(cumsum(post.H0[order(1-assocM)])) / (Nb.Cell - sum(1-post.H0))
    }

    if (assocMeasure=="LB") {
      FDR <- (cumsum(post.H0[order(assocM,decreasing=TRUE)]) / (1:length(post.H0)))
      FNR <- rev(cumsum((1-post.H0)[order(1-assocM,decreasing=TRUE)])) / (Nb.Cell - 1:length(post.H0))
      Se <- cumsum((1-post.H0)[order(assocM,decreasing=TRUE)]) / (sum(1-post.H0))
      Sp <- rev(cumsum(post.H0[order(1-assocM,decreasing=TRUE)])) / (Nb.Cell - sum(1-post.H0))
    }
    
   
    
    if (detectCriter =="FDR") nbSig<- sum(FDR <= criterThres)
    if (detectCriter == "nbSig") nbSig <- min(criterThres,Nb.Cell)
    if (detectCriter == "assocMeasure") {
      if (assocMeasure=="post.H0") nbSig <- sum(assocM <= criterThres, na.rm = TRUE)
      if (assocMeasure=="LB") nbSig <- sum(assocM >= criterThres, na.rm = TRUE)
    }
    res <- vector(mode="list")
    if(allRes==T){
      ############################ Output #############################
      res$call <- match.call()
      if (assocMeasure=="post.H0") {
        res$allSig <- data.frame( drugLab[order(assocM)],
                                  aeLab[order(assocM)],
                                  n[order(assocM)],
                                  expN[order(assocM)],
                                  assocM[order(assocM)],
                                  (n/expN)[order(assocM)],
                                  dMargin[order(assocM)],
                                  aeMargin[order(assocM)],
                                  FDR, FNR, Se, Sp )
        colnames(res$allSig) <- c("drug","event","count","expected count","postH0",
                                  "n/expN","drug margin","event margin","FDR","FNR","Se","Sp")
      }
      if (assocMeasure=="LB" ) {
        res$allSig <- data.frame( drugLab[order(assocM, decreasing=TRUE)],
                                  aeLab[order(assocM, decreasing=TRUE)],
                                  n[order(assocM, decreasing=TRUE)],
                                  expN[order(assocM,decreasing=TRUE)],
                                  assocM[order(assocM,decreasing=TRUE)],
                                  (n/expN)[order(assocM,decreasing=TRUE)],
                                  dMargin[order(assocM,decreasing=TRUE)],
                                  aeMargin[order(assocM,decreasing=TRUE)],
                                  FDR, FNR, Se, Sp,
                                  post.H0[order(assocM,decreasing=TRUE)] )
    colnames(res$allSig) <- c("drug","event","count","expected count",
                                                          "Q_0.025(log(IC))","n/expN","drug margin","event margin","FDR","FNR","Se","Sp","postH0")
      
  }
      
      res$sig <- res$allSig[1:nbSig,]
      
      
      # Number of signals
      res$nSig <- nbSig
      res
    }
    else
    {
      res$call <- match.call() 
      if (assocMeasure=="post.H0") {
        res$sig <- data.frame( drugLab[order(assocM)],
                               aeLab[order(assocM)],
                               n[order(assocM)],
                               expN[order(assocM)],
                               assocM[order(assocM)],
                               (n/expN)[order(assocM)],
                               dMargin[order(assocM)],
                               aeMargin[order(assocM)],
                               FDR, FNR, Se, Sp )
        colnames(res$sig) <- c("drug","event","count","expected count","postH0",
                               "n/expN","drug margin","event margin","FDR","FNR","Se","Sp")
      }   
      
      
      if (assocMeasure=="LB" ) {
        res$sig <- data.frame( drugLab[order(assocM, decreasing=TRUE)],
                               aeLab[order(assocM, decreasing=TRUE)],
                               n[order(assocM, decreasing=TRUE)],
                               expN[order(assocM,decreasing=TRUE)],
                               assocM[order(assocM,decreasing=TRUE)],
                               (n/expN)[order(assocM,decreasing=TRUE)],
                               dMargin[order(assocM,decreasing=TRUE)],
                               aeMargin[order(assocM,decreasing=TRUE)],
                               FDR, FNR, Se, Sp,
                               post.H0[order(assocM,decreasing=TRUE)] )
 colnames(res$sig) <- c("drug","event","count","expected count",
                                                       "Q_0.025(log(IC))","n/expN","drug margin","event margin","FDR","FNR","Se","Sp","postH0")
        
        
      }
      
      # List of Signals generated according to the method 
      res$sig <- res$sig[1:nbSig,]
      
      res
      
    }
}

   
