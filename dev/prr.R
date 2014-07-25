#' @encoding UTF-8
#' @name prr
#' @title Proportional Reporting Ratio
#' @param object An object of class PvInd or PvCont. 
#' @param rr0 The tested relative risk. By default, rr0=1.
#' @param assocMeasure Statistic used to order the drug-events pairs: 
#' pvalue: P-value 
#' Lb95: Lower bound of the 95\% two sided confidence interval of log(prr).
#' @param detectCriter Decision rule for the signal generation based on:
#' \itemize{
#' \item FDR: a prespecified value (\code{detectCriter}) for the false discovery rate (Default value) 
#' \item nbSig: a prespecified number (\code{detectCriter}) of signals
#' \item assocMeasure: a prespecified value (\code{detectCriter}) for the Ranking statistic (\code{assocMeasure})
#' }
#' @param  criterThres Threshold for \code{detectCriter}. Ex 0.05 for FDR.
#' @param nMin Minimum number of spontaneous reports for a drug-event pair to be potentially considered as a signal. By default, \code{nMin}=1.
#' @param strat Character vector containing the name of the covariates to be used for a stratified analysis
#' @param allRes Logical value indicating whether all drug-event combination results must be provided.
#' @description Proportional Reporting Ratio initially proposed by Evans et al. (2001) extended to the multiple comparison framework. Note that the computed variance is different from the one used in van Puijenbroek et al. (2002). The rule proposed by Evans et al. is not implemeted.
#' @author Youness Ergaibi & Ismaïl Ahmed
#' @return allSig Data.frame summarizing the results for all drug-event combinations with at least \code{nMin} spontaneous reports ordered by \code{assocMeasure}. This output is provided if \code{allRes=TRUE}. Operating characteristics are estimated according to Ahmed et al (Biometrics 2010).
#' @return sig Same as \code{allSig} but restricted to the list of generated signals.
#' @return nSig Number of generated signals.
#' @return call Arguments entered when calling \code{prr}.
#' @keywords prr
`prr` <-
  function(object, rr0 = 1, assocMeasure=c("pvalue","Lb95"), detectCriter=c("FDR","nbSig","assocMeasure"), criterThres = 0.05, nMin=1, strat=NULL,allRes=F) {
    
    if(!is.PvCont(object)&!is.PvInd(object)) stop("object must be of pvInd or pvCont class")
    if(is.null(strat)){
      
      if(is.PvCont(object))object<-object
      if(is.PvInd(object))object<-pvInd2pvCont(object)
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
      
      
      logPRR<-log(n)-log(n+n10)-log(n01)+log(n01+n00)
      id<-which(logPRR==Inf)
      logPRR[id] <- NA
      
      varlogPRR <- 1/n - 1/(n10+n) + 1/n01 - 1/(n01+n00)
      id<-which(varlogPRR==Inf)
      varlogPRR[id] <- NA
      
      
    }
    else
    {
      
      if(is.PvCont(object))stop("on peut pas stratifier car l'objet est de type pvCont")
      if(is.PvInd(object)){
        temp<-pvInd2pvContStrat(object,strat)
        #print(temp)
        L<-length(temp)
        v=0
        w=0
        for(l in 1:L){
          n<-temp[[l]]@n
          #print(n)
          dMargin<-temp[[l]]@dMargin
          # print(dMargin)
          aeMargin<-temp[[l]]@aeMargin
          expN<-temp[[l]]@expN
          drugLab<-temp[[l]]@drugLab
          aeLab<-temp[[l]]@aeLab
          N <- sum(n)
          # print(N)
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
          
          
          
          logPRRs<-log(n)-log(n+n10)-log(n01)+log(n01+n00)
          id<-which(logPRRs==Inf)
          logPRRs[id] <- NA
          
          
          varlogPRRs<- 1/n - 1/(n10+n) + 1/n01 - 1/(n01+n00)
          id<-which(varlogPRRs==Inf)
          varlogPRRs[id] <- NA
          v<-((1/varlogPRRs )*logPRRs)+v
          w<-(1/varlogPRRs)+w #sum w
        }
        varlogPRR<-1/w
        logPRR<- varlogPRR *v
        
        
        
      }
      
      
    }
    
    
    pval.logPRR.uni <- 1-pnorm(logPRR,log(or0),sqrt(varlogPRR))
    stat <- (logPRR-log(or0))/sqrt(varlogPRR) # on va trier les signaux par rapport aux valeurs centrées réduites
    pval.uni <- pval.logPRR.uni
    
    pval.uni[pval.uni>1] <-1
    pval.uni[pval.uni<0] <-0
    id<-which(is.na(pval.uni)==T)
    pval.uni[id] <- 1
    
    LBE.res <- LBE(2 * apply(cbind(  pval.uni, 1- pval.uni),1,min),plot.type="none")
    pi.c <- LBE.res$pi0
    fdr <- pi.c  * sort(pval.uni[pval.uni <= .5]) / (c(1:sum(pval.uni <= .5)) / Nb.Cell)
    fdr <- c(fdr,
             pi.c /(2 *((sum(pval.uni <= .5)+1) : Nb.Cell)/ Nb.Cell)
             + 1 
             - sum(pval.uni <= .5)/((sum(pval.uni <= .5)+1):Nb.Cell)
    )
    FDR <- apply(cbind(fdr,1),1,min)
    FDR[id] <- NA
    #id<-which(assocMeasure=="Lb95")
    if(assocMeasure=="Lb95") FDR <- rep(NA,length(n))
    
    # Calculation of the Lower Bound
    LB <- qnorm(0.025,logPRR,sqrt(varlogPRR))
    
    id<-which(assocMeasure=="Lb95")
    idx<-which(assocMeasure=="pvalue")
    if(length(idx)==1) stat <- pval.uni
    if(length(id)==1) stat<- LB
    
    # Calculation of the number of signals according to the decision rule (pval/FDR/Nb of Signals)
    
    if((detectCriter=="FDR") && (assocMeasure=="pvalue"))
      nbSig <- sum(FDR <= criterThres, na.rm=T)
    if(detectCriter == "nbSig") nbSig <- min(criterThres,Nb.Cell)
    if(detectCriter == "assocMeasure") {
      if(assocMeasure=="pvalue") nbSig <- sum(stat <= criterThres, na.rm = TRUE)
      if(assocMeasure=="Lb95") nbSig <- sum(stat >= criterThres, na.rm = TRUE)
    }
    ############################ SORTIE DE LA FONCTION #############################
    res <- vector(mode="list")
    if(allRes==T){
      res$call <- match.call()
      
      
      decrease <- ifelse(assocMeasure=="Lb95", TRUE, FALSE)
      
      res$allSig <- data.frame( drugLab[order(stat, decreasing=decrease)],
                                aeLab[order(stat, decreasing=decrease)],
                                n[order(stat, decreasing=decrease)],
                                expN[order(stat, decreasing=decrease)],
                                stat[order(stat, decreasing=decrease)],
                                exp(logPRR)[order(stat, decreasing=decrease)],
                                dMargin[order(stat, decreasing=decrease)],
                                aeMargin[order(stat, decreasing=decrease)],
                                FDR )
      colnames(res$allSig) <- c("drug code","event effect","count","expected count","p-value", "ROR","drug margin","event margin","FDR")
      if (assocMeasure=="Lb95") colnames(res$allSig)[5] <- "Lb95logPRR" 
      
      
      res$sig<- res$allSig[1:nbSig,]
      
      
      
      # Number of signals
      res$nSig<- nbSig
      res
    }
    else
    {
      decrease <- ifelse(assocMeasure=="Lb95", TRUE, FALSE)
      #res$INPUT.PARAM <- data.frame(or0, nMin,  assocMeasure, criterThres, detectCriter)
      res$call <- match.call()
      res$sig<-data.frame( drugLab[order(stat, decreasing=decrease)],
                           aeLab[order(stat, decreasing=decrease)],
                           n[order(stat, decreasing=decrease)],
                           expN[order(stat, decreasing=decrease)],
                           stat[order(stat, decreasing=decrease)],
                           exp(logPRR)[order(stat, decreasing=decrease)],
                           dMargin[order(stat, decreasing=decrease)],
                           aeMargin[order(stat, decreasing=decrease)],
                           FDR ) 
      colnames(res$sig) <- c("drug code","event effect","count","expected count","p-value", "ROR","drug margin","event margin","FDR")
      res$sig <-  res$sig[1:nbSig,]
      res
    }
  }

