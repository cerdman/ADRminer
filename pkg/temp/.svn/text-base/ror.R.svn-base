#' @title Reporting Odds Ratio
#' @param object An object of class PvInd or PvCont. 
#' @param or0 The tested odds ratio. By default, or0 = 1.
#' @param assocMeasure Statistic used to order the drug-events pairs: 
#' \itemize{
#' \item pvalue: P-value 
#' \item Lb95: Lower bound of the 95\% two sided confidence interval of log(ror).
#' }
#' @param detectCriter Decision rule for the signal generation based on
#' \itemize{
#' \item FDR: a prespecified value (\code{detectCriter}) for the false discovery rate (Default value) 
#' \item nbSig: a prespecified number (\code{detectCriter}) of signals
#' \item assocMeasure: a prespecified value (\code{detectCriter}) for the Ranking statistic (\code{assocMeasure})
#' }
#' @param  criterThres Threshold for \code{detectCriter}. Ex 0.05 for FDR.
#' @param nMin Minimum number of spontaneous reports for a drug-event pair to be potentially considered as a signal. By default, \code{nMin}=1.
#' @param strat Character vector containing the name of the covariates to be used for a stratified analysis
#' @param allRes Logical value indicating whether all drug-event combination results must be provided.
#' @description Reporting Odds Ratio extended to the multiple comparison framework.
#' @author Youness Ergaibi & Ismaïl Ahmed
#' @return allSig Data.frame summarizing the results for all drug-event combinations with at least \code{nMin} spontaneous reports ordered by \code{assocMeasure}. This output is provided if \code{allRes=TRUE}. Operating characteristics are estimated according to Ahmed et al (Biometrics 2010).
#' @return sig Same as \code{allSig} but restricted to the list of generated signals.
#' @return nSig Number of generated signals.
#' @return call Arguments entered when calling \code{ror}.
#' @keywords ror

setGeneric(
  name="ror",
  def=function(object, ...){standardGeneric("ror")}
)


`ror` <-
  function(object, or0 = 1, assocMeasure=c("pvalue","Lb95"), detectCriter=c("FDR","nbSig","assocMeasure"), criterThres = 0.05, nMin=1, strat=NULL, allRes=F) {
    
    if(!is.PvCont(object) & !is.PvInd(object)) stop("object must be of PvInd or PvCont class")
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
      logROR<-log(n)+log(n00)-log(n10)-log(n01)
      id<-which(logROR==Inf)
      logROR[id] <- NA
      
      varlogROR <- 1/n + 1/n10 + 1/n01 + 1/n00
      id<-which(varlogROR==Inf)
      varlogROR[id] <- NA
    }
    else{
      if(is.PvCont(object))stop("on peut pas stratifier car l'objet est de type PvCont")
      if(is.PvInd(object)){
        temp<-pvInd2pvContStrat(object,strat)
        #print(temp)
        L<-length(temp)
        v<-0
        w<-0
        # print(L)
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
          #print(Nb.Cell)
          logRORs<-log(n)+log(n00)-log(n10)-log(n01)
          id<-which(logRORs==Inf)
          logRORs[id] <- NA
          #print(logRORs)
          
          varlogRORs<- 1/n + 1/n10 + 1/n01 + 1/n00
          id<-which(varlogRORs==Inf)
          varlogRORs[id] <- NA
          v<-((1/varlogRORs )*logRORs)+v
          w<-(1/varlogRORs)+w #sum wi
        }
        varlogROR<-1/w
        logROR<- varlogROR *v
      
        
        
      }
    }
    pval.logOR.uni <- 1-pnorm(logROR,log(or0),sqrt(varlogROR))
    stat <- (logROR-log(or0))/sqrt(varlogROR) # on va trier les signaux par rapport aux valeurs centrées réduites
    pval.uni <- pval.logOR.uni
    
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
    LB <- qnorm(0.025,logROR,sqrt(varlogROR))
    
    
    if(assocMeasure=="pvalue") stat <- pval.uni
    if(assocMeasure=="Lb95") stat<- LB
    
    # Calculation of the number of signals according to the decision rule (pval/FDR/Nb of Signals)
    
    if((detectCriter=="FDR") && (assocMeasure=="pvalue"))
      nbSig <- sum(FDR <= criterThres, na.rm=T)
    if(detectCriter == "nbSig") nbSig <- min(criterThres,Nb.Cell)
    if(detectCriter == "assocMeasure") {
      if(assocMeasure=="pvalue") nbSig <- sum(stat <= criterThres, na.rm = TRUE)
      if(assocMeasure=="Lb95") nbSig <- sum(stat >= criterThres, na.rm = TRUE)
    }
    res <- vector(mode="list")
    if(allRes==T){
    ############################ SORTIE DE LA FONCTION #############################
   
    #res$INPUT.PARAM <- data.frame(or0, nMin,  assocMeasure, criterThres, detectCriter)
    res$call <- match.call()
    
    
    decrease <- ifelse(assocMeasure=="Lb95", TRUE, FALSE)
    
    res$allSig <- data.frame( drugLab[order(stat, decreasing=decrease)],
                                  aeLab[order(stat, decreasing=decrease)],
                                  n[order(stat, decreasing=decrease)],
                                  expN[order(stat, decreasing=decrease)],
                                  stat[order(stat, decreasing=decrease)],
                                  exp(logROR)[order(stat, decreasing=decrease)],
                                  dMargin[order(stat, decreasing=decrease)],
                                  aeMargin[order(stat, decreasing=decrease)],
                                  FDR )
    colnames(res$allSig) <- c("drug code","event effect","count","expected count","p-value", "ROR","drug margin","event margin","FDR")
    if (assocMeasure=="Lb95") colnames(res$allsig)[5] <- "Lb95LogROR" 
    
    
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
                  exp(logROR)[order(stat, decreasing=decrease)],
                  dMargin[order(stat, decreasing=decrease)],
                  aeMargin[order(stat, decreasing=decrease)],
                  FDR ) 
      colnames(res$sig) <- c("drug code","event effect","count","expected count","p-value", "ROR","drug margin","event margin","FDR")
      res$sig <-  res$sig[1:nbSig,]
      res
    }
  }
