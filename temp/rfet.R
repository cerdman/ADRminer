#' \code{rfet} function
#' @title Reporting Fisher's Exact Test
#' @param object is an object(pvInd or pvCont)
#' @param  or0 is a value of the tested odds ratio. By default  or0=1
#' @param detectCriter is a Decision rule for the signal generation based on:
#' FDR: (Default value) 
#' nbSig:Number of signals 
#' pvalue:Ranking statistic. See \code{pvalue}
#' @param  criterThres is a Threshold for \code{detectCriter}. Ex 0.05 for FDR.
#' @param MID.PVAL the statistic of interest becomes the mid-P-values instead of the P-values resulting from the Fisher's exact test. By default \code{MID.PVAL=FALSE}.
#' @param nMin is a Minimum number of notifications for a couple to be potentially considered as a signal. By default, \code{nMin=1}.
#' @param allRes is a logical value
#' @description Reporting Odds Ratio proposed by van Puijenbroak et al. (2002) extended to the multiple comparison framework.
#' @author Youness Ergaibi & Ismaïl Ahmed
#' @return allSig Data.frame summarizing the results of all couples with at least \code{nMin} notifications ordered by \code{assocMeasure}. It contains notably the labels, the cell counts, the expected counts n1. * n.1 / N
#' @return sig Same Data.frame as \code{allSig} but restricted to the list of generated signals.
#' @return nbSig Number of generated signals.
#' @return call Parameters entered in the function.
#' @keywords RFET

`rfet` <-
  function(object,  or0 = 1, detectCriter=c("FDR","nbSig","pvalue"), criterThres = 0.05, MID.PVAL=FALSE, nMin=1, allRes=F) {
    
 

if(!is.PvCont(object) & !is.PvInd(object)) stop("object must be of PvInd or PvCont class")
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
  logRFET<-log(n)+log(n00)-log(n10)-log(n01)
  
 
  pval.fish.uni <- vector(length=Nb.Cell)
  for (p in 1 : Nb.Cell) {
    pval.fish.uni[p] <-  fisher.test(matrix(c(n[p],n10[p],n01[p],n00[p]),ncol=2,byrow=TRUE),or= or0,alternative="g")$p.value
  }
  
  if (MID.PVAL == TRUE) { # option MID.PVAL
    for (p in 1 : Nb.Cell) {
      pval.fish.uni[p] <- pval.fish.uni[p] - 0.5 * dnoncenhypergeom(x = n[p], n1 = n[p] + n01[p] , n2 = n10[p] + n00[p], m1 = n[p] + n10[p], psi =  or0)
    }
  }
  
  pval.uni <- pval.fish.uni
  pval.uni[pval.uni>1] <-1
  pval.uni[pval.uni<0] <-0
  
  stat <- pval.uni
  LBE.res <- LBE(2 * apply(cbind(stat, 1-stat),1,min),plot.type="none")
  pi.c <- LBE.res$pi0
  
  fdr <- pi.c  * sort(stat[stat <= .5]) / (c(1:sum(stat <= .5)) / Nb.Cell)
  fdr <- c(fdr,
           pi.c /(2 *((sum(stat <= .5)+1) : Nb.Cell)/ Nb.Cell)
           + 1 
           - sum(stat <= .5)/((sum(stat <= .5)+1):Nb.Cell)
  )
  FDR <- apply(cbind(fdr,1),1,min)


  if(detectCriter=="FDR")  nbSig <- sum(FDR <= criterThres, na.rm=T)
    if(detectCriter == "nbSig") nbSig <- min(criterThres,Nb.Cell)
     if(detectCriter == "pvalue") nbSig <- sum(stat <= criterThres, na.rm = TRUE)
 
     
     ############################ SORTIE DE LA FONCTION #############################
     res <- vector(mode="list")
     if(allRes==T){
       res$call <- match.call()
         res$allSig <- data.frame( drugLab[order(stat)],
                                 aeLab[order(stat)],
                                 n[order(stat)],
                                 expN[order(stat)],
                                 stat[order(stat)],
                                 exp(logRFET)[order(stat)],
                                 dMargin[order(stat)],
                                 aeMargin[order(stat)],
                                 FDR )
       colnames(res$allSig) <- c("drug code","event effect","count","expected count","p-value", "RFET","drug margin","event margin","FDR")
       res$sig<- res$allSig[1:nbSig,]
       res$nSig<- nbSig
       res
      
     }
     else
     {
       res$call <- match.call()
       res$sig <- data.frame( drugLab[order(stat)],
                                 aeLab[order(stat)],
                                 n[order(stat)],
                                 expN[order(stat)],
                                 stat[order(stat)],
                                 exp(logRFET)[order(stat)],
                                 dMargin[order(stat)],
                                 aeMargin[order(stat)],
                                 FDR )
       colnames(res$sig) <- c("drug code","event effect","count","expected count","p-value", "RFET","drug margin","event margin","FDR")
       res$sig <-  res$sig[1:nbSig,]
       res
     }
  }