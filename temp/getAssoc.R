#' \code{getAssoc} function
#' @usage  getAssoc(res,ATC=NULL,meddra=NULL)
#' @param res  must be the output of one the signal detection method function GPS.
#' @param ATC The label of the drug. By default, ATC=FALSE.
#' @param meddra The label of the adverse event. By default, meddra=FALSE.
#' @return ATC Recalls the label of the drug.
#' @return meddra Recalls the label of the event.
#' @return EXIST_ATC  Indicates if the label of the drug exists in the database.
#' @return EXIST_meddra Indicates if the label of the adverse event exists in the database.
#' @return  EXIST_COUPLE Indicates if the couple is present in the database.
#' @return  LIST  It is a dataframe that contains the labels, the counts, the expected counts, the value of the statistic of interest, the rank and the estimated FDR for each couple.
#' @description This function makes possible to extract some information from the output of the ADRminer functions for a given couple adverse event-drug, for a drug or for an adverse event.
#' @author Youness Ergaibi & Ismaïl Ahmed & Antoine Poncet 
#' @keywords getAssoc


`getAssoc` <-
  function(res,ATC=NULL,meddra=NULL) {
    # fonction permettant de rechercher les résultats pour un médicament, un effet indésirable ou les deux
    
    
    res$allSig[,1]<-as.character(res$allSig[,1])
    res$allSig[,2]<-as.character(res$allSig[,2])
  
   
    if (res$call$assocMeasure== "LB" | res$call$assocMeasure == "postE")
      rang <- rank(-res$allSig[,5])
    else rang <- rank(res$allSig[,5])
    
    EXIST_ATC <- TRUE
    EXIST_meddra <- TRUE
    if (is.null(ATC)==FALSE) EXIST_ATC <- length(res$allSig[,1][res$allSig[,1]==ATC])!=0
    if (is.null(meddra)==FALSE) EXIST_meddra <- length(res$allSig[,2][res$allSig[,2]==meddra])!=0
    
    
    if (is.null(meddra)==FALSE & is.null(ATC)==FALSE) {
      Lmedic   <- res$allSig[,1][res$allSig[,1]==ATC & res$allSig[,2]==meddra]
      Leffect  <- res$allSig[,2][res$allSig[,1]==ATC & res$allSig[,2]==meddra]
      effectif <- res$allSig[,3][res$allSig[,1]==ATC & res$allSig[,2]==meddra]
      expect_E <- res$allSig[,4][res$allSig[,1]==ATC & res$allSig[,2]==meddra]
      stat     <- res$allSig[,5][res$allSig[,1]==ATC & res$allSig[,2]==meddra]
      r        <- rang[res$allSig[,1]==ATC & res$allSig[,2]==meddra]
    
      fdr      <- res$allSig[,9][res$allSig[,1]==ATC & res$allSig[,2]==meddra]
    }
    
    if (is.null(meddra) & is.null(ATC)==FALSE) {
      Lmedic   <- res$allSig[,1][res$allSig[,1]==ATC]
      Leffect  <- res$allSig[,2][res$allSig[,1]==ATC]
      effectif <- res$allSig[,3][res$allSig[,1]==ATC]
      expect_E <- res$allSig[,4][res$allSig[,1]==ATC]
      stat     <- res$allSig[,5][res$allSig[,1]==ATC]
      r        <- rang[res$allSig[,1]==ATC]
      fdr      <- res$allSig[,9][res$allSig[,1]==ATC]
    }
    
    if (is.null(meddra)==FALSE & is.null(ATC)) {
      Lmedic   <- res$allSig[,1][res$allSig[,2]==meddra]
      Leffect  <- res$allSig[,2][res$allSig[,2]==meddra]
      effectif <- res$allSig[,3][res$allSig[,2]==meddra]
      expect_E <- res$allSig[,4][res$allSig[,2]==meddra]
      stat     <- res$allSig[,5][res$allSig[,2]==meddra]
      r        <- rang[res$allSig[,2]==meddra]
      #fdr      <- res$OpChar[,1][res$allSig[,2]==meddra]
      fdr      <- res$allSig[,9][res$allSig[,2]==meddra]
    }
    
    
    RES <- vector(mode="list")
    RES$ATC <- ATC
    RES$meddra <- meddra
    RES$EXIST_meddra <- EXIST_meddra
    RES$EXIST_ATC  <- EXIST_ATC
    RES$EXIST_COUPLE <- length(Lmedic)!=0
    RES$LIST <- data.frame(Lmedic,Leffect,effectif,expect_E,stat,r,fdr)
    colnames(RES$LIST) <- c("drug name","event","count","expected count","stat","rank","fdr")
    RES
  }

