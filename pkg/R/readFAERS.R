#' @encoding UTF-8
#' @title read data from FAERS
#' @description \code{readFAERS} makes it possible to import and convert data from the AERS into a pvInd object. Data are available on the website of the FDA. At present, we only consider initial report (I_F_COD==I). We also keep only data for which drugs do not correspond to a valid trade name (VAL_VBM==1). For more details, please refer to the documentation provided with the FAERS files (ASC_NTS.Doc)
#' @param drugFile path to the DRUG***.txt ASCII file 
#' @param reacFile path to the REAC***.txt ASCII file
#' @param demoFile path to the DEMO***.txt ASCII file
#' @param ROLE_COD a character vector corresponding to the drug's reported role in event. Only reports for drug with the reported role are kept in the pvInd object.
#' \itemize{
#' \item PS: Primary Suspect Drug
#' \item SS: Secondary Suspect Drug
#' \item C: Concomitant
#' \item I: Interacting
#' }
#' @param OCCP_COD a character vector corresponding to abbreviation for the reporter's type of occupation. Only reports for drug with the reported role are kept in the pvInd object.
#' \itemize{
#' \item MD: Physician
#' \item PH: Pharmacist
#' \item OT: Other health-professional
#' \item LW: Lawyer
#' \item CN: Consumer
#' }
#' The default value is set to c("MD","OT", "PH"). 
# @param id identifier to be used in the ASCII files
#' @param aeMarginMin Minimum number of reports in which an ae must be involved in order to be kept
#' @param drugMarginMin Minimum number of reports in which a drug must be involved in order to be kept
#' @author Ismaïl Ahmed
#' @export
#' 

readFAERS <-function(drugFile, reacFile, demoFile, ROLE_COD=c("PS", "SS"), 
                     OCCP_COD=c("MD","OT", "PH"), aeMarginMin = 5, drugMarginMin = 5){
  drug <- as.data.frame(fread(drugFile, sep="$", header=T, colClasses = "character"))
  reac <- as.data.frame(fread(reacFile, sep="$", header=T, colClasses = "character"))
  demo <- as.data.frame(fread(demoFile, sep="$", header=T, colClasses = "character"))
  
  i_f_code <- demo$i_f_code
  demo <- demo[, c("primaryid", "age", "age_cod", "gndr_cod", "reporter_country", "occp_cod")] 
  #print(dim(demo))
  demo <- demo[i_f_code =="I",]
  #print(dim(demo))
  demo <- demo[demo$occp_cod %in% OCCP_COD,]
  #print(dim(demo))
   
  val_vbm <- drug$val_vbm
  drug <- drug[, c("primaryid", "role_cod", "drugname")]
  #print(dim(drug))
  drug <- drug[val_vbm == "1",]
  #print(dim(drug))
  drug <- drug[drug$role_cod %in% ROLE_COD, ]
  #print(dim(drug))
  drug <- drug[drug$primaryid %in% demo$primaryid, ]
  #print(dim(drug))
  drug <- drug[, c("primaryid", "drugname")]
  reac <- reac[, c("primaryid", "pt")]
  #print(dim(reac))
  reac <- reac[reac$primaryid %in% demo$primaryid, ]
  #print(dim(reac))
  
  demo$age <- as.numeric(demo$age)
  demo$age[demo$age_cod == "DEC"] <- demo$age * 10
  demo$age[demo$age_cod == "DY"] <- demo$age / 365  
  demo$age[demo$age_cod == "HR"] <- demo$age / (365 *24)
  demo$age[demo$age_cod == "MON"] <- demo$age / 12
  demo$age[demo$age_cod == "WK"] <- demo$age / 52
  
  res <- pvInd(drug = drug, ae = reac, cov = demo)
  if ((aeMarginMin > 1) || (drugMarginMin > 1)) 
      res <- pvIndResize(res, aeMarginMin = aeMarginMin, drugMarginMin = drugMarginMin)
  res
}