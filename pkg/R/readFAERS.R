#' @encoding UTF-8
#' @title read data from FAERS
#' @description \code{readFAERS} makes it possible to import and convert data from the AERS into a pvInd object. Data are available on the website of the FDA. Currently, we only consider initial report (I_F_COD==I). We also remove data for which drugs do not correspond to a valid trade name (VAL_VBM==1). For more details, please refer to the documentation provided with the FAERS files (ASC_NTS.Doc). First, we try to read the data with the \code{fread} function of the \code{data.table} package. If the latter doesn't work, the data are read with the \code{read.table} function. At present, the cleaning process for the age covariate is a bit drastic: all value < 0 are considered as NA, value associated with missing AGE_COD are considered as NA and ages entered in second are also considered as missing. Then the age is converted into year.
#' @param drugFile path to the DRUG***.txt ASCII file (required)
#' @param reacFile path to the REAC***.txt ASCII file (required)
#' @param demoFile path to the DEMO***.txt ASCII file (required)
#' @param ROLE_COD a character vector corresponding to the drug's reported role in event. Only reports for drug with the reported role are kept in the pvInd object.
#' \itemize{
#' \item PS: Primary Suspect Drug
#' \item SS: Secondary Suspect Drug
#' \item C: Concomitant
#' \item I: Interacting
#' }
#' The default value is to kept only the PS and SS drugs (\code{ROLE_COD=c("PS", "SS")})
#' @param OCCP_COD a character vector corresponding to abbreviation for the reporter's type of occupation. Only reports for drug with the reported role are kept in the pvInd object.
#' \itemize{
#' \item MD: Physician
#' \item PH: Pharmacist
#' \item OT: Other health-professional
#' \item LW: Lawyer
#' \item CN: Consumer
#' }
#' @return an object of class pvInd
#' The default value is set to c("MD","OT", "PH"). 
# @param id identifier to be used in the ASCII files
# @param aeMarginMin Minimum number of reports in which an ae must be involved in order to be kept
# @param drugMarginMin Minimum number of reports in which a drug must be involved in order to be kept
#' @author Ismaïl Ahmed
#' @export
#' 

readFAERS <-function(drugFile, reacFile, demoFile, ROLE_COD=c("PS", "SS"), 
                     OCCP_COD=c("MD","OT", "PH")){ #, aeMarginMin = 5, drugMarginMin = 5){
  drug <- try(as.data.frame(fread(drugFile, sep="$", header=T, colClasses = "character")))
  if (class(drug) == "try-error")
    drug <- read.table(drugFile, sep="$", quote="", comment.char = "", header=T, colClasses = "character")
  reac <- try(as.data.frame(fread(reacFile, sep="$", header=T, colClasses = "character")))
  if (class(reac) == "try-error")
    reac <- read.table(reacFile, sep="$", quote="", comment.char = "", header=T, colClasses = "character")
  demo <- try(as.data.frame(fread(demoFile, sep="$", header=T, colClasses = "character")))
  if (class(demo) == "try-error")
    demo <- read.table(demoFile, sep="$", quote="", comment.char = "", header=T, colClasses = "character", )
  
  #print(head(drug))
  #print(head(reac))
  #print(head(demo))
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
  reac <- reac[reac$primaryid %in% drug$primaryid, ]
  drug <- drug[drug$primaryid %in% reac$primaryid, ]
  demo <- demo[demo$primaryid %in% drug$primaryid, ]
  #print(dim(reac))
  
  demo$age[demo$age==""] <- NA
  #demo$age_cod[demo$age_cod ==""] <- NA
  demo$age[demo$age_cod==""] <- NA
  demo$age[demo$age_cod=="SEC"] <- NA
  demo$age <- as.numeric(demo$age)
  demo$age[demo$age < 0] <- NA
  demo$age[demo$age_cod == "DEC"] <- demo$age[demo$age_cod == "DEC"] * 10
  demo$age[demo$age_cod == "DY"] <- demo$age[demo$age_cod == "DY"] / 365  
  demo$age[demo$age_cod == "HR"] <- demo$age[demo$age_cod == "HR"] / (365 *24)
  demo$age[demo$age_cod == "MON"] <- demo$age[demo$age_cod == "MON"] / 12
  demo$age[demo$age_cod == "WK"] <- demo$age[demo$age_cod == "WK"] / 52
  demo$age[demo$age > 120] <- NA
  demo <- demo[, -3] ## remove age_cod
  
  res <- pvInd(drug = drug, ae = reac, cov = demo)
  #if ((aeMarginMin > 1) || (drugMarginMin > 1)) 
  #    res <- pvIndResize(res, aeMarginMin = aeMarginMin, drugMarginMin = drugMarginMin)
  res
}