#' @encoding UTF-8
#' @name pvCont-class
#' @aliases show,pvCont-method
#' @aliases names,pvCont-method
#' @docType class
#' @title pvCont class
#' @author Youness Ergaibi & Ismaïl Ahmed
#' @exportClass pvCont

# pvCont Class definition --------------------------------------------------
setClass(
  "pvCont",
  representation(drugLab="factor", aeLab="factor", n="matrix", drugMargin="matrix", aeMargin="matrix", expN="matrix", N="numeric", strat="charOrN", coord="matrix"),
  prototype(drugLab=factor(), aeLab=factor(), n=matrix(), drugMargin=matrix(0), aeMargin=matrix(0), expN=matrix(0), N=0, strat=NULL, matrix(0))
)

# names -------------------------------------------------------------------

setMethod("names", signature(x = "pvCont"), function(x) return(slotNames(x)))# end names method pvCont


# show ------------------------------------------------
setMethod(
  "show",
  "pvCont",
  function (object){
    cat("S4 class: ", as.character(class(object)),"\n")
    cat(nlevels(object@drugLab), " drugs \n" )
    cat(nlevels(object@aeLab), " adverse events \n" )
    if (!is.null(object@strat)) cat("Strates: ", paste (object@strat, " ; "), "\n" )
    #     cat("adr pair counts: ")
    #     cat("@n:\n")
    #     print(head(data.frame(drugLab=object@drugLab, aeLab=object@aeLab, round(object@n,2))))
    #     cat("\n")
    #     cat("expected number: ")
    #     cat("@expN:\n")
    #     print(head(data.frame(drugLab=object@drugLab, aeLab=object@aeLab, round(object@expN,2))))
    #     cat("\n")
    #     cat("@drugMargin:\n")
    #     print(head(data.frame(drugLab=object@drugLab, aeLab=object@aeLab, object@drugMargin)))
    #     cat("\n")
    #     cat("@aeMargin:\n")
    #     print(head(data.frame(drugLab=object@drugLab, aeLab=object@aeLab, object@aeMargin)))
    #     cat("@N:\n")
    #     print(str(object@N))    
  }
)# end show method for pvCont


# # summary method definition ----------------------------------------------------
# setMethod ("summary", signature(object="pvCont"), function(object, ...){
#   nD<-nlevels(object@drugLab)
#   nAe<-nlevels(object@aeLab)
#   N<-sum(object@n)
#   
#   object@drugMargin <- pvCont(object@drugLab,object@aeLab,object@n)@drugMargin    # marges des lignes
#   object@aeMargin <- pvCont(object@drugLab,object@aeLab,object@n)@aeMargin   # marges des colonnes
#   object@expN<-pvCont(object@drugLab,object@aeLab,object@n)@expN
#   cat("nD =", nD ,"\n" )  
#   cat("nAe =", nAe , "\n" ) 
#   cat("N =",N,"\n" )
#   cat("drugMargin=",object@drugMargin,"\n" )   
#   cat("\naeMargin=",object@aeMargin,"\n" ) 
#   cat("\nexpN=", object@expN, "\n" ) 
# })#end  summary method for pvCont

# getdrugLab method definition ----------------------------------------------------
#setGeneric ( "getdrugLab", function(object){standardGeneric("getdrugLab")})
#setMethod("getdrugLab","pvCont",function(object )
#{return (object@drugLab) } )#end getdrugLab method for pvCont

# getaeLab method definition ----------------------------------------------------
#setGeneric ( "getaeLab", function(object){standardGeneric("getaeLab")})

#setMethod("getaeLab","pvCont",function(object )
#{return (object@aeLab) } )#end getaeLab method for pvCont
#####################################################################################
