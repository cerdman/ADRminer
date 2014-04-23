#' @encoding UTF-8
#' @name pvCont-class
#' @aliases show,pvCont-method names,pvCont-method
#' @docType class
#' @title pvCont class
#' @author Youness Ergaibi & Ismaïl Ahmed
#' @exportClass pvCont

# pvCont Class definition --------------------------------------------------
setClass(
  "pvCont",
  representation(drugLab="factor", aeLab="factor", n="matrix", drugMargin="matrix", aeMargin="matrix", expN="matrix", N="numeric", strat="character", coord="matrix"),
  prototype(drugLab=factor(), aeLab=factor(), n=matrix(), drugMargin=matrix(0), aeMargin=matrix(0), expN=matrix(0), N=0, strat=character(), matrix(0))
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
    if (length(object@strat)) cat("Strates: ", paste (object@strat, " ; "), "\n" )
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

