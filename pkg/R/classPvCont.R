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
  representation(dLab="factor", aeLab="factor", n="matrix", dMargin="matrix", aeMargin="matrix", expN="matrix", N="numeric", strat="charOrN", coord="matrix"),
  prototype(dLab=factor(), aeLab=factor(), n=matrix(), dMargin=matrix(0), aeMargin=matrix(0), expN=matrix(0), N=0, strat=NULL, matrix(0))
)



# # summary method definition ----------------------------------------------------
# setMethod ("summary", signature(object="pvCont"), function(object, ...){
#   nD<-nlevels(object@drugLab)
#   nAe<-nlevels(object@aeLab)
#   N<-sum(object@n)
#   
#   object@dMargin <- pvCont(object@drugLab,object@aeLab,object@n)@dMargin    # marges des lignes
#   object@aeMargin <- pvCont(object@drugLab,object@aeLab,object@n)@aeMargin   # marges des colonnes
#   object@expN<-pvCont(object@drugLab,object@aeLab,object@n)@expN
#   cat("nD =", nD ,"\n" )  
#   cat("nAe =", nAe , "\n" ) 
#   cat("N =",N,"\n" )
#   cat("dMargin=",object@dMargin,"\n" )   
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
