# pvInd ---------------------------------------------------------------
#' @encoding UTF-8
#' @name pvInd-class
#' @docType class
#' @title pvInd class
#' @author Youness Ergaibi & Ismaïl Ahmed
#' @exportClass pvInd

# pvInd Class definition --------------------------------------------------
setClass(
  "pvInd",
  representation(drug="dgCMatrix", ae="dgCMatrix", dMargin="numOrN", aeMargin="numOrN", cov="dfOrN"),
  prototype(drug=sparseMatrix(1,1,x=0), ae=sparseMatrix(1,1,x=0), dMargin=numeric(), aeMargin=numeric(), cov=NULL)
)

.validPvInd <- function(object){
  if (nrow(object@drug) != nrow(object@ae)){
    cat("\n The Number of observations in drug and ae matrix are not equal \n")
    return(FALSE)
  }
  if (length(object@dMargin) != ncol(object@drug)){
    cat("\n Length of @dMargin has to equal ncol(@drug) \n")
    return(FALSE)
  }
  if (length(object@aeMargin) != ncol(object@ae)){
    cat("\n Length of @aeMargin has to equal ncol(@ae) \n")
    return(FALSE)
  }
  return(TRUE)
}
setValidity("pvInd", .validPvInd)

# 
# # pvIndTime ---------------------------------------------------------------
# #' @encoding UTF-8
# #' @name pvIndTime-class
# #' @aliases pvIndTime
# #' @aliases names,pvIndTime-method
# #' @docType class
# #' @title pvIndTime class
# #' @author Ismaïl Ahmed
# #' @exportClass pvIndTime
# #'  
# setClass(
#   "pvIndTime",
#   contains="pvInd",
#   representation(time="ANY")
# )
# 
# .validPvIndTime <- function(object){
#   if (nrow(object@drug) != nrow(object@ae)){
#     cat("\n The numbers of observations in the drug and ae matrices differ \n")
#     return(FALSE)
#   }
#   return(TRUE)
# }
# setValidity("pvIndTime", .validPvIndTime)
# 
# 
# 


# getDrug method definition -----------------------------------------------
# setGeneric ( "getDrug", function(object){standardGeneric("getDrug")})
# setMethod("getDrug","pvInd", function(object){return(object@drug)})#end getdrug method for pvInd

# getAe method definition -------------------------------------------------
# setGeneric ( "getAe", function(object){standardGeneric("getAe")})
# setMethod("getAe","pvInd", function(object){return(object@ae)})#end getAe method for pvInd

# getAe method definition -------------------------------------------------
# setGeneric ( "getDrugMargin", function(object){standardGeneric("getDrugMargin")})
# setMethod("getDrugMargin","pvInd", function(object){return(object@dMargin)})#end getDrugMargin method for pvInd

# getAe method definition -------------------------------------------------
# setGeneric ( "getAeMargin", function(object){standardGeneric("getAeMargin")})
# setMethod("getAeMargin","pvInd", function(object){return(object@ae)})#end getAeMargin method for pvInd


#### Isma : la suite est a revoir ens

# summary method definition -------------------------------------------------
# setMethod ("summary", signature(object="pvInd"), function(object, ...){
#   x <- object
#   if(!is.pvInd(x)) stop("Provided object is not a valid pvInd.")
#   nD <- length(colnames(x@drug))
#   nAe <- length(colnames(x@ae))
#   nInd <- length(rownames(x@drug))#il faut que le nombre individu dans la table drugs soit égale au #nombre  d'individu dans la table  ae
#   nMd <- round(mean(apply(x@drug,1,sum)),digits=2)
#   nMae <- round(mean(apply(x@ae,1,sum)),digits=2)
#   sex <- table(x@sex)
#   age <- summary(x@age)
#   cat("nD =", nD , "nombre des drugs \n" )  
#   cat("nAe =", nAe , "nombre des adverse events \n" ) 
#   cat("nInd =", nInd, "nombre des individus \n" )
#   cat("nMd =", nMd, "nombre moyen des drugs par individu \n" )
#   cat("nMae =", nMae, "nombre moyen des ae par individu \n" )
#   cat("nombre d\'individu par chaque sex est:" )
#   print(sex)
#   cat("le resume d\'age des individus est:\n" )
#   print(age)
#   
# })#end  summary method for pvInd
# 
# # summary method definition -------------------------------------------------
# setMethod ("plot" , signature (x ="pvInd" , y ="ANY") ,
#            function (x , y=NULL , ...){
#              par(mfrow=c(2,1))
#              boxplot(x@age, main="age des individus")
#              pie(table(x@sex),main="le pourcentage de sex des individues")
#            })#end plot method for pvInd
# 
# 
# 
