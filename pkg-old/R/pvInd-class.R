# pvInd ---------------------------------------------------------------
#' @encoding UTF-8
#' @title Formal class "pvInd"
#' @description The class pvInd is a class (S4) for storing individual spontaneous report data. 
#' @name pvInd-class
#' @docType class
#' @title pvInd class
#' @author Youness Ergaibi & Ismaïl Ahmed
#' @exportClass pvInd
#' @aliases names,pvInd-method
#' @aliases shows,pvInd-method
#' @aliases getDrug,pvInd-method
#' @aliases getAe,pvInd-method
#' @aliases getDrugMargin,pvInd-method
#' @aliases getAeMargin,pvInd-method
#' @aliases getCov,pvInd-method
#' @aliases $,pvInd-method
#' @aliases $<-,pvInd-method
#' @slot drug Drug matrix in sparse format \link{\code{Matrix}}.
#' @slot ae Adverse Event matrix in sparse format \link{\code{Matrix}}.
#' @slot drugMargin Vector of the marginal counts for each drug.
#' @slot aeMargin Vector of the marginal counts for each ae.
#' @slot cov \code{data.frame} containing individual covariates.
#' @section Methods:
#' \describe{
#'  \item{names}{\code{signature(x = "pvInd")}: returns the names of the slots of the object.}
#'  \item{show}{\code{signature(x = "pvInd")}: printing of the object.}
#'  \item{$}{\code{signature(x = "pvInd")}: similar to the  @@ operator; used to access the content of slots of the object.}
#'  \item{$<-}{\code{signature(x = "pvInd")}: similar to the @@ operator; used to replace the content of slots of the object.}
#'  \item{getDrug}{\code{signature(x = "pvInd")}: returns the names of the slots of the object.}
#'  \item{getAe}{\code{signature(x = "pvInd")}: returns the adverse event matrix.}
#'  \item{getDrugMargin}{\code{signature(x = "pvInd")}: returns the Drug margin counts.}
#'  \item{getAeMargin}{\code{signature(x = "pvInd")}: returns the adverse event margin counts.}
#'  \item{getCov}{\code{signature(x = "pvInd")}: returns the covariate data.frame.}
#' }
#### to be changed if applied to several classes
#' @aliases getDrug getAe getDrugMargin getAeMargin getCov

# pvInd Class definition --------------------------------------------------
setClass(
  "pvInd",
  representation(drug="dgCMatrix", ae="dgCMatrix", drugMargin="numOrN", aeMargin="numOrN", cov="dfOrN"),
  prototype(drug=sparseMatrix(1,1,x=0), ae=sparseMatrix(1,1,x=0), drugMargin=numeric(), aeMargin=numeric(), cov=NULL)
)

.validPvInd <- function(object){
  if (nrow(object@drug) != nrow(object@ae)){
    cat("\n The Number of observations in drug and ae matrix are not equal \n")
    return(FALSE)
  }
  if (length(object@drugMargin) != ncol(object@drug)){
    cat("\n Length of @drugMargin has to equal ncol(@drug) \n")
    return(FALSE)
  }
  if (length(object@aeMargin) != ncol(object@ae)){
    cat("\n Length of @aeMargin has to equal ncol(@ae) \n")
    return(FALSE)
  }
  return(TRUE)
}
setValidity("pvInd", .validPvInd)


# names -------------------------------------------------------------------
setMethod("names", signature(x = "pvInd"), function(x) return(slotNames(x)))# end names method for pvInd


# show --------------------------------------------------------------------
setMethod(
  "show",
  "pvInd",
  function (object){
    cat("S4 class:", as.character(class(object)), "\n")    
    cat("\n@drug: ", nrow(object@drug), "x", ncol(object@drug), ", Drug sparse matrix:\n" )
    print(as.matrix(head(object@drug[,1:6])))        
    cat("@ae: ", nrow(object@ae), "x", ncol(object@ae), ", AE sparse matrix:\n" )
    print(as.matrix(head(object@ae[,1:6])))  
    #cat("@drugMargin: (length=", length(object@drugMargin), ")", head(object@drugMargin), "\n" , sep="")   
    #cat("@aeMargin: (length=", length(object@aeMargin), ")", head(object@aeMargin), "\n" , sep="")          
    if(!is.null(object@cov)) {
      cat("@cov: Covariate data.frame:\n")  
      print(head(object@cov))
    }
  }
)# end show method for pvInd


# 
# # Accessors ---------------------------------------------------------------

setGeneric("getDrug", function(x, ...) standardGeneric("getDrug"))
setGeneric("getAe", function(x, ...) standardGeneric("getAe"))
setGeneric("getDrugMargin", function(x, ...) standardGeneric("getDrugMargin"))
setGeneric("getAeMargin", function(x, ...) standardGeneric("getAeMargin"))
setGeneric("getCov", function(x, ...) standardGeneric("getCov"))
# 
## getDrug
#' @export
setMethod("getDrug","pvInd", function(x,...){
  return(x@drug)
})

## getAe
#' @export
setMethod("getAe","pvInd", function(x,...){
  return(x@ae)
})

## getDrugMargin
#' @export
setMethod("getDrugMargin","pvInd", function(x,...){
  return(x@drugMargin)
})

## getAeMargin
#' @export
setMethod("getAeMargin","pvInd", function(x,...){
  return(x@aeMargin)
})

## getCov
#' @export
setMethod("getCov","pvInd", function(x,...){
  return(x@cov)
})
# 
# ## $
#' @export
setMethod("$","pvInd",function(x, name) {
return(slot(x,name))
})
# 

## $<-
#' @export
setMethod("$<-","pvInd",function(x, name, value) {
  slot(x,name,check=TRUE) <- value
  return(x)
})
# 


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
# setMethod("getDrugMargin","pvInd", function(object){return(object@drugMargin)})#end getDrugMargin method for pvInd

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
