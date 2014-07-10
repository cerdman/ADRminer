# pvInd ---------------------------------------------------------------
#' @encoding UTF-8
#' @title Formal class "pvInd"
#' @description The class pvInd is a class (S4) for storing individual spontaneous report data. 
#' @name pvInd-class
#' @docType class
#' @author Youness Ergaibi & Ismaïl Ahmed
#' @exportClass pvInd
#' @aliases names,pvInd-method 
#' show,pvInd-method
#' getDrug,pvInd-method 
#' getAe,pvInd-method 
#' getDrugMargin,pvInd-method 
#' getAeMargin,pvInd-method 
#' getCov,pvInd-method 
#' $,pvInd-method 
#' $<-,pvInd-method 
#' naRm 
#' naRm,pvInd-method
#' getDrug
#' getAe
#' getDrugMargin
#' getAeMargin
#' getCov 
#' @slot drug Drug matrix in sparse format \code{Matrix}.
#' @slot ae Adverse Event matrix in sparse format \code{Matrix}.
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
#'  \item{naRm}{\code{signature(x = "pvInd")}: remove spontaneous reports with NA values in the covariate slot.}
#' }

# pvInd Class definition --------------------------------------------------
setClass(
  "pvInd", 
  representation(
    drug="Matrix", 
    ae="Matrix", 
    drugMargin="numeric", 
    aeMargin="numeric", 
    cov="data.frame"
  ),
  prototype(drug=Matrix(0, sparse = T), ae=Matrix(0, sparse = T), drugMargin=numeric(), aeMargin=numeric(), cov=data.frame())
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


# na.omit -----------------------------------------------------------------
#' @export
setGeneric("naRm", function(x, ...) standardGeneric("naRm"))

setMethod(
  "naRm", 
  signature(x = "pvInd"), 
  function(x, ...){
    if (nrow(x@cov)==0) return(x)
    idxRow <- which(is.na(x@cov), arr.ind=T)[,1]
    if (length(idxRow>0)) {
      idxRow <- unique(idxRow)
      x@drug <- x@drug[-idxRow,]
      x@ae <- x@ae[-idxRow,]
      x@cov <- x@cov[-idxRow,]
      nSp <- nrow(x@drug)
      x@drugMargin <- as.numeric(Matrix(1, nrow=1, ncol=nSp) %*% x@drug)
      x@aeMargin <- as.numeric(Matrix(1, nrow=1, ncol=nSp) %*% x@ae)
    }
    return(x)
  }
) # end naRm method for pvInd



# show --------------------------------------------------------------------
setMethod(
  "show",
  "pvInd",
  function (object){
    cat("S4 class:", as.character(class(object)), "\n")    
    cat("\n@drug: ", nrow(object@drug), "x", ncol(object@drug), ", Drug sparse matrix:\n" )
#    show(object@drug)        
    cat("@ae: ", nrow(object@ae), "x", ncol(object@ae), ", AE sparse matrix:\n" )
#    show(object@ae)  
    #cat("@drugMargin: (length=", length(object@drugMargin), ")", head(object@drugMargin), "\n" , sep="")   
    #cat("@aeMargin: (length=", length(object@aeMargin), ")", head(object@aeMargin), "\n" , sep="")          
    cat("@cov: Covariate data.frame:", nrow(object@cov), "x", ncol(object@cov), "\n" )
#    print(head(object@cov))
  }
)# end show method for pvInd


# 
# # Accessors ---------------------------------------------------------------
#' @export
setGeneric("getDrug", function(x, ...) standardGeneric("getDrug"))
#' @export
setGeneric("getAe", function(x, ...) standardGeneric("getAe"))
#' @export
setGeneric("getDrugMargin", function(x, ...) standardGeneric("getDrugMargin"))
#' @export
setGeneric("getAeMargin", function(x, ...) standardGeneric("getAeMargin"))
#' @export
setGeneric("getCov", function(x, ...) standardGeneric("getCov"))

# 
## getDrug

setMethod("getDrug","pvInd", function(x,...){
  return(x@drug)
})

## getAe

setMethod("getAe","pvInd", function(x,...){
  return(x@ae)
})

## getDrugMargin

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

setMethod("$","pvInd",function(x, name) {
  return(slot(x,name))
})
# 

## $<-

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
