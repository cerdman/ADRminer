#  pvInd ------------------------------------------------
setMethod(
  "show",
  "pvInd",
  function (object){
    cat("S4 class:", as.character(class(object)), "\n")    
    cat("\n@drug: ", nrow(object@drug), "x", ncol(object@drug), ", Drug sparse matrix:\n" )
    print(as.matrix(head(object@drug[,1:6])))        
    cat("@ae: ", nrow(object@ae), "x", ncol(object@ae), ", AE sparse matrix:\n" )
    print(as.matrix(head(object@ae[,1:6])))  
    #cat("@dMargin: (length=", length(object@dMargin), ")", head(object@dMargin), "\n" , sep="")   
    #cat("@aeMargin: (length=", length(object@aeMargin), ")", head(object@aeMargin), "\n" , sep="")          
    if(!is.null(object@cov)) {
      cat("@cov: Covariate data.frame:\n")  
      print(head(object@cov))
    }
  }
)# end show method for pvInd

# pvCont ------------------------------------------------
setMethod(
  "show",
  "pvCont",
  function (object){
    cat("S4 class: ", as.character(class(object)),"\n")
    cat(nlevels(object@dLab), " drugs \n" )
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
#     cat("@dMargin:\n")
#     print(head(data.frame(drugLab=object@drugLab, aeLab=object@aeLab, object@dMargin)))
#     cat("\n")
#     cat("@aeMargin:\n")
#     print(head(data.frame(drugLab=object@drugLab, aeLab=object@aeLab, object@aeMargin)))
#     cat("@N:\n")
#     print(str(object@N))    
  }
)# end show method for pvCont
