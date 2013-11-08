# Virtual classes ---------------------------------------------------------
#' @name virtualClass
#' @docType class
#' @aliases factOrN-class numOrN-class dfOrN-class charOrN-class
#' @title Virtual classes for ADRminer
#' @description These virtual classes are only for internal use 

setClassUnion("factOrN", c("factor","NULL"))
setClassUnion("numOrN", c("integer","numeric","NULL"))
setClassUnion("dfOrN", c("data.frame","NULL"))
setClassUnion("charOrN", c("character","NULL"))
