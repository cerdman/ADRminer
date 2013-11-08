
# pvInd -------------------------------------------------
setMethod("names", signature(x = "pvInd"), function(x) return(slotNames(x)))# end names method for pvInd

# pvIndTime -------------------------------------------------
#setMethod("names", signature(x = "pvIndTime"), function(x) return(slotNames(x)))# end names method for pvInd

# PvCont -------------------------------------------------
setMethod("names", signature(x = "pvCont"), function(x) return(slotNames(x)))# end names method pvCont

