data(sr1)
data(covSr1)
head(covSr1)
covSr1 <- data.frame(covSr1, sexFactor=factor(covSr1$sex, exclude = NULL))  ## exclude=NULL allows to consider NA values as a factor level
summary(covSr1$sexFactor)
ageFactor <- cut(covSr1$age, breaks = c(0, 3, 10, 18, 40, 65, 110), include.lowest = T)
ageFactor <- factor(ageFactor, exclude = NULL)
covSr1 <- data.frame(covSr1, ageFactor)
summary(covSr1$ageFactor)

pvIndSr1 <- pvInd(adr = sr1, cov = covSr1)
pvIndSr1


pvIndNaRm <- naRm(pvIndSr1)

pvIndSr1Red <- pvIndResize(pvIndNaRm, aeMarginMin = 100, drugMarginMin = 50)

pvPen2(pvIndSr1Red, aeId = 1, parallel=F)
