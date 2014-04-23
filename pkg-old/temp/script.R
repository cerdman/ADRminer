dataR <- read.table("~/Desktop/adrminerMat/95-03ATC-pt-indivV2.txt", header=T)
                    
pvD <- pvInd(dataR[,1:3])
resGPS <- gps(pvD, criterThres = 0.20)

#pvD2 <- pvIndResize(pvD, 50, 50)

selAe <- unique(as.character(resGPS$sig[,2]))

selD <- vector("list", length=length(selAe))
postH0 <- vector("list", length=length(selAe))
for (i in 1:length(selAe)){
  print(i)
  selD[[i]] <- unique(as.character(resGPS$sig[resGPS$sig[,2] == selAe[i],1]))
  postH0[[i]] <- resGPS$sig$postH0[resGPS$sig[,2] == selAe[i]]
}

resCvGlmnetPos <- vector("list", length=length(selAe))
resCvGlmnetPos2 <- vector("list", length=length(selAe))
resCvGlmnet <- vector("list", length=length(selAe))
for (i in 1:10){
  print(i)
  if (length(selD[[i]]) >= 10){
    print(length(selD[[i]]))
    x <- pvD@drug[,selD[[i]]]
    y <- pvD@ae[, selAe[i]]
    foldid <- sample(1:10, length(y), replace=T)
    resCvGlmnetPos[[i]] <- cv.glmnet(x, y, family="binomial", lower.limits=0, foldid= foldid)
    #resCvGlmnetPos2[[i]] <- cv.glmnet(x, y, family="binomial", lower.limits=0, foldid= foldid, penalty.factor=postH0[[i]])
    resCvGlmnet[[i]] <- cv.glmnet(x, y, family="binomial", foldid=foldid)
  }else{
    resCvGlmnetPos[[i]] <- NA
  }
}


ptime <- system.time({

  foldid <- sample(1:10, 117516, replace=T)
cvG <- cv.glmnet(x=pvD2@drug, y=pvD2@ae[,1], family="binomial", dfmax=50, lower.limits=0, foldid=foldid)
cvGno <- cv.glmnet(x=pvD2@drug, y=pvD2@ae[,1], family="binomial", dfmax=50, foldid=foldid)
idx <- which(cvG$lambda== cvG$lambda.min)
adaPen <- cvG$glmnet.fit$beta[,idx]
xRed <- pvD2@drug[,-which(is.infinite(adaPen))]

adaPenRed <- adaPen[-which(is.infinite(adaPen))]
  xRed <- scale(xRed, center=F, scale=1/adaPenRed)
cvG2 <- cv.glmnet(x=xRed, y=pvD2@ae[,1], family="binomial")
})
toto <- glmnet(x=xRed, y=pvD2@ae[,1], family="binomial", lower.limits=0,  penalty.factor=adaPenRed)

require(randomForest )
toto <- randomForest(x=as.matrix(xRed), y=as.vector(pvD2@ae[,1]))


xRed <- pvD@drug[,as.character((resGPS$sig[which(resGPS$sig[,2]=="10015832"),1]))]
y <- pvD@ae[,"10015832"]
cvG <- cv.glmnet(x=pvD2@drug, y=pvD2@ae[,1], family="binomial", dfmax=50, lower.limits=0, foldid=foldid)
cvGno <- cv.glmnet(x=pvD2@drug, y=pvD2@ae[,1], family="binomial", dfmax=50, foldid=foldid)

## different selection stratégies
## gps pour faire un premier screening assez lache
## rafinement avec le lasso
## utiliser les postH0 de GPS comme poids pour le weighted lasso
## adaptative lasso en deux étapes mais je n'arrive pas à le faire fonctionner
