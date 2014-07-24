drugFile <- c("/Users/isma/Dropbox/Work/Calculs/Pharmaco/Data/FAERS/faers_ascii_2012q4/ascii/drug12q4.txt")
reacFile <- c("/Users/isma/Dropbox/Work/Calculs/Pharmaco/Data/FAERS/faers_ascii_2012q4/ascii/reac12q4.txt")
demoFile <- c("/Users/isma/Dropbox/Work/Calculs/Pharmaco/Data/FAERS/faers_ascii_2012q4/ascii/demo12q4.txt")

faers12q4 <- readFAERS(drugFile, reacFile, demoFile)

drugFile <- c("/Users/isma/Dropbox/Work/Calculs/Pharmaco/Data/FAERS/faers_ascii_2012q4/ascii/drug12q4.txt.zip")
reacFile <- c("/Users/isma/Dropbox/Work/Calculs/Pharmaco/Data/FAERS/faers_ascii_2012q4/ascii/reac12q4.txt.zip")
demoFile <- c("/Users/isma/Dropbox/Work/Calculs/Pharmaco/Data/FAERS/faers_ascii_2012q4/ascii/demo12q4.txt.zip")

faers12q4 <- readFAERS(drugFile, reacFile, demoFile)

fread(input = "https://www.dropbox.com/s/x04kkfh8vegrm22/drug12q4.txt")

drugFile <- c("/Users/isma/Dropbox/Work/Calculs/Pharmaco/Data/FAERS/faers_ascii_2013q1/asci/DRUG13Q1.txt")
reacFile <- c("/Users/isma/Dropbox/Work/Calculs/Pharmaco/Data/FAERS/faers_ascii_2013q1/asci/REAC13Q1.txt")
demoFile <- c("/Users/isma/Dropbox/Work/Calculs/Pharmaco/Data/FAERS/faers_ascii_2013q1/asci/DEMO13Q1.txt")

faers13q1 <- readFAERS(drugFile, reacFile, demoFile)

drugFile <- c("/Users/isma/Dropbox/Work/Calculs/Pharmaco/Data/FAERS/faers_ascii_2013q2/asii/DRUG13Q2.txt")
reacFile <- c("/Users/isma/Dropbox/Work/Calculs/Pharmaco/Data/FAERS/faers_ascii_2013q2/asii/REAC13Q2.txt")
demoFile <- c("/Users/isma/Dropbox/Work/Calculs/Pharmaco/Data/FAERS/faers_ascii_2013q2/asii/DEMO13Q2.txt")

faers13q2 <- readFAERS(drugFile, reacFile, demoFile)

drugFile <- c("/Users/isma/Dropbox/Work/Calculs/Pharmaco/Data/FAERS/AERS_ASCII_2013q3/ASCII/DRUG13Q3.txt")
reacFile <- c("/Users/isma/Dropbox/Work/Calculs/Pharmaco/Data/FAERS/AERS_ASCII_2013q3/ASCII/REAC13Q3.txt")
demoFile <- c("/Users/isma/Dropbox/Work/Calculs/Pharmaco/Data/FAERS/AERS_ASCII_2013q3/ASCII/DEMO13Q3.txt")

faers13q3 <- readFAERS(drugFile, reacFile, demoFile)
all <- rbind(faers12q4, faers13q1, faers13q2, faers13q3)

allResize <- pvIndResize(all, aeMarginMin = 50, drugMarginMin = 50)

ageFac <- cut(allResize$cov$age, c(0,1,10,25,50, 120), include.lowest = T)
ageFac <- addNA(ageFac)
allResize$cov$ageFac <- ageFac
allResize$cov$gndr_cod <- factor(allResize$cov$gndr_cod, )
levels(allResize$cov$gndr_cod)
levels(allResize$cov$gndr_cod) <- c("UNK", "F", "M", "UNK", "UNK")

#allResizeNA <- naRm(allResize) 

resGPS <- gps(allResize, strat = c("gndr_cod", "ageFac"), )
resPvPen <- pvPen2(object = allResize, aeId = "Death", nDrugMax = 100)
resPvPen2 <- pvPen(object = allResize, aeId = "Death", nDrugMax = 100,covId = c("ageFac", "gndr_cod"))



