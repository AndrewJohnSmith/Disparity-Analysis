#Set your working directory
setwd("/Users/Home/Dropbox/Amherst/Publications/Evolution2014/Dryad/Scripts")

##Load in a few key packages...
#...or install them using install.packages(c("abind", "reshape", "vegan", "spam", "fields", "plotrix", "gdata"))
require(abind)
require(reshape)
require(vegan)
require(spam)
require(fields)
require(plotrix)
require(gdata)


###########################
#Correction Factor Scripts#
###########################
#Import correction factor functions.
source("CorrectionFunctions.Source.R")

#Read in your data - must follow the same format.
#If taxa span multiple bins use | symbol.
	#For example, if a taxon resides in bin 1, 2 and 3, then use 1|2|3 in the bin row (row 2).
LuoExtinct<-read.csv("Luo2007Example.csv",header=F)

##Diparity Correction factor##
correctionFactor(Data=LuoExtinct, numSims=100, axes=10)

##Disparity without correction##
DisparityNoCorr(Data=LuoExtinct, axes=10, iterations=10)

###T-test###
PCOMatrix<-read.csv("PCOMatrixExample.csv", header = TRUE)
DisparityTTest(PCOMatrix, replicates=1000)


######################
#Data Removal Scripts#
######################
#Import data removal functions.
 #Note that the CorrectionFunctions.R source cannot also be loaded.
source("Ordination-Disparity.Source.R")

#Import the extinct and extant data matricies.
importExtant <- as.matrix(read.csv("LuoMatrixExtant.csv", na.strings="NA"))
importExtinct <- as.matrix(read.table("LuoMatrixExtinct.txt", header=T, sep="\t"))

#If you want to simulate random loss switch linkage to F and use set all the characters to 0.
importExtinct<-as.matrix(read.table("LuoMatrixExtinctRandom.txt", header=T, sep="\t"))

##Ordination Spaces##
#Calculate PCO ordination space
importPCO <- calculate.pco(importExtant, axes=10, "manhattan")

#Save the PCO matrix?
write.csv(importPCO, "PCO.NOLoss.csv")
#Once saved it can be read in at any time using the line below.
importPCO<-read.csv("PCO.NOLoss.csv", header=T, row.names=1)

#Calculate the difference between the orginal and 'data-loss' ordination spaces.
PCOsimulateAndanalyze(extinctDataset=importExtinct, extantDataset=importExtant, importPCO=importPCO, numSims=100, charsToLose=5, softChars=2, raw.dist.method="manhattan", pco.dist.method="euclidean", axes=10)

#Plot the differences in a 2D morphospace. The PCO.Sim data is obtained by running the previous function.
plot.pco(pco.data=PCO.Sim, no.loss.pco=importPCO, Add=F, xlimit=c(-100,100), ylimit=c(-100,100), colour="red")

##Disparity##
DisparitysimulateAndanalyze(extinctDataset=importExtinct, extantDataset=importExtant, numSims=100, charsToLose=5, softChars=2, axes=10)