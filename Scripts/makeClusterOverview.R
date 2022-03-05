library(pheatmap)
library(stringr)
library(RColorBrewer)
library(igraph)
library(beeswarm)

rm(list=ls())


source("Scripts/functions.R")

writeTable=T
aniPath="Data/ani/ANIm_hadamard.tab"
ucPath="Data/all_phages.uc60_def2"
meshPath="Data/meshout70.txt"
minCluster=3

k=dir(path = "Data/",pattern = "uc", full.names = T)[1]
for(k in dir(path = "Data/",pattern = "uc", full.names = T)) {
  
  uc=read.table(k)
  
  colnames(uc)=c("Type","ClusterID", "Length", "ID%", "Strand", "-","-","Alignment","SequenceID","TargetID" )
  
  uc1=uc[which(uc$Type=="H" | uc$Type=="S"),]
  
  uc2=uc1[order(uc1$ClusterID),]
  
  uc3=uc[which(uc$Type=="C"),]
  
  annot=AnnotateClusters()
  
  #uc4=cbind(uc3,annot)
  
  if(writeTable) {
    write.table(uc2,file=paste("Data/ClusterMembership", gsub("Data/all_phages.uc(..)_def2","_\\1",k),".txt",sep=""), row.names = F, col.names = T)
    
    write.table(uc3,file=paste("Data/ClusterOverview", gsub("Data/all_phages.uc(..)_def2","_\\1",k),".txt",sep=""), row.names = F, col.names = T)
  }
  
}


uc5=uc2
uc5$`ID%`[uc5$`ID%`=="*"] = 100
uc5$`ID%`=100-as.numeric(uc5$`ID%`)

pdf(file = "Figs/ClusterDissimilarities.pdf", width = 16, height = 12, onefile = T )
par(mar=c(2.5,2.5,0.5,0.5), mgp=c(1.4,0.4,0), font.lab=2, mfrow=c(2,1))

plot(sort(uc5$`ID%`), ylab="Dissimilarity", pch=21, bg="grey", cex=.8)

beeswarm( uc5$`ID%`~ uc5$ClusterID, pch=21,las=2, pwbg=color.gradient(uc5$Length),cex=.9, corral="wrap", xlab="ClusterID",ylab="Dissimilarity")

legend("topright",legend = round(seq(min(uc5$Length),max(uc5$Length), length.out=7),0), title = "Length", col = color.gradient(1:7), pch=16)
dev.off()
