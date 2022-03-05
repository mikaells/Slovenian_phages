library(pheatmap)
library(stringr)
library(RColorBrewer)
library(igraph)
library(naturalsort)

rm(list=ls())


source("Scripts/functions.R")

printPdf=T
aniPath="Data/ani/ANIm_hadamard.tab"
ucPath="Data/all.uc70"
minCluster=3


g(ANI, meta) %=% BuildData(aniPath = aniPath, ucPath = ucPath,minCluster = minCluster)#, subs=1:100)



isoDF_forScoary=makeIsoDF(meta = meta)

rownames(isoDF_forScoary)=ifelse(nchar(rownames(isoDF_forScoary))<5, paste("PS.",rownames(isoDF_forScoary), sep=""),rownames(isoDF_forScoary))


GPA=read.table("Data/gene_presence_absence.csv", sep=",", header = T)


colnames(GPA)

GPA1=GPA[,1:14]

GPA2=GPA[,15:NCOL(GPA)]

keepRows = (rownames(isoDF_forScoary) %in% colnames(GPA2))

isoDF_forScoary2=makeBin(isoDF_forScoary[keepRows,])

#rownames(isoDF_forScoary2) =gsub("PS\\.","PS-",rownames(isoDF_forScoary2))

data.frame(c(" ",naturalsort(rownames(isoDF_forScoary2))), naturalsort(colnames(GPA2)))

rownames(isoDF_forScoary2) %in% colnames(GPA2)

GPA3=GPA2[,(colnames(GPA2) %in% rownames(isoDF_forScoary2))]

all(colnames(GPA3) %in% rownames(isoDF_forScoary2))

all(rownames(isoDF_forScoary2) %in% colnames(GPA3) )

isoDF_forScoary3=data.frame("Name"=rownames(isoDF_forScoary2),isoDF_forScoary2[,-1])


GPA4=data.frame(GPA1, GPA3)

write.csv(GPA4, file = "Data/gene_presence_absence_fixed.csv",  row.names = F)

write.csv(isoDF_forScoary3, "Data/PhagePresAbs.csv", row.names = F)
