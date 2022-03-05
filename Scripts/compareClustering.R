library(pheatmap)
library(stringr)
library(RColorBrewer)
library(igraph)
library(naturalsort)

rm(list=ls())


source("Scripts/functions.R")

clus60=read.table("Data/all_phages.uc60_def2")
clus70=read.table("Data/all_phages.uc70_def2")

colnames(clus60)=c("Type","ClusterID", "Length", "ID%", "Strand", "-","-","Alignment","SequenceID","TargetID" )
colnames(clus70)=c("Type","ClusterID", "Length", "ID%", "Strand", "-","-","Alignment","SequenceID","TargetID" )

clus60_seq=clus60[which(clus60$Type=="H" | clus60$Type=="S"),-8]
clus70_seq=clus70[which(clus70$Type=="H" | clus70$Type=="S"),-8]

clus60_cen=clus60[which(clus60$Type=="C"),]
clus70_cen=clus70[which(clus70$Type=="C"),]

clus60_seq_name=clus60_seq[naturalorder(clus60_seq$SequenceID),]
clus70_seq_name=clus70_seq[naturalorder(clus70_seq$SequenceID),]

both=merge(clus60_seq_name, clus70_seq_name, by = "SequenceID")

both_sort=both[order(both$ClusterID.x),]






