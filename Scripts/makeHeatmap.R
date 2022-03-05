library(pheatmap)
library(stringr)
library(RColorBrewer)
library(igraph)

rm(list=ls())


source("Scripts/functions.R")

printPdf=T
aniPath="Data/ani/ANIm_hadamard.tab"
ucPath="Data/all_phages.uc60_def2"
meshPath="Data/meshout70.txt"
minCluster=3

k=dir(path = "Data/",pattern = "uc", full.names = T)[1]
for(k in dir(path = "Data/",pattern = "uc", full.names = T)) {
  
  g(ANI, meta) %=% BuildData(aniPath = aniPath, ucPath = k,minCluster = minCluster, meshPath = meshPath, useMesh = F)#, subs=1:200)

  cols=brewer.pal(12, "Paired")
  
  #make a list for variables and their colors 
  
  #first define the colors to use
  annot_cols=list(
    Type = 1:2,
    Date = 1:2,
    #there are N colors to make, so we make a gradient across here from the 12 colors defined in cols
    ClusterFix = color.gradient(1:length(unique(meta$ClusterFix)),cols)
  )
  
  #add names to the colors so pheatmap knows what to pick
  names(annot_cols$Type)=c("plasmid","phages")
  names(annot_cols$Date)=c("NCBI","River bank")
  names(annot_cols$ClusterFix)=naturalsort::naturalsort(unique(meta$ClusterFix))
  
  whichSingl=which(names(annot_cols$ClusterFix)=="singleton")
  #flip around so that singletons are first
  
  annot_cols$ClusterFix=c(annot_cols$ClusterFix[whichSingl], annot_cols$ClusterFix[-whichSingl])
  
  # annot_cols$ClusterFix[c(1,whichSingl)] = annot_cols$ClusterFix[c(whichSingl,1)]
  # names(annot_cols$ClusterFix)[c(1,whichSingl)]=names(annot_cols$ClusterFix)[c(whichSingl,1)]
  
  
  
  #make singletons grey
  annot_cols$ClusterFix[1]="grey70"
  
  
  #create cluster objects for heatmap on binary data
  binClusterMat=ifelse(ANI>0,1,0)
  
  colClusterMat=hclust(dist(t(binClusterMat)),"ward.D2")
  rowClusterMat=hclust(dist(binClusterMat),"ward.D2")
  

  
  if(printPdf) pdf(paste("Figs/PhageHeatmapBinClust", gsub("Data/all_phages.uc(..)_def2","_\\1",k),".pdf",sep=""), width=16,height = 16)
  PH=pheatmap(ANI,  fontsize_col = 6,fontsize_row = 6, fontsize = 10, cluster_rows = rowClusterMat, cluster_cols = colClusterMat,
              display_numbers = F,labels_row = meta$ID,  annotation_colors = annot_cols, cutree_cols = 20,cutree_rows = 20,
              number_format = "%.0f", annotation_row = subset(meta, select = c("logLen","Start", "Type", "Date", "ClusterFix"))  )
  if(printPdf) dev.off()
  
}
