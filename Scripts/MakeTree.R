library(pheatmap)
library(stringr)
library(RColorBrewer)
library(igraph)
library(ape)
library(ggtree)
library(naturalsort)
library(ggplot2)
library(ggnewscale)

rm(list=ls());gc()

dev.off()
source("Scripts/functions.R")

printPdf=F
aniPath="Data/ani/ANIm_hadamard.tab"
treePath="Data/roary_out/Bsub_reduced.tree"
ucPath="Data/all_phages.uc70_def2"
meshPath="Data/meshout70.txt"
minCluster=3

k=dir(path = "Data/",pattern = "uc", full.names = T)[1]

for(k in dir(path = "Data/",pattern = "uc", full.names = T)) {
  
  #read in all data
  g(ANI, meta) %=% BuildData(aniPath = aniPath, ucPath = k,meshPath = meshPath,minCluster = minCluster, useMesh = F)#, subs=1:100)
  
  
  
  #make presence/absence matrix for all isolates
  isoDF=makeIsoDF(meta)
  meta$ClusterFix
  
  #make annotation matrix mapping phages to phylogeny
  annot=AnnotateClusters(blastDir = "Data/blast_out_genomes70/")
  
  #read in raw tree
  tree0=read.tree(treePath)
  
  #tips that are removed in tree, keep for transparency
  badTips=c("CP018295.1", "CP009611.1", "CP009749.1", "CP009748.1", "CP014840.1")
  #tree=drop.tip(tree0, badTips)
  
  #overwrite tree
  tree=tree0
  
  #fix labels to match with meta
  tree$tip.label=gsub("PS-","",tree$tip.label)
  
  
  #note that tree is smaller than isolates 
  checkDF1=data.frame(naturalsort(rownames(isoDF)), c(naturalsort(tree$tip.label),rep("",17)))
  
  #sort isoDF by tree labels
  #may be handled by ggtree, but for safety
  isoDF1=(isoDF[rownames(isoDF) %in% tree$tip.label,])
  
  missingIso=rownames(isoDF[!(rownames(isoDF) %in% tree$tip.label),])
  
  #missingIso[!(missingIso %in% badTips)]
  
  
  isoDF2=isoDF1[,order(colSums(isoDF1),decreasing = T)]
  
  isoDF3=(makeBin(isoDF2))
  
  HC=hclust(dist(t((isoDF3))),"ward.D2")
  
  isoDF4=isoDF3[,HC$order]
  
  whichSingl=which(colnames(isoDF4)=="singleton")
  
  isoDF5=isoDF4[,c(whichSingl,(1:NCOL(isoDF4))[-whichSingl])]
  
  totPha=data.frame(rowSums(isoDF2))
  colnames(totPha) ="TotalPhage"
  
  notInTree=tree$tip.label [!(tree$tip.label %in% rownames(isoDF1)) ]
  
  checkDF2=data.frame(c(naturalsort(rownames(isoDF1)),""), c(naturalsort(tree$tip.label)))
  
  tree1=drop.tip(tree,notInTree)
  
  tree1=root(tree1, "CP041015.1")
  
  isoMeta=data.frame("Source"=factor(ifelse(nchar(tree1$tip.label)>5,"NCBI","River")))
  rownames(isoMeta)= tree1$tip.label
  
  #Tree
  basicTree = ggtree(tree1,layout = "rect") + geom_treescale() +
    geom_tiplab(colour = factor(ifelse(nchar(tree1$tip.label)>5,1,2)),
                size = 2,
                linesize = 0.005,
                align = TRUE, linetype = "dashed")
  
  treeHeat1 <- gheatmap(basicTree, data = (isoDF5),width = 2, offset = 0.01 ,font.size = 3, low = "white",high = "red",
                        color = "black", colnames = T) + scale_colour_viridis_d(option="D", name="discrete\nvalue") + theme(legend.position = "none")#+ xlim(NA, 24)
  
  treeHeat2 <- treeHeat1 + new_scale_fill()
  
  treeHeat3=gheatmap(treeHeat2, data = subset(meta,select =  Start), offset=.15, width=.1, 
                     colnames_angle=0)  + theme(legend.position = "none")# +scale_fill_manual(values=startCols)
  
  
  
  treeHeat4=gheatmap(treeHeat2, totPha, offset=.15, width=.1, 
                     colnames_angle=0)  + theme(legend.position = "none")
  
  
  if(printPdf) pdf( paste("Figs/BacTree", gsub("Data/all_phages.uc(..)_def2","_\\1",k),".pdf",sep=""), width=16,height = 16)
  treeHeat3
  if(printPdf) dev.off()
  
}
