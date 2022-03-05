library(pheatmap)
library(stringr)
library(RColorBrewer)
library(igraph)

rm(list=ls());gc()


source("Scripts/functions.R")

printPdf=T
aniPath="Data/ani/ANIm_hadamard.tab"
ucPath="Data/all_phages.uc60_def2"
blastDir="Data/blast_out_genomes/"
minCluster=3
addAnnotation=F
meshPath="Data/meshout70.txt"
useMesh=T


g(ANI, meta) %=% BuildData(aniPath = aniPath, ucPath = ucPath,minCluster = minCluster,meshPath = meshPath, useMesh = F)#, subs=1:100)


annot=AnnotateClusters(blastDir = blastDir )

meta$Annotation=""

i=2
for(i in 1:NROW(meta)) {
  
  if(length(which(meta$Cluster[i]==annot$cluster))==0) {
    meta$Annotation[i] =""
  } else {
    meta$Annotation[i] = matchAnnot=annot[which(meta$Cluster[i]==annot$cluster),"bestVote"]
  }
}

meta$AnnotationFixed=gsub("Bacteriophage ","",gsub(" complete genome.","",gsub("[A-Za-z]* phage","",gsub(">[A-Z]{2}[0-9]* ","",meta$Annotation))))

#par(mar=c(20,2,1,1))

#plot(sort(table(meta$AnnotationFixed)), las=2)

#plot(sort(table(meta$ClusterFix)))

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
names(annot_cols$Date)=c("Old","New")
names(annot_cols$ClusterFix)=unique(meta$ClusterFix)

#flip around so that singletons are first
annot_cols$ClusterFix[c(1,2)] = annot_cols$ClusterFix[c(2,1)]
names(annot_cols$ClusterFix)[c(1,2)]=names(annot_cols$ClusterFix)[c(2,1)]

#make singletons grey
annot_cols$ClusterFix[1]="grey70"


startCols=color.gradient(meta$Start, rainbow(6))

clusCols=color.gradient(1:length(unique(meta$ClusterFix)),cols) [factor(meta$ClusterFix)]

clusCols[meta$ClusterFix=="singleton"]="grey70"


#hist(unlist(ANI))

#ANI[ANI<0.8]=0

# Make an Igraph object from ANI:
network <- graph_from_adjacency_matrix(as.matrix(ANI) , weighted=T, mode="undirected", diag=F )

#Create layout
# l=layout_nicely(network) 
# write.table(l,"good_layout1.txt")



l=as.matrix(read.table("good_layout1.txt"))

#plot(l,col=clusCols, pch=16, cex=0.7)


l2=spreadClus(l, minSize = 10, spreadVal = 0.05, topOutlier = 10, moveVal = 4.5, outlierscale =2.5, useDBscan = F)

#l2=spreadClus(l, minSize = 10, spreadVal = 0.1, topOutlier = 15, moveVal = 0.1, outlierscale = 0.9, plotPoints = T)

par(mfrow=c(1,1))
# plot(l,col=clusCols, pch=16, cex=0.7)
# plot(l2,col=clusCols, pch=16, cex=0.7)
#plot(l2,bg=clusCols,  pch=ifelse(meta$Date=="NCBI", 21,22), cex=sqrt(meta$Length)/200, main=ucPath)
#plot(l, main=ucPath)

plot(l2, col=0, main=ucPath)
text(l2, labels = meta$ClusterFix, col=clusCols)


# plot(l_60, bg="red", xlim=c(-20,90), ylim=c(-50,90), pch=21)
# points(l_70, bg="green", pch=21)
# legend("topleft",legend = c("60","70"), pch=21,pt.bg=c("red", "green"))
# 
# 
# for(k in 1:NROW(bothLayouts)){
#   lines(x = bothLayouts[k,c(1,3)], y = bothLayouts[k,c(2,4)], col="grey90")
# }


#plot(sqrt(meta$Length)/200, col=startCols)

i="45"
for(i in unique(meta$ClusterFix)) {
  
  if(i=="singleton") {
    next
  }
  clusterMeta=meta[which(meta$ClusterFix==i),]
  if(NROW(clusterMeta)<10) {
    next
  }
  
  clusterCoords=l2[which(meta$ClusterFix==i),]
  clusterMeans=colMeans(clusterCoords)  
  text(clusterMeans[1],clusterMeans[2], unique(clusterMeta$AnnotationFixed))
  
}

# 
# text(l2,bg=startCols,  labels = meta$Date)


#plot(seq(0,1E6, 1E4),y=rep(1,length(seq(0,1E6, 1E4))), col=color.gradient(seq(0,1E6, 1E4), rainbow(6)))
#legend("bottomleft",legend = c("0","","","","","","5*10^6" ), pch=15,col= color.gradient(seq(1,6,1), rainbow(6)),y.intersp = .4 )

#table(ifelse(meta$Date=="Old", "square","circle"))

#table(meta$Date)

if(printPdf) pdf("Figs/PhageNetwork_Position70.pdf", width=16,height = 16)

allPlot=plot(x = network,  
             axes=F,rescale=T,layout=l2, 
             # === vertex
             vertex.color = startCols,# clusCols,#c("black","red")[factor(meta$Date)],   # Node color
             #vertex.frame.color = ifelse(meta$Date=="Old", "grey70","black"),# Node border color
             vertex.shape=ifelse(meta$Date=="NCBI", "square","circle"),                        # One of “none”, “circle”, “square”, “csquare”, “rectangle” “crectangle”, “vrectangle”, “pie”, “raster”, or “sphere”
             vertex.size=sqrt(meta$Length)/50,                               # Size of the node (default is 15)
             vertex.size2=NA,                               # The second size of the node (e.g. for a rectangle)
             
             # === vertex label
             vertex.label="",#translationTable2$nodeNames,                 # Character vector used to label the nodes
             vertex.label.color="black",
             vertex.label.font=2,                          # Font: 1 plain, 2 bold, 3, italic, 4 bold italic, 5 symbol
             vertex.label.cex=.8,                         # Font size (multiplication factor, device-dependent)
             vertex.label.dist=0,                          # Distance between the label and the vertex
             vertex.label.degree=0 ,                      # The position of the label in relation to the vertex (use pi)
             
             edge.color=color.gradient(log(E(network)$weight+1),colors = c("grey90","yellow","red")),#"grey50",                           # Edge color
             edge.width=E(network)$weight,#edge.betweenness(network)*0.01,                                 # Edge width, defaults to 1
             edge.arrow.size=1,                            # Arrow size, defaults to 1
             edge.arrow.width=1,                           # Arrow width, defaults to 1
             edge.lty="solid",                             # Line type, could be 0 or “blank”, 1 or “solid”, 2 or “dashed”, 3 or “dotted”, 4 or “dotdash”, 5 or “longdash”, 6 or “twodash”
             edge.curved=0.3      ,                        # Edge curvature, range 0-1 (FALSE sets it to 0, TRUE to 0.5)
             #axes=T
             #xlim=c(min(l),max(l)),ylim=c(min(l),max(l))#, asp = 0
)
# plot(l, col= factor(meta$Date))
legend("topleft",legend = c( "NCBI","River bank") , col = 1, pch=c(15,16), cex = 1.5, title = "Source", bty = "n",)

legend("bottomleft",legend = c("","","0","","","","","","5*10^6" ), pch=15,title="Insertion site", bty="n",
       col= c(0,0,color.gradient(seq(1,6,1), rainbow(6))),y.intersp = .4,cex=1.5 )


legend("bottomright",legend = c("","",min(meta$Length),"","","", max(meta$Length)),cex=1.5, pch=16, col="grey70", 
       pt.cex = c(0,0,seq(0.5, 1.5, length.out=5)), bty="n", title= "Length",y.intersp = .4 )
#legend("bottomright",legend = c("phage", "plasmid"), col = 1, pch=c(16,15), cex=1.5)

if(addAnnotation) {
  for(i in unique(meta$ClusterFix)) {
    
    if(i=="singleton") {
      next
    }
    clusterMeta=meta[which(meta$ClusterFix==i),]
    if(NROW(clusterMeta)<10) {
      next
    }
    
    clusterCoords=norm_coords(l2)[which(meta$ClusterFix==i),]
    clusterMeans=colMeans(clusterCoords)
    text(clusterMeans[1],clusterMeans[2]-.2, unique(clusterMeta$AnnotationFixed))
    
  }
}

if(printPdf) dev.off()


################


if(printPdf) pdf("Figs/PhageNetwork_Clusters70.pdf", width=16,height = 16)

allPlot=plot(x = network,  
             axes=F,rescale=T,layout=l2, 
             # === vertex
             vertex.color = clusCols,# clusCols,#c("black","red")[factor(meta$Date)],   # Node color
             #vertex.frame.color = ifelse(meta$Date=="Old", "grey70","black"),# Node border color
             vertex.shape=ifelse(meta$Date=="NCBI", "square","circle"),                        # One of “none”, “circle”, “square”, “csquare”, “rectangle” “crectangle”, “vrectangle”, “pie”, “raster”, or “sphere”
             vertex.size=sqrt(meta$Length)/50,                               # Size of the node (default is 15)
             vertex.size2=NA,                               # The second size of the node (e.g. for a rectangle)
             
             # === vertex label
             vertex.label="",#translationTable2$nodeNames,                 # Character vector used to label the nodes
             vertex.label.color="black",
             vertex.label.font=2,                          # Font: 1 plain, 2 bold, 3, italic, 4 bold italic, 5 symbol
             vertex.label.cex=.8,                         # Font size (multiplication factor, device-dependent)
             vertex.label.dist=0,                          # Distance between the label and the vertex
             vertex.label.degree=0 ,                      # The position of the label in relation to the vertex (use pi)
             
             edge.color=color.gradient(log(E(network)$weight+1),colors = c("grey90","yellow","red")),#"grey50",                           # Edge color
             edge.width=E(network)$weight,#edge.betweenness(network)*0.01,                                 # Edge width, defaults to 1
             edge.arrow.size=1,                            # Arrow size, defaults to 1
             edge.arrow.width=1,                           # Arrow width, defaults to 1
             edge.lty="solid",                             # Line type, could be 0 or “blank”, 1 or “solid”, 2 or “dashed”, 3 or “dotted”, 4 or “dotdash”, 5 or “longdash”, 6 or “twodash”
             edge.curved=0.3      ,                        # Edge curvature, range 0-1 (FALSE sets it to 0, TRUE to 0.5)
             #axes=T
             #xlim=c(min(l),max(l)),ylim=c(min(l),max(l))#, asp = 0
)
# plot(l, col= factor(meta$Date))
# text(l, labels = meta$Date)
legend("topleft",legend = c( "NCBI","River bank") , col = 1, pch=c(15,16), cex = 1.5, title = "Source", bty = "n",)
#legend("bottomleft",legend = c("0","","","","","","5*10^6" ), pch=15,title="Insertion site", bty="n", col= color.gradient(seq(1,6,1), rainbow(6)),y.intersp = 1,cex=1.5 )

legend("bottomright",legend = c("","",min(meta$Length),"","","", max(meta$Length)),cex=1.5, pch=16, col="grey70", 
       pt.cex = c(0,0,seq(0.5, 1.5, length.out=5)), bty="n", title= "Length",y.intersp = .4 )

legCols=c("grey70",color.gradient(1:24,cols))
legend("bottomleft",legend = c("","","","", rep("",25)),col = c(0,0,0,0,legCols), pch=15,title.adj = 5,
       cex = 1.5,y.intersp = .2, tit, title = "Cluster ID", bty = "n")


if(printPdf) dev.off()

