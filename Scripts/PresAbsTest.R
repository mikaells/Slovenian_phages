library(pheatmap)
library(stringr)
library(RColorBrewer)
library(igraph)

rm(list=ls())


source("Scripts/functions.R")

printPdf=F
aniPath="Data/ani/ANIm_hadamard.tab"
ucPath="Data/all.uc70"
minCluster=3


g(ANI, meta) %=% BuildData(aniPath = aniPath, ucPath = ucPath,minCluster = minCluster)#, subs=1:100)

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

plot(meta$Length~meta$Start, col=annot_cols$ClusterFix[meta$ClusterFix], pch=as.numeric(factor(meta$ClusterFix)))

findEpsi(as.matrix(meta$Start),maxY = 200, maxRange = 50000, steps=100, minP = 20)
abline(v=31500)

DBC=dbscan::dbscan(as.matrix(meta$Start),eps = 31500, minPts = 20)


plot(meta$Length~meta$Start, col=DBC$cluster+1, corral="wrap")

plot(meta$Length~meta$Start)

annot_cols$ClusterFix[meta$ClusterFix]

df=makeIsoDF(meta)




chisq.test(df[,2], df[,3])

which("25"== colnames(df))

which("31"== colnames(df))

plot(jitter(makeBin(df[,10])), jitter(makeBin(df[,16])), xlab="25", ylab="31")

i=26;j=1
MIN=1
MIN_indx=c("")
par(mfrow=c(4,4), mar=c(2.5,2.5,0,0), mgp=c(1.4,0.4,0))
for(i in 2:NCOL(df)) {
  for(j in 2:NCOL(df)) {
        
    if(i==j) next
    
    eitherPres=df[df[,i]>0 | df[,j]>0,c(i,j)]
    
  
    TAB=table(data.frame(makeBin(eitherPres)))
    
    CHI=chisq.test(TAB)
    
    if(CHI$p.value<MIN) {
      cat(paste("Previous min:", MIN,"\n"))
      MIN_indx=paste(i,j)
      MIN=CHI$p.value
      cat(paste("New min", MIN, MIN_indx,"\n"))
      
    }
    
    if(CHI$p.value<1E-30) {
      pheatmap(data.frame(makeBin(eitherPres)), cluster_rows = T, cluster_cols = T)
      
      #plot(jitter(eitherPres[,1])~ jitter(eitherPres[,2]))
      #cat(paste(i, j,"|",CHI$p.value),"\n")
    }
  }
}
    #GLM=glm(makeBin(df[,i]) ~ makeBin(df[,j]), family="binomial")
    
    #AOV=anova(GLM, test="Chisq")
    
    TAB=table(data.frame(makeBin(df[,c(i,j)])))
    
    CHI=chisq.test(TAB)
    
    #CHI=chisq.test(makeBin(df[,i]), makeBin(df[,j]))
    if(CHI$p.value<0.01) {
      
      cat(paste(i, j,"|",CHI$p.value),"\n")
      
      
      plot(jitter(makeBin(df[,i])), jitter(makeBin(df[,j])), main=paste(i,j),xlab=i, ylab=j)
      mosaicplot(TAB)
    }
  }
}

chisq.test(makeBin(df[,10]), makeBin(df[,16]))


df2=df[,c(1,1+order(colSums(df[,-1])))]
df3=df2[order(rowSums(df2[,-1])),]

(colSums(df3[,-1]))

(rowSums(df3[,-1]))


pheatmap(df3[,-1], cluster_cols = F, cluster_rows = F)
