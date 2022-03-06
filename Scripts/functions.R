makeBin=function(x) {
  
  return(ifelse(x==0,0,1))
  
}

color.gradient <- function(x, colors=c("red","yellow","green"), colsteps=100) {
  return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )
}



findEpsi=function(L, minRange=0,maxRange=3, steps=0.1, maxY=200, minP=10) {
  plot(0,0, col=0,xlim=c(minRange,maxRange), ylim=c(0,maxY))
  legend("topright", legend = c("nClust","nOutliers"), col=1:2, pch=16)
  for(i in seq(minRange,maxRange,steps)) {
    
    DB=dbscan::dbscan(L,i,minP)
    
    points(i,length(unique(DB$cluster)), pch=16, col=1)
    points(i,length(which(DB$cluster==0)), pch=16,col=2)
    
  }
}

makeIsoDF=function(meta) {
  
  isolates=naturalsort::naturalsort(unique(meta$ID))
  uniqClus=naturalsort::naturalsort(unique(meta$ClusterFix))
  
  df <- data.frame(matrix(0,ncol =1+length(uniqClus) , nrow = length(isolates)))
  colnames(df) <- c("Isolates", uniqClus)
  
  df$Isolates=isolates
  
  rownames(df)= isolates
  
  i=1
  for(i in 1:NROW(isolates)) {
    isoName=df$Isolates[i]
    isoIndx=which(isoName == meta$ID)
    iso=meta[isoIndx,]
    
    j=1
    for(j in 1:NROW(iso)) {
      df[i,which(iso$ClusterFix[j] == colnames(df))]=df[i,which(iso$ClusterFix[j] == colnames(df))] +1
    }
  }
  
  rownames(df)=df$Isolates
  df=df[,-1]
  
  return(df)
}

#l: layout to mess with
#transTab: basically metadata, dont actually use it for anything
#eps: epsilon value for dbscan. Use findEpsi() to work that out
#minSize: minimum size of a cluster before you spread it
#spreadVal: how much to spread a cluster
#topoutlier: outliers are calculated as the N most distant points from the center
#how much to move the outliers, .75 will move them to 75% of original distance
#moveVal: how much to move clusters away from center
#plotPoints: to plot or not, useful for finding values


spreadClus=function(L,transTab, eps=0.5,  minSize=16, spreadVal=0.1, topOutlier=20, outlierscale=0.75,moveVal=1.1, plotPoints=F, useDBscan=T) {
  
  
  if(useDBscan) {
    #find major clusters
    dbclus=dbscan::dbscan(L, eps =eps)
    
    #save them with main table
    clusterID=dbclus$cluster+1
  } else {
    clusterID=meta$Cluster+1
  }
  
  #have a look
  if(plotPoints) plot(L, col=0)#translationTable2$dbcan+1)
  
  #make a new layout to overwrite
  l2=L
  
  if(plotPoints) plot(l2, col=clusterID, pch =16, cex =.1)
  
  #plot again
  if(plotPoints) plot(L, col=1,pch =16, cex =.1)
  
  j=1
  
  
  allDist=sqrt(l2[,1]^2+l2[,2]^2)
  
  bigDist=order(allDist, decreasing = T)[1:topOutlier]
  
  l2[bigDist,]=l2[bigDist,]*outlierscale
  
  
  #spreading out clusters
  #for all clusters in graph, do
  for(j in unique(clusterID)) {
    if(j==0) next #if cluster is just the non-clustered
    #find index of cluster members and calculate centroid
    dbIndx=which(clusterID==j)
    dbCent=c(mean(l[dbIndx,1]),mean(l[dbIndx,2]))
    
    #if cluster is too small, next
    if(length(dbIndx)<minSize) next
    
    #for all cluster members, do
    i=1
    for(i in (dbIndx)){
      #find relative position and distance to center
      xx=(l[i,][1]-dbCent[1])
      yy=(l[i,][2]-dbCent[2])  
      dd=sqrt(xx^2+yy^2)
      
      #to each point, add max three but less if cluster is small 
      #to pull away from center
      #l2[i,]=l[i,]+(3-(3/(1*(length(dbIndx))^(1/10))))*c(xx,yy)
      l2[i,]=l[i,]+(spreadVal*sqrt(length(dbIndx)))*c(xx,yy)
      #points(x = l2[i,1],y=l2[i,2], pch=16,cex=0.1,col=meta$dbcan[i]+1)
    }
    
    allDist=sqrt(l2[,1]^2+l2[,2]^2)
    
    l2[dbIndx,]=l2[dbIndx,]*moveVal
    
  }
  
  
  
  if(plotPoints) points(l2, col=clusterID, pch =16, cex =.1)
  
  
  return(l2)
  
}

#aniPath = aniPath; ucPath = "Data/all_phages.uc60_def2";minCluster = minCluster; meshPath=meshPath

BuildData=function(aniPath,ucPath, minCluster,subs, meshPath, useMesh=F ) {
  
  #Read identity matrix
  ANI=read.table(aniPath)
  
  if(!missing(subs)) {
    ANI=ANI[subs,subs]
  }
  
  #Read vsearch clustering data
  uc=read.table(ucPath, sep="\t")
  colnames(uc) = c("Record","Cluster","Size","ID","Strand","U1","U2","Aln","Q","T")
  
  #read meshcluster data
  mesh=read.table(meshPath)
  mesh$V2=gsub(">","",mesh$V2)
  
  #store raw names before cleaning them up
  rawNames0=rownames(ANI)
  rawNames=rawNames0
  
  #manually fix genomes with weird names
  if(missing(subs)) {
    rawNames[335] = "25_phage_plasmid_1a__length_86702_de..."
    rawNames[666] = "25_phage_plasmid_1b_length_86702_de..."
    rawNames[818] = "24_phage_plasmid_2a_length_14181_de..."
    rawNames[837] = "24_phage_plasmid_2b_length_14181_de..."
  }
  
  #Fix odd names to make them splittable into meta-data
  fixedNames=gsub("__","_",gsub("phage_plasmid_(.)", "plasmid_\\1",gsub(pattern = ":","_phages_1_", rawNames)))
  
  #make meta-data by splitting on "_"
  meta=data.frame(str_split_fixed(string = fixedNames, pattern = "_|:", 6), stringsAsFactors = F)
  
  #make empty variables for length, starting positions and vsearch cluster
  meta$len=0
  meta$start=0
  meta$Cluster=-1
  meta$meshClust=-1
  meta$ClusterFix=-1
  
  #match rownames for meta as the raw names to match back to vsearch cluster
  rownames(meta)=rawNames0
  
  table(rawNames0 %in% uc$Q)
  
  #i=1
  #run through meta-file
  for(i in 1:NROW(meta)) {
    
    #if genome name implies it is a full plasmid
    #length is in the name and there is no start
    if(meta[i,4]=="length") {
      meta$len[i] = as.numeric(meta[i,5])
      meta$start[i]=1
      #else calculate length from region
      #and likewise for start
    } else {
      temp=as.numeric(str_split_fixed(meta[i,4],"-",2)  )
      meta$len[i]=temp[2]-temp[1]
      meta$start[i]=temp[1]
    }
    
    #match genome name with the vsearch cluster data
    UC_temp=uc[grep(rawNames0[i], uc$Q),]  
    
    #fish the cluster out and add to meta-data
    #some genomes are both centroids and hits, so unique handles this
    meta$Cluster[i]=unique(UC_temp$Cluster)
    
    
    meshTemp=mesh[grep(rawNames0[i], mesh$V2),]  
    meta$meshClust[i]=unique(meshTemp$V1)
    
  }
  
  if(useMesh) {
    
    #work out which clusters are very small
    singleClust=names(which(table(meta$meshClust)<minCluster))
    
    #overwrite small clusters with just 'singleton'
    meta$ClusterFix=ifelse(meta$meshClust %in% singleClust,"singleton",meta$meshClust)
    
  } else {
    #work out which clusters are very small
    singleClust=names(which(table(meta$Cluster)<minCluster))
    
    #overwrite small clusters with just 'singleton'
    meta$ClusterFix=ifelse(meta$Cluster %in% singleClust,"singleton",meta$Cluster)
    
  }
  
  #clean up metadata  
  meta=meta[,-c(5,6)]
  
  #add nice names to the meta-data
  colnames(meta)=c("ID","Type","Index","Region","Length","Start","Cluster","meshClust", "ClusterFix")
  
  #fix data types and add log lengths and whether or not genomes are new or old
  meta$Type=factor(meta$Type)
  meta$logLen=log(meta$Length)
  meta$Date=ifelse(grepl("CP",meta$ID), "NCBI","River bank")
  
  return(list(ANI, meta))
}


# Generic form
'%=%' = function(l, r, ...) UseMethod('%=%')

# Binary Operator
'%=%.lbunch' = function(l, r, ...) {
  Envir = as.environment(-1)
  
  if (length(r) > length(l))
    warning("RHS has more args than LHS. Only first", length(l), "used.")
  
  if (length(l) > length(r))  {
    warning("LHS has more args than RHS. RHS will be repeated.")
    r <- extendToMatch(r, l)
  }
  
  for (II in 1:length(l)) {
    do.call('<-', list(l[[II]], r[[II]]), envir=Envir)
  }
}

# Used if LHS is larger than RHS
extendToMatch <- function(source, destin) {
  s <- length(source)
  d <- length(destin)
  
  # Assume that destin is a length when it is a single number and source is not
  if(d==1 && s>1 && !is.null(as.numeric(destin)))
    d <- destin
  
  dif <- d - s
  if (dif > 0) {
    source <- rep(source, ceiling(d/s))[1:d]
  }
  return (source)
}

# Grouping the left hand side
g = function(...) {
  List = as.list(substitute(list(...)))[-1L]
  class(List) = 'lbunch'
  return(List)
}

AnnotateClusters=function(blastDir="Data/blast_out_genomes70/",headerPath="Data/genome_headers.txt", verbose=F,plotAlligns=F ) {
  library(stringr)
  headers=readLines(headerPath)
  headerDF=cbind(headers,str_split_fixed(headers, " " , 2))
  headerDF[,2]=gsub(">","", headerDF[,2])
  
  blastCols=c( "qseqid", "sseqid", "pident", "length", "qcovs", "qcovhsp",  "mismatch", 
               "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
  
  blastFiles=naturalsort::naturalsort(dir(blastDir, full.names = T, pattern = ".blast$"))
  
  blast=blastFiles[1]
  
  summaryDF=data.frame(clusterFile=basename(blastFiles),cluster=gsub(".*uc.._([0-9]*).blast$","\\1", blastFiles),bestVote="None",voteFrac=0, maxCov=0)
  
  blastCounter=1
  blast=blastFiles[1]
  for(blast in blastFiles){
    
    # if(blast=="blast_out_genomes/all.uc70_21.blast"){
    #   break
    # }
    
    if(file.info(blast)$size==0){
      if(verbose) {
        cat(paste(basename(blast), "has no hits.\n\n______\n\n"))
      }
      next
    }
    
    blastFile=read.table(blast)
    colnames(blastFile)=blastCols
    
    if(plotAlligns) {
      plot(0, 0, xlim=c(0,max(blastFile$qstart)),ylim=c(0,length(unique(blastFile$sseqid))), col=0,main=basename(blast))
      
      m=unique(blastFile$sseqid)[1]
      subCounter=1
      for(m in unique(blastFile$sseqid)) {
        idCols=color.gradient(blastSub$pident)
        blastSub=subset(blastFile, sseqid == m)
        
        for(k in 1:NROW(blastSub)){
          lines(x=c(blastSub$qstart[k],blastSub$qend[k]), y=c(subCounter,subCounter), col=idCols[k], lwd=2 )
        }
        subCounter=subCounter+1
      }
    }
    
    blastFile$weightId=(blastFile$pident/100) * blastFile$qcovhsp
    
    uniqGenos=unique(blastFile$qseqid) 
    
    bestHits=data.frame(uniqGenos,best="")
    j_counter=1
    for(j in uniqGenos ) {
      uniqGenoIndx=which(blastFile$qseqid==j)
      
      uniqGenoBlast=blastFile[uniqGenoIndx,]
      
      uniqMatches=unique(uniqGenoBlast$sseqid)
      
      sums=data.frame(uniqMatches,sums=0)
      k_counter=1
      for(k in uniqMatches) {
        
        uniqMatchIndx=which(uniqGenoBlast$sseqid==k)
        
        sums$sums[k_counter]=sum(uniqGenoBlast[uniqMatchIndx,"weightId"])
        k_counter=k_counter+1
      }
      
      topHit0=sums$uniqMatches[which.max(sums$sums)]
      #topHit=grep(topHit0, headers, value = T, fixed = T)
      topHit = headers[which(topHit0==headerDF[,2])]
      
      bestHits$best[j_counter]=topHit
      
      j_counter=j_counter+1
    }
    
    bestTab=sort(table(bestHits$best), decreasing = T)
    
    summaryDF$bestVote[blastCounter]=names(bestTab)[1]
    summaryDF$voteFrac[blastCounter]=(bestTab/sum(bestTab))[1]
    summaryDF$maxCov[blastCounter]=max(blastFile$qcovs)
    
    if(verbose) {
      cat(paste(basename(blast),"| max =", max(blastFile$qcovs),"%\n"))
      print(table(bestHits$best))
      cat(paste("________\n\n"))
    }
    
    blastCounter=blastCounter+1
  }
  
  return(summaryDF)
}
