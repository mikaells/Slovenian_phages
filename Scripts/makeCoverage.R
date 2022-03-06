source("Scripts/functions.R")
aa=AnnotateClusters(blastDir = "Data/blast_out_genomes70/", headerPath = "Data/genome_headers.txt", plotAlligns = F)

which.max(maxes)

pdf(file = "coveragePerCluster.pdf")
par(mar=c(2.5,2.5,0.5,0.5), mgp=c(1.4,0.4,0), font.lab=2)

plot(sort(aa$maxCov), pch=16, cex=.5, ylab="Best overlap [%]")
abline(v=53.5, col="grey")
abline(h=seq(0,100, 10), col="gray")
dev.off()

hist(maxes, breaks = 100)

