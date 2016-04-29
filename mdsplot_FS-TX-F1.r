#Lets make MDS plot of three breeds based on their genotypes. 
#The file comes from /media/ejo129/somics/Data/SNP_genotype_data/Release_020414/Oct15/FS-TX-F1/
#the plot was made using qqman package https://cran.r-project.org/web/packages/qqman/vignettes/qqman.html 
library(ggplot2)
library(qqman)
require(graphics)
gwasResults<-read.table("IBS-MDS/FS-TX-F1_assoc.qassoc", header=T)
gwasResults<-subset(gwasResults, P>0)
manhattan(gwasResults, genomewideline=F, suggestiveline=F, ylim=c(0,10), cex=0.25, cex.axis=0.8, col=palette())
#for specific chromosome
#manhattan(subset(gwasResults, CHR==2), genomewideline=F, suggestiveline=F, ylim=c(0,10), main="Chr 2",cex=0.25 )
qq(gwasResults$P)


mds<-read.table("IBS-MDS/MDS_FS-TX-F1_mibs.mds", header=T)
head(mds)
ids<-read.table("IBS-MDS/familyid.csv", header=T)
ids
colnames(ids)<-c("FID", "FS", "TX", "F1", "Breed")
Breed<-ids$Breed
FID<-ids$FID
FID
mds1<-data.frame(mds, Breed)
qplot(C1, C2, data=mds1, color=Breed, alpha=I(0.5), size=I(3))
#qplot(C1, C2, data=mds1, color=Breed, alpha=I(0.5), size=I(3), label=FID, geom = "text") #to get labels for all samples


mds_king<-read.table("/media/ejo129/somics/Data/SNP_genotype_data/Release_020414/Oct15/FS-TX-F1/kingpc1.csv", header=T)
head(mds_king)
mds_king1<-data.frame(mds_king, Breed)
qplot(C1, C2, data=mds_king1, color = Breed, alpha=I(0.5), size=I(3))

FS<-1:11
TX<-12:22
F1<-23:31



dev.off()


