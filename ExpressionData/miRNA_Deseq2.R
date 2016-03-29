#Now we don't need to specity the working directory
#setwd("/media/somics/Data/gitrepo/mydiary")
library("DESeq2")
library("BiocParallel")
#library("edgeR")
register(MulticoreParam(4))

#filenames <- file.path(dir, paste0(sampleTable$Run, "_subset.bam"))
#Let's take breed as group and diet as condition. Now, we want to find the differences between breed groups
#influenced by diet conditions. Therefore, our aim is to see how flushing diet has effect within breed
#ie. FINF vs FINC or between breed i.e FINF vs TEXF and FINC vs TEXC and so on.

#Features (exons) for each sample were counted using HTSeq-count. We will use those counts for the analysis.
#directory<-"/media/ejo129/hd2/mRNA/DifferentialExpression/DESeq_htseq_Modified_using_only_37C6/counts"
#sampleFiles<-list.files(directory)
#sampleFiles
micountdata<-read.csv("ExpressionData/novel_known_miRNA.csv",row.names = 1, header=TRUE)
head(micountdata)
#The sample table was constructed using sample names with respective groups and conditions.
misampleinfo<-read.csv("ExpressionData/SampleTable_miRNA2S.csv", header=TRUE)
head(misampleinfo)
misample<-misampleinfo$Sample
mibreed<-misampleinfo$Breed
midiet<-misampleinfo$Diet
misampletable<-data.frame(sample=misample, breed=mibreed, diet=midiet)
migroup<-factor(paste(mibreed, midiet, sep="."))
misampletable<-cbind(misampletable, migroup=migroup)
misampletable
colnames(micountdata)<-misampletable$sample
middsMat<-DESeqDataSetFromMatrix(countData=micountdata, colData=misampletable, design=~migroup)
midds<-DESeq(middsMat)
miexpre_all<-results(midds)
miexpre_bm5<-subset(midds,miexpre_all$baseMean >=5 )
miexpre_bm50<-subset(midds, miexpre_all$baseMean >=50)
miexpre_bm100<-subset(midds, miexpre_all$baseMean >= 100)
nrow(miexpre_bm5)#391 (281 novel) 
nrow(miexpre_bm50)# 236 (153 novel) 
nrow(miexpre_bm100)#191 (123 novel)
miexpre_bm100_data<-mcols(miexpre_bm100, use.names=TRUE)
#write.csv(expre_bm100_data, "novelmiRNAsBaseMeangeq100.csv")
miexpre_bm10<-subset(midds, miexpre_all$baseMean>=10)
nrow(miexpre_bm10)# 342 (238 novel)
miexpre_bm10_data<-mcols(miexpre_bm10, use.names=TRUE)
write.csv(miexpre_bm10_data, "ExpressionData/miRNA_out/known_novelmiRNAsBaseMeangeq10.csv")

# miknFS_TX<-results(midds, contrast=c("breed","FS","TX"))
# miknFS_TX_sig<-subset(miknFS_TX, padj < 0.1)
# nrow(miknFS_TX_sig)#38 (18 novel)
# write.csv(miknFS_TX_sig, "mikn_FS_TX_sig.csv")
# miknFS_TX_sig
# 
# miknFS_F1<-results(midds, contrast=c("breed", "FS", "F1"))
# miknFS_F1_sig<-subset(miknFS_F1, padj < 0.1)
# nrow(miknFS_F1_sig)#35(15 novel)
# write.csv(miknFS_F1_sig, "mikn_FS_F1_sig.csv")
# miknFS_F1_sig
# 
# miknTX_F1<-results(midds, contrast=c("breed", "TX", "F1"))
# miknTX_F1_sig<-subset(miknTX_F1, padj< 0.1)
# nrow(miknTX_F1_sig)#5(3 novel)
# write.csv(miknTX_F1_sig, "mikn_TX_F1_sig.csv")
# 
# miknF_C<-results(midds, contrast=c("diet", "F","C"))
# miknF_C_sig<-subset(miknF_C, padj<0.1)
# nrow(miknF_C_sig)#0
# miknF_C_sig
#write.csv(knF_C_sig, "kn_F_C_sig.csv")


# group<-factor(paste(breed, diet, sep="."))
# sampletable1<-cbind(sampletable, group=group)
# sampletable1
# ddsMat1<-DESeqDataSetFromMatrix(countData=countdata, colData=sampletable1, design=~group)
# dds1<-DESeq(ddsMat1)
# dds1
# res<-results(dds1)
# res
# rld<-rlog(dds1)
# plotPCA(rld, intgroup=c("breed", "diet"))
# plotDispEsts(dds1)
# plotMA(dds1)

midknFS_Diet<-results(midds, contrast=c("migroup", "FS.F", "FS.C"))
midknFS_Diet_sig<-subset(midknFS_Diet, padj < 0.05)
nrow(midknFS_Diet_sig)#0

midknTX_Diet<-results(midds, contrast=c("migroup", "TX.F", "TX.C"))
midknTX_Diet_sig<-subset(midknTX_Diet, padj< 0.05)
nrow(midknTX_Diet_sig)#0



midknF1_Diet<-results(midds, contrast=c("migroup", "F1.F", "F1.C"))
midknF1_Diet_sig<-subset(midknF1_Diet, padj<0.05)
nrow(midknF1_Diet_sig)#0
#write.csv(midknF1_Diet_sig, "miRNA_out/midknF1_Diet_sig.csv")


midknFS_TX<-results(midds, contrast=list(c("migroupFS.F", "migroupTX.C"), c("migroupFS.C", "migroupTX.F")))
midknFS_TX_sig<-subset(midknFS_TX, padj<0.05)
nrow(midknFS_TX_sig)#0



midknFS_F1<-results(midds, contrast=list(c("migroupFS.F", "migroupF1.C"), c("migroupFS.C", "migroupF1.F")))
midknFS_F1_sig<-subset(midknFS_F1, padj<0.05)
nrow(midknFS_F1_sig) #0



midknTX_F1<-results(midds, contrast=list(c("migroupTX.F", "migroupF1.C"), c("migroupTX.C", "migroupF1.F")))
midknTX_F1_sig<-subset(midknTX_F1, padj<0.05)
nrow(midknTX_F1_sig)#0


midknfFS_TX<-results(midds, contrast=c("migroup", "FS.F", "TX.F"))
midknfFS_TX_sig<-subset(midknfFS_TX, padj<0.05)
nrow(midknfFS_TX_sig) #0
#summary(midknfFS_TX_sig)
#write.csv(midknfFS_TX_sig, "ExpressionData/miRNA_out/midknfFS_TX_sig.csv")

midkncFS_TX<-results(midds, contrast=c("migroup", "FS.C", "TX.C"))
midkncFS_TX_sig<-subset(midkncFS_TX, padj < 0.05)
nrow(midkncFS_TX_sig)#3
summary(midkncFS_TX_sig)
write.csv(midkncFS_TX_sig, "ExpressionData/miRNA_out/midkncFS_TX_sig.csv")

midknfFS_F1<-results(midds, contrast=c("migroup", "FS.F", "F1.F"))
midknfFS_F1_sig<-subset(midknfFS_F1, padj<0.05)
nrow(midknfFS_F1_sig) #0
#summary(midknfFS_F1_sig)
#write.csv(midknfFS_F1_sig, "ExpressionData/miRNA_out/midknfFS_F1_sig.csv")

midncFS_F1<-results(midds, contrast=c("migroup", "FS.C", "F1.C"))
midncFS_F1_sig<-subset(midncFS_F1, padj < 0.05)
nrow(midncFS_F1_sig) #0
summary(midncFS_F1_sig)

midknfTX_F1<-results(midds, contrast=c("migroup", "TX.F", "F1.F"))
midknfTX_F1_sig<-subset(midknfTX_F1, padj<0.05)
nrow(midknfTX_F1_sig)#0
summary(midknfTX_F1_sig)
#write.csv(midknfTX_F1_sig, "ExpressionData/miRNA_out/midknfTX_F1_sig.csv")

midkncTX_F1<-results(midds, contrast=c("migroup", "TX.C", "F1.C"))
midkncTX_F1_sig<-subset(midkncTX_F1, padj < 0.05)
nrow(midkncTX_F1_sig) #0

