#Differential expression analysis on within  and between breeds that are affected by diet treatment. We have three breed groups: Finnsheep (high prolific), Texel (low prolific)
# and F1-cross of Finnsheep and Texel. In all breed groups, half of the ewes were kept in flushing diet (additional nutrient). Thus, here we find how diet is affecting within the breed but also how 
#different breeds react differently to diet (i.e genes that are differentially expressed with diet as a contributing factor).
#setwd("/media/somics/Data/mRNA/2/Aug15/DiffExp/39Samples_Jan16/")
library("biomaRt")
library("DESeq2")
library("BiocParallel")
library("RColorBrewer")
library("gplots")
register(MulticoreParam(6))
#listMarts(host='www.ensembl.org')
ensembl = useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", dataset="oaries_gene_ensembl")
#ensembl = useDataset("oaries_gene_ensembl", mart = ensembl)
#ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="oaries_gene_ensembl", host = "jul2015.archive.ensembl.org")
#Features (exons) for each sample were counted using HTSeq-count. We will use those counts for the analysis.

##################################################################################################
##################################################################################################
directory<-"ExpressionData/mRNA_counts/"
sampleFiles<-grep("counts", list.files(directory), value = TRUE)
sampleFiles<-list.files(directory, pattern = "*_counts")
sampleFiles
sampleinfo<-read.csv("ExpressionData/SampleTable_mRNA_all.csv", header=TRUE)
row.names(sampleinfo)
sample<-sampleinfo$Sample
breed<-sampleinfo$Breed
diet<-sampleinfo$Diet
id<-sampleinfo$SampleID
sampletable<-data.frame(sampleName=sample, fileName=sampleFiles,breed=breed, diet=diet, id=id)
group<-factor(paste(breed, diet, sep="."))
sampletable<-cbind(sampletable, group=group)
sampletable

#****************************************************************************************
#Comparisions including diet as second factor
#****************************************************************************************
ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable = sampletable, directory = directory, design = ~group)
ddsHTSeqcol<-collapseReplicates(ddsHTSeq, sampletable$id)
dds<-DESeq(ddsHTSeqcol)
dds
res<-results(dds)
res_baseMean5 <- subset(dds, res$baseMean > 5)
nrow(res_baseMean5)#16 402
#The effect of diet in Finnsheep. i.e Finnsheep with flushing diet vs Finnsheep with control diet.
dFS_Diet<-results(dds, contrast=c("group", "FS.F", "FS.C"))
dFS_Diet_sig<-subset(dFS_Diet, padj < 0.05)
nrow(dFS_Diet_sig)#0 i.e no effect of diet in Finnsheep

##########################################################################################
#The effect of diet in Texel
##########################################################################################
dTX_Diet<-results(dds, contrast=c("group", "TX.F", "TX.C"))
dTX_Diet_sig<-subset(dTX_Diet, padj< 0.05)
nrow(dTX_Diet_sig)#117
summary(dTX_Diet_sig)
write.csv(dTX_Diet_sig, "mRNA_out/dTX_Diet_sig.csv")
rndTX_Diet_sig<-rownames(dTX_Diet_sig)
#Retrive GO IDs although GO annotations and KEGG pathway analysis are done in Cytoscape using ClueGO plugin
godTX_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id", "go_id"), values=rndTX_Diet_sig, mart=ensembl)
write.csv(godTX_sig, "mRNA_out/godTX_sig.csv")
#We retrieve some basic information for each differentially expressed genes using biomart.
bmdTX_Diet_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id",
                                                                                        "external_gene_name","chromosome_name", "description"), values=rndTX_Diet_sig, mart=ensembl)
#In case some unknown genes in sheep are annotated in other species, it has also been useful to do the GO annotation and KEGG pathways using cow homologs as it's genome is much more annotated.
ortho_bmdTX_Diet_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id", "external_gene_name",
                                                                                              "btaurus_homolog_ensembl_gene", "sscrofa_homolog_ensembl_gene", "mmusculus_homolog_ensembl_gene", "hsapiens_homolog_ensembl_gene"), values=rndTX_Diet_sig, mart=ensembl)

write.csv(ortho_bmdTX_Diet_sig, "mRNA_out/ortho_bmdTX_diet_sig.csv")
sortbmdTX_Diet_sig<-cbind(bmdTX_Diet_sig[match(rndTX_Diet_sig, bmdTX_Diet_sig$ensembl_gene_id),])
bmdTX_Diet<-cbind(dTX_Diet_sig, bmdTX_Diet_sig)
write.csv(bmdTX_Diet, "mRNA_out/bmdTX_Diet_sig.csv")

###################################################################################
#The effect of diet in F1-cross
###################################################################################
dF1_Diet<-results(dds, contrast=c("group", "F1.F", "F1.C"))
dF1_Diet_sig<-subset(dF1_Diet, padj<0.05)
nrow(dF1_Diet_sig)#26
summary(dF1_Diet_sig)
write.csv(dF1_Diet_sig, "mRNA_out/dF1_Diet_sig.csv")
rndF1_Diet_sig<-rownames(dF1_Diet_sig)
godF1_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id", "go_id"), values=rndF1_Diet_sig, mart=ensembl)
write.csv(godF1_sig, "mRNA_out/godF1_sig.csv")
bmdF1_Diet_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id",
                                                                                        "external_gene_name", "chromosome_name", "description"), values=rndF1_Diet_sig, mart=ensembl)
ortho_bmdF1_Diet_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id", "external_gene_name",
                                                                                              "btaurus_homolog_ensembl_gene", "sscrofa_homolog_ensembl_gene", "mmusculus_homolog_ensembl_gene", "hsapiens_homolog_ensembl_gene"), values=rndF1_Diet_sig, mart=ensembl)
write.csv(ortho_bmdF1_Diet_sig, "mRNA_out/ortho_bmdF1_diet_sig.csv")
sortbmdF1_Diet_sig<-cbind(bmdF1_Diet_sig[match(rndF1_Diet_sig, bmdF1_Diet_sig$ensembl_gene_id),])
bmdF1_Diet<-cbind(dF1_Diet_sig, bmdF1_Diet_sig)
write.csv(bmdF1_Diet, "mRNA_out/bmdF1_Diet_sig.csv")

###################################################################################
#The difference in Finnsheep and Texel due to diet effect
###################################################################################
dFS_TX<-results(dds, contrast=list(c("groupFS.F", "groupTX.C"), c("groupFS.C", "groupTX.F")))
dFS_TX_sig<-subset(dFS_TX, padj<0.05)
nrow(dFS_TX_sig)#35
summary(dFS_TX_sig)
write.csv(dFS_TX_sig, "mRNA_out/dFS_TX_sig.csv")
rndFS_TX_sig<-rownames(dFS_TX_sig)
godFS_TX_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id", "go_id"), values=rndFS_TX_sig, mart=ensembl)
write.csv(godFS_TX_sig, "mRNA_out/godFS_TX_sig.csv")
bmdFS_TX_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id",
                                                                                      "external_gene_name", "chromosome_name", "description"), values=rndFS_TX_sig, mart=ensembl)
ortho_bmdFS_TX_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id", "external_gene_name",
                                                                                            "btaurus_homolog_ensembl_gene", "sscrofa_homolog_ensembl_gene", "mmusculus_homolog_ensembl_gene", "hsapiens_homolog_ensembl_gene"), values=rndFS_TX_sig, mart=ensembl)
write.csv(ortho_bmdFS_TX_sig, "mRNA_out/ortho_bmdFS_TX_sig.csv")
sortbmdFS_TX_sig<-cbind(bmdFS_TX_sig[match(rndFS_TX_sig, bmdFS_TX_sig$ensembl_gene_id),])
bmdFS_TX<-cbind(dFS_TX_sig, bmdFS_TX_sig)
write.csv(bmdFS_TX, "mRNA_out/bmdFS_TX_sig.csv")

###################################################################################
#The difference in Finsheep and F1 due to diet effect
###################################################################################
dFS_F1<-results(dds, contrast=list(c("groupFS.F", "groupF1.C"), c("groupFS.C", "groupF1.F")), test = "wald")
dFS_F1_sig<-subset(dFS_F1, padj<0.05)
nrow(dFS_F1_sig) #5
summary(dFS_F1_sig)
write.csv(dFS_F1_sig, "mRNA_out/dFS_F1_sig.csv")
rndFS_F1_sig<-rownames(dFS_F1_sig)
godFS_F1_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id", "go_id"), values=rndFS_F1_sig, mart=ensembl)
write.csv(godFS_F1_sig, "mRNA_out/godFS_F1_sig.csv")
bmdFS_F1_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id",
                                                                                      "external_gene_name", "chromosome_name", "description"), values=rndFS_F1_sig, mart=ensembl)
ortho_bmdFS_F1_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id", "external_gene_name",
                                                                                            "btaurus_homolog_ensembl_gene", "sscrofa_homolog_ensembl_gene", "mmusculus_homolog_ensembl_gene", "hsapiens_homolog_ensembl_gene"), values=rndFS_F1_sig, mart=ensembl)
write.csv(ortho_bmdFS_F1_sig, "mRNA_out/ortho_bmdFS_F1_sig.csv")
sortbmdFS_TX_sig<-cbind(bmdFS_TX_sig[match(rndFS_TX_sig, bmdFS_TX_sig$ensembl_gene_id),])
sortbmdFS_F1_sig<-cbind(bmdFS_F1_sig[match(rndFS_F1_sig, bmdFS_F1_sig$ensembl_gene_id),])
bmdFS_F1<-cbind(dFS_F1_sig, bmdFS_F1_sig)
write.csv(bmdFS_F1, "mRNA_out/bmdFS_F1_sig.csv")

###################################################################################
#The difference in Texel and F1 due to diet effect
###################################################################################
dTX_F1<-results(dds, contrast=list(c("groupTX.F", "groupF1.C"), c("groupTX.C", "groupF1.F")))
dTX_F1_sig<-subset(dTX_F1, padj<0.05)
nrow(dTX_F1_sig)#73
summary(dTX_F1_sig)
write.csv(dTX_F1_sig, "mRNA_out/dTX_F1_sig.csv")
rndTX_F1_sig<-rownames(dTX_F1_sig)
godTX_F1_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id", "go_id"), values=rndTX_F1_sig, mart=ensembl)
write.csv(godTX_F1_sig, "mRNA_out/godTX_F1_sig.csv")
bmdTX_F1_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id",
                                                                                      "external_gene_name", "chromosome_name", "description"), values=rndTX_F1_sig, mart=ensembl)
ortho_bmdTX_F1_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id", "external_gene_name",
                                                                                            "btaurus_homolog_ensembl_gene", "sscrofa_homolog_ensembl_gene", "mmusculus_homolog_ensembl_gene", "hsapiens_homolog_ensembl_gene"), values=rndTX_F1_sig, mart=ensembl)
write.csv(ortho_bmdTX_F1_sig, "mRNA_out/ortho_bmdTX_F1_sig.csv")
sortbmdTX_F1_sig<-cbind(bmdTX_F1_sig[match(rndTX_F1_sig, bmdTX_F1_sig$ensembl_gene_id),])
bmdTX_F1<-cbind(dTX_F1_sig, bmdTX_F1_sig)
write.csv(bmdTX_F1, "mRNA_out/bmdTX_F1_sig.csv")

#######################################################################################
#The difference between Finnsheep and Texel with Control diet
#######################################################################################
cFS_TX <- results(dds, contrast=c("group", "FS.C", "TX.C"))
cFS_TX_sig<-subset(cFS_TX, padj < 0.05)
nrow(cFS_TX_sig) #3
summary(cFS_TX_sig)
write.csv(cFS_TX_sig, "mRNA_out/cFS_TX_sig.csv")
rncFS_TX_sig<-rownames(cFS_TX_sig)
gocFS_TX_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id", "go_id"), values=rncFS_TX_sig, mart=ensembl)
write.csv(gocFS_TX_sig, "mRNA_out/gocFS_TX_sig.csv")
bmcFS_TX_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id",
                                                                                      "external_gene_name", "chromosome_name", "description"), values=rncFS_TX_sig, mart=ensembl)
ortho_bmcFS_TX_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id", "external_gene_name",
                                                                                            "btaurus_homolog_ensembl_gene", "sscrofa_homolog_ensembl_gene", "mmusculus_homolog_ensembl_gene", "hsapiens_homolog_ensembl_gene"), values=rncFS_TX_sig, mart=ensembl)
write.csv(ortho_bmcFS_TX_sig, "mRNA_out/ortho_bmcFS_TX_sig.csv")
sortbmcFS_TX_sig<-cbind(bmcFS_TX_sig[match(rncFS_TX_sig, bmcFS_TX_sig$ensembl_gene_id),])
bmcFS_TX<-cbind(cFS_TX_sig, bmcFS_TX_sig)
write.csv(bmcFS_TX, "mRNA_out/bmcFS_TX_sig.csv")

#######################################################################################
#The difference between Finnsheep and Texel with flushing diet
########################################################################################
fFS_TX<-results(dds, contrast=c("group", "FS.F", "TX.F"))
fFS_TX_sig<-subset(fFS_TX, padj < 0.05)
nrow(fFS_TX_sig) #621
summary(fFS_TX_sig)
write.csv(fFS_TX_sig, "mRNA_out/fFS_TX_sig.csv")
rnfFS_TX_sig<-rownames(fFS_TX_sig)
gofFS_TX_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id", "go_id"), values=rnfFS_TX_sig, mart=ensembl)
write.csv(gofFS_TX_sig, "mRNA_out/gofFS_TX_sig.csv")
bmfFS_TX_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id",
                                                                                      "external_gene_name", "chromosome_name", "description"), values=rnfFS_TX_sig, mart=ensembl)
ortho_bmfFS_TX_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id", "external_gene_name",
                                                                                            "btaurus_homolog_ensembl_gene", "sscrofa_homolog_ensembl_gene", "mmusculus_homolog_ensembl_gene", "hsapiens_homolog_ensembl_gene"), values=rnfFS_TX_sig, mart=ensembl)
write.csv(ortho_bmfFS_TX_sig, "mRNA_out/ortho_bmfFS_TX_sig.csv")
sortbmfFS_TX_sig<-cbind(bmfFS_TX_sig[match(rnfFS_TX_sig, bmfFS_TX_sig$ensembl_gene_id),])
bmfFS_TX<-cbind(fFS_TX_sig, bmfFS_TX_sig)
write.csv(bmfFS_TX, "mRNA_out/bmfFS_TX_sig.csv")

##############################################################################################
#The difference between Finnsheep and F1 with control diet
##############################################################################################
cFS_F1 <- results(dds, contrast=c("group", "FS.C", "F1.C"))
cFS_F1_sig<-subset(cFS_F1, padj < 0.05)
nrow(cFS_F1_sig) #2
summary(cFS_F1_sig)
write.csv(cFS_F1_sig, "mRNA_out/cFS_F1_sig.csv")
rncFS_F1_sig<-rownames(cFS_F1_sig)
gocFS_F1_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id", "go_id"), values=rncFS_F1_sig, mart=ensembl)
write.csv(gocFS_F1_sig, "mRNA_out/gocFS_F1_sig.csv")
bmcFS_F1_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id",
                                                                                      "external_gene_name", "chromosome_name", "description"), values=rncFS_F1_sig, mart=ensembl)
ortho_bmcFS_F1_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id", "external_gene_name",
                                                                                            "btaurus_homolog_ensembl_gene", "sscrofa_homolog_ensembl_gene", "mmusculus_homolog_ensembl_gene", "hsapiens_homolog_ensembl_gene"), values=rncFS_F1_sig, mart=ensembl)
write.csv(ortho_bmcFS_F1_sig, "mRNA_out/ortho_bmcFS_F1_sig.csv")
sortbmcFS_F1_sig<-cbind(bmcFS_F1_sig[match(rncFS_F1_sig, bmcFS_F1_sig$ensembl_gene_id),])
bmcFS_F1<-cbind(cFS_F1_sig, bmcFS_F1_sig)
write.csv(bmcFS_F1, "mRNA_out/bmcFS_F1_sig.csv")

##############################################################################################
#The difference between Finnsheep and F1 with flushing diet
##############################################################################################
fFS_F1 <- results(dds, contrast=c("group", "FS.F", "F1.F"))
fFS_F1_sig<-subset(fFS_F1, padj < 0.05)
nrow(fFS_F1_sig) #51
summary(fFS_F1_sig)
write.csv(fFS_F1_sig, "mRNA_out/fFS_F1_sig.csv")
rnfFS_F1_sig<-rownames(fFS_F1_sig)
gofFS_F1_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id", "go_id"), values=rnfFS_F1_sig, mart=ensembl)
write.csv(gofFS_F1_sig, "mRNA_out/gofFS_F1_sig.csv")
bmfFS_F1_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id",
                                                                                      "external_gene_name", "chromosome_name", "description"), values=rnfFS_F1_sig, mart=ensembl)
ortho_bmfFS_F1_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id", "external_gene_name",
                                                                                            "btaurus_homolog_ensembl_gene", "sscrofa_homolog_ensembl_gene", "mmusculus_homolog_ensembl_gene", "hsapiens_homolog_ensembl_gene"), values=rnfFS_F1_sig, mart=ensembl)
write.csv(ortho_bmfFS_F1_sig, "mRNA_out/ortho_bmfFS_F1_sig.csv")
sortbmfFS_F1_sig<-cbind(bmfFS_F1_sig[match(rnfFS_F1_sig, bmfFS_F1_sig$ensembl_gene_id),])
bmfFS_F1<-cbind(fFS_F1_sig, bmfFS_F1_sig)
write.csv(bmfFS_F1, "mRNA_out/bmfFS_F1_sig.csv")

################################################################################################
#The difference between Texel and F1 with control diet
################################################################################################
cTX_F1 <- results(dds, contrast=c("group", "TX.C", "F1.C"))
cTX_F1_sig<-subset(cTX_F1, padj < 0.05)
nrow(cTX_F1_sig) #60
summary(cTX_F1_sig)
write.csv(cTX_F1_sig, "mRNA_out/cTX_F1_sig.csv")
rncTX_F1_sig<-rownames(cTX_F1_sig)
gocTX_F1_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id", "go_id"), values=rncTX_F1_sig, mart=ensembl)
write.csv(gocTX_F1_sig, "mRNA_out/gocTX_F1_sig.csv")
bmcTX_F1_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id",
                                                                                      "external_gene_name", "chromosome_name", "description"), values=rncTX_F1_sig, mart=ensembl)
ortho_bmcTX_F1_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id", "external_gene_name",
                                                                                            "btaurus_homolog_ensembl_gene", "sscrofa_homolog_ensembl_gene", "mmusculus_homolog_ensembl_gene", "hsapiens_homolog_ensembl_gene"), values=rncTX_F1_sig, mart=ensembl)
write.csv(ortho_bmcTX_F1_sig, "mRNA_out/ortho_bmcTX_F1_sig.csv")
sortbmcTX_F1_sig<-cbind(bmcTX_F1_sig[match(rncTX_F1_sig, bmcTX_F1_sig$ensembl_gene_id),])
bmcTX_F1<-cbind(cTX_F1_sig, bmcTX_F1_sig)
write.csv(bmcTX_F1, "mRNA_out/bmcTX_F1_sig.csv")

################################################################################################
#The difference between Texel and F1 with flushing diet
################################################################################################
fTX_F1 <- results(dds, contrast=c("group", "TX.F", "F1.F"))
fTX_F1_sig<-subset(fTX_F1, padj<0.05)
nrow(fTX_F1_sig) #311
summary(fTX_F1_sig)
write.csv(fTX_F1_sig, "mRNA_out/fTX_F1_sig.csv")
rnfTX_F1_sig<-rownames(fTX_F1_sig)
gofTX_F1_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id", "go_id"), values=rnfTX_F1_sig, mart=ensembl)
write.csv(gofTX_F1_sig, "mRNA_out/gofTX_F1_sig.csv")
bmfTX_F1_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id",
                                                                                      "external_gene_name", "chromosome_name", "description"), values=rnfTX_F1_sig, mart=ensembl)
ortho_bmfTX_F1_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id", "external_gene_name",
                                                                                            "btaurus_homolog_ensembl_gene", "sscrofa_homolog_ensembl_gene", "mmusculus_homolog_ensembl_gene", "hsapiens_homolog_ensembl_gene"), values=rnfTX_F1_sig, mart=ensembl)
write.csv(ortho_bmfTX_F1_sig, "mRNA_out/ortho_bmfTX_F1_sig.csv")
sortbmfTX_F1_sig<-cbind(bmfTX_F1_sig[match(rnfTX_F1_sig, bmfTX_F1_sig$ensembl_gene_id),])
bmfTX_F1<-cbind(fTX_F1_sig, bmfTX_F1_sig)
write.csv(bmfTX_F1, "mRNA_out/bmfTX_F1_sig.csv")


#Merge all DE tables, collect unique IDs, fetch to cluego app and find out how may IDs are not available to be annotated there
library(plyr)
df1<-data.frame(dFS_Diet_sig)
df1$rn <-rownames(df1)
df2<-data.frame(dTX_Diet_sig)
df2$rn<-rownames(df2)
df3<-data.frame(dF1_Diet_sig)
df3$rn<-rownames(df3)
df4<-data.frame(dFS_TX_sig)
df4$rn<-rownames(df4)
df5<-data.frame(dFS_F1_sig)
df5$rn<-rownames(df5)
df6<-data.frame(dTX_F1_sig)
df6$rn<-rownames(df6)
df7<-data.frame(cFS_TX_sig)
df7$rn<-rownames(df7)
df8<-data.frame(fFS_TX_sig)
df8$rn<-rownames(df8)
df9<-data.frame(cFS_F1_sig)
df9$rn<-rownames(df9)
df10<-data.frame(fFS_F1_sig)
df10$rn<-rownames(df10)
df11<-data.frame(cTX_F1_sig)
df11$rn<-rownames(df11)
df12<-data.frame(fTX_F1_sig)
df12$rn<-rownames(df12)
all_rownames<-join_all(list(df2, df3, df4, df5, df6, df7, df8, df9, df10, df11, df12), by = 'rn', type = 'full')
nrow(all_rownames)
all_rn<-all_rownames$rn
write.csv(all_rn, "mRNA_out/all_rownames.csv")

bm_all_de_entrez<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("entrezgene", "ensembl_gene_id",
                                                                                          "ensembl_transcript_id", "ensembl_peptide_id"), values=all_rn, mart=ensembl)
write.csv((bm_all_de_entrez), "mRNA_out/bm_all_de_entrez.csv")

#It was found out that 
#all_cluego<-read.csv("All_de_to_be_annotated.csv")
#bm_all_cluego_entrez<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("entrezgene", "ensembl_gene_id",
                                                                                              "ensembl_transcript_id", "ensembl_peptide_id"), values=all_cluego, mart=ensembl)

#write.csv(na.omit(bm_all_cluego_entrez), "mRNA_out/bm_all_cluego_entrez.csv")


attr<-listAttributes(ensembl)
head(attr, n=50)
head(dFS_TX)
head(dFS_F1)
all_de <- merge(dFS_TX, dFS_F1)
nrow(all_de)
colnames(dFS_TX_sig)
#Dispersion plot
res <- results(dds)
png("qc-dispersions.png", 1000, 1000, pointsize=20)
plotDispEsts(dds, main="Dispersion plot")
dev.off()

# Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds)
head(assay(rld))
hist(assay(rld))

# Colors for plots below
## Ugly:
## (mycols <- 1:length(unique(condition)))
## Use RColorBrewer, better
require(graphics); require(grDevices)
library(RColorBrewer)
(mycols <- brewer.pal(8, "Dark2")[1:length(unique(group))])
# Sample distance heatmap
sampleDists <- as.matrix(dist(t(assay(rld))))
sampleDists
library(gplots)
png("qc-heatmap-samples.png", w=1000, h=1000, pointsize=10)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors = mycols[group], RowSideColors = mycols[group],
          margin=c(10, 10), main="Sample Distance Matrix")
dev.off()
#heatmap can be made more colorful with ColSideColors=mycols[group], RowSideColors=mycols[group],

# Principal components analysis
rld_pca <- function (rld, intgroup = "group", ntop = 500, colors=NULL, legendpos="bottomleft", main="PCA Biplot", textcx=1, ...) {
  require(genefilter)
  require(calibrate)
  require(RColorBrewer)
  rv = rowVars(assay(rld))
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld)[select, ]))
  fac = factor(apply(as.data.frame(colData(rld)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
  if (is.null(colors)) {
    if (nlevels(fac) >= 3) {
      colors = brewer.pal(nlevels(fac), "Paired")
    }   else {
      colors = c("black", "red")
    }
  }
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC1 (",as.character(pc2var),"%)")
  plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main=main, ...)
  with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))
  legend(legendpos, legend=levels(fac), col=colors, pch=20)
  #     rldyplot(PC2 ~ PC1, groups = fac, data = as.data.frame(pca$rld),
  #            pch = 16, cerld = 2, aspect = "iso", col = colours, main = draw.key(key = list(rect = list(col = colours),
  #                                                                                         terldt = list(levels(fac)), rep = FALSE)))
}
png("qc-pca.png", 1000, 1000, pointsize=20)
rld_pca(rld, colors=mycols, intgroup="group", xlim=c(-75, 35))
dev.off()


# Get differential expression results
res <- results(dds)
res
table(res$padj<0.05)
## Order by adjusted p-value
res <- res[order(res$padj), ]
## Merge with normalized count data
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)
## Write results
write.csv(resdata, file="diffexpr-results.csv")

## Examine plot of p-values
hist(res$pvalue, breaks=50, col="grey")

## Examine independent filtering
attr(res, "filterThreshold")
plot(attr(res,"filterNumRej"), type="b", xlab="quantiles of baseMean", ylab="number of rejections")

## MA plot
## Could do with built-in DESeq2 function:
## DESeq2::plotMA(dds, ylim=c(-1,1), cex=1)
## I like mine better:
library(calibrate)
maplot <- function (res, thresh=0.05, labelsig=TRUE, textcx=1, ...) {
  with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
  with(subset(res, padj<thresh), points(baseMean, log2FoldChange, col="red", pch=20, cex=1.5))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<thresh), textxy(baseMean, log2FoldChange, labs=Gene, cex=textcx, col=2))
  }
}
png("diffexpr-maplot.png", 1500, 1000, pointsize=20)
maplot(resdata, main="MA Plot")
dev.off()

## Volcano plot with "significant" genes labeled, one example
volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}
png("diffexpr-volcanoplot.png", 1200, 1000, pointsize=20)
volcanoplot(resdata, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-2.3, 2))
dev.off()

select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:30]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(counts(dds,normalized=TRUE)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10,6))

heatmap.2(assay(rld)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10, 6))
