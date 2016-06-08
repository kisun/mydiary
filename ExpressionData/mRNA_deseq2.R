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
ensembl2 = useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "btaurus_gene_ensembl")
ensembl3 = useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
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
sampleinfo
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
head(res_baseMean5)
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
nrow(dTX_Diet_sig)#118
summary(dTX_Diet_sig)
write.csv(dTX_Diet_sig, "ExpressionData/mRNA_out/dTX_Diet_sig.csv")
rndTX_Diet_sig<-rownames(dTX_Diet_sig)
#Retrive GO IDs although GO annotations and KEGG pathway analysis are done in Cytoscape using ClueGO plugin
godTX_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id", "go_id"), values=rndTX_Diet_sig, mart=ensembl)
write.csv(godTX_sig, "ExpressionData/mRNA_out/godTX_sig.csv")
#We retrieve some basic information for each differentially expressed genes using biomart.
bmdTX_Diet_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id",
                                                                                        "external_gene_name","chromosome_name", "description"), values=rndTX_Diet_sig, mart=ensembl)
#In case some unknown genes in sheep are annotated in other species, it has also been useful to do the GO annotation and KEGG pathways using cow homologs as it's genome is much more annotated.
ortho_bmdTX_Diet_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id", "external_gene_name",
                                                                                              "btaurus_homolog_ensembl_gene", "sscrofa_homolog_ensembl_gene", "mmusculus_homolog_ensembl_gene", "hsapiens_homolog_ensembl_gene"), values=rndTX_Diet_sig, mart=ensembl)


write.csv(ortho_bmdTX_Diet_sig, "ExpressionData/mRNA_out/ortho_bmdTX_diet_sig.csv")
sortbmdTX_Diet_sig<-cbind(bmdTX_Diet_sig[match(rndTX_Diet_sig, bmdTX_Diet_sig$ensembl_gene_id),])
bmdTX_Diet<-cbind(dTX_Diet_sig, bmdTX_Diet_sig)
write.csv(bmdTX_Diet, "ExpressionData/mRNA_out/bmdTX_Diet_sig.csv")

###################################################################################
#The effect of diet in F1-cross
###################################################################################
dF1_Diet<-results(dds, contrast=c("group", "F1.F", "F1.C"))
dF1_Diet_sig<-subset(dF1_Diet, padj<0.05)
nrow(dF1_Diet_sig)#25
summary(dF1_Diet_sig)
write.csv(dF1_Diet_sig, "ExpressionData/mRNA_out/dF1_Diet_sig.csv")
rndF1_Diet_sig<-rownames(dF1_Diet_sig)
godF1_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id", "go_id"), values=rndF1_Diet_sig, mart=ensembl)
write.csv(godF1_sig, "ExpressionData/mRNA_out/godF1_sig.csv")
bmdF1_Diet_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id",
                                                                                        "external_gene_name", "chromosome_name", "description"), values=rndF1_Diet_sig, mart=ensembl)
ortho_bmdF1_Diet_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id", "external_gene_name",
                                                                                              "btaurus_homolog_ensembl_gene", "sscrofa_homolog_ensembl_gene", "mmusculus_homolog_ensembl_gene", "hsapiens_homolog_ensembl_gene"), values=rndF1_Diet_sig, mart=ensembl)
write.csv(ortho_bmdF1_Diet_sig, "ExpressionData/mRNA_out/ortho_bmdF1_diet_sig.csv")
sortbmdF1_Diet_sig<-cbind(bmdF1_Diet_sig[match(rndF1_Diet_sig, bmdF1_Diet_sig$ensembl_gene_id),])
bmdF1_Diet<-cbind(dF1_Diet_sig, bmdF1_Diet_sig)
write.csv(bmdF1_Diet, "ExpressionData/mRNA_out/bmdF1_Diet_sig.csv")

###################################################################################
#The difference in Finnsheep and Texel due to diet effect
###################################################################################
dFS_TX<-results(dds, contrast=list(c("groupFS.F", "groupTX.C"), c("groupFS.C", "groupTX.F")))
dFS_TX_sig<-subset(dFS_TX, padj<0.05)
nrow(dFS_TX_sig)#38
summary(dFS_TX_sig)
write.csv(dFS_TX_sig, "ExpressionData/mRNA_out/dFS_TX_sig.csv")
rndFS_TX_sig<-rownames(dFS_TX_sig)
godFS_TX_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id", "go_id"), values=rndFS_TX_sig, mart=ensembl)
write.csv(godFS_TX_sig, "ExpressionData/mRNA_out/godFS_TX_sig.csv")
bmdFS_TX_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id",
                                                                                      "external_gene_name", "chromosome_name", "description"), values=rndFS_TX_sig, mart=ensembl)
ortho_bmdFS_TX_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id", "external_gene_name",
                                                                                            "btaurus_homolog_ensembl_gene", "sscrofa_homolog_ensembl_gene", "mmusculus_homolog_ensembl_gene", "hsapiens_homolog_ensembl_gene"), values=rndFS_TX_sig, mart=ensembl)
write.csv(ortho_bmdFS_TX_sig, "ExpressionData/mRNA_out/ortho_bmdFS_TX_sig.csv")
sortbmdFS_TX_sig<-cbind(bmdFS_TX_sig[match(rndFS_TX_sig, bmdFS_TX_sig$ensembl_gene_id),])
bmdFS_TX<-cbind(dFS_TX_sig, bmdFS_TX_sig)
write.csv(bmdFS_TX, "ExpressionData/mRNA_out/bmdFS_TX_sig.csv")

###################################################################################
#The difference in Finsheep and F1 due to diet effect
###################################################################################
dFS_F1<-results(dds, contrast=list(c("groupFS.F", "groupF1.C"), c("groupFS.C", "groupF1.F")))
dFS_F1_sig<-subset(dFS_F1, padj<0.05)
nrow(dFS_F1_sig) #5
summary(dFS_F1_sig)
write.csv(dFS_F1_sig, "ExpressionData/mRNA_out/dFS_F1_sig.csv")
rndFS_F1_sig<-rownames(dFS_F1_sig)
godFS_F1_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id", "go_id"), values=rndFS_F1_sig, mart=ensembl)
write.csv(godFS_F1_sig, "ExpressionData/mRNA_out/godFS_F1_sig.csv")
bmdFS_F1_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id",
                                                                                      "external_gene_name", "chromosome_name", "description"), values=rndFS_F1_sig, mart=ensembl)
ortho_bmdFS_F1_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id", "external_gene_name",
                                                                                            "btaurus_homolog_ensembl_gene", "sscrofa_homolog_ensembl_gene", "mmusculus_homolog_ensembl_gene", "hsapiens_homolog_ensembl_gene"), values=rndFS_F1_sig, mart=ensembl)
write.csv(ortho_bmdFS_F1_sig, "ExpressionData/mRNA_out/ortho_bmdFS_F1_sig.csv")
sortbmdFS_TX_sig<-cbind(bmdFS_TX_sig[match(rndFS_TX_sig, bmdFS_TX_sig$ensembl_gene_id),])
sortbmdFS_F1_sig<-cbind(bmdFS_F1_sig[match(rndFS_F1_sig, bmdFS_F1_sig$ensembl_gene_id),])
bmdFS_F1<-cbind(dFS_F1_sig, bmdFS_F1_sig)
write.csv(bmdFS_F1, "ExpressionData/mRNA_out/bmdFS_F1_sig.csv")

###################################################################################
#The difference in Texel and F1 due to diet effect
###################################################################################
dTX_F1<-results(dds, contrast=list(c("groupTX.F", "groupF1.C"), c("groupTX.C", "groupF1.F")))
dTX_F1_sig<-subset(dTX_F1, padj<0.05)
nrow(dTX_F1_sig)#68
summary(dTX_F1_sig)
write.csv(dTX_F1_sig, "ExpressionData/mRNA_out/dTX_F1_sig.csv")
rndTX_F1_sig<-rownames(dTX_F1_sig)
godTX_F1_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id", "go_id"), values=rndTX_F1_sig, mart=ensembl)
write.csv(godTX_F1_sig, "ExpressionData/mRNA_out/godTX_F1_sig.csv")
bmdTX_F1_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id",
                                                                                      "external_gene_name", "chromosome_name", "description"), values=rndTX_F1_sig, mart=ensembl)
ortho_bmdTX_F1_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id", "external_gene_name",
                                                                                            "btaurus_homolog_ensembl_gene", "sscrofa_homolog_ensembl_gene", "mmusculus_homolog_ensembl_gene", "hsapiens_homolog_ensembl_gene"), values=rndTX_F1_sig, mart=ensembl)
write.csv(ortho_bmdTX_F1_sig, "ExpressionData/mRNA_out/ortho_bmdTX_F1_sig.csv")
sortbmdTX_F1_sig<-cbind(bmdTX_F1_sig[match(rndTX_F1_sig, bmdTX_F1_sig$ensembl_gene_id),])
bmdTX_F1<-cbind(dTX_F1_sig, bmdTX_F1_sig)
write.csv(bmdTX_F1, "ExpressionData/mRNA_out/bmdTX_F1_sig.csv")

#######################################################################################
#The difference between Finnsheep and Texel with Control diet
#######################################################################################
cFS_TX <- results(dds, contrast=c("group", "FS.C", "TX.C"))
cFS_TX_sig<-subset(cFS_TX, padj < 0.05)
nrow(cFS_TX_sig) #3
summary(cFS_TX_sig)
write.csv(cFS_TX_sig, "ExpressionData/mRNA_out/cFS_TX_sig.csv")
rncFS_TX_sig<-rownames(cFS_TX_sig)
gocFS_TX_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id", "go_id"), values=rncFS_TX_sig, mart=ensembl)
write.csv(gocFS_TX_sig, "ExpressionData/mRNA_out/gocFS_TX_sig.csv")
bmcFS_TX_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id",
                                                                                      "external_gene_name", "chromosome_name", "description"), values=rncFS_TX_sig, mart=ensembl)
ortho_bmcFS_TX_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id", "external_gene_name",
                                                                                            "btaurus_homolog_ensembl_gene", "sscrofa_homolog_ensembl_gene", "mmusculus_homolog_ensembl_gene", "hsapiens_homolog_ensembl_gene"), values=rncFS_TX_sig, mart=ensembl)
write.csv(ortho_bmcFS_TX_sig, "ExpressionData/mRNA_out/ortho_bmcFS_TX_sig.csv")
sortbmcFS_TX_sig<-cbind(bmcFS_TX_sig[match(rncFS_TX_sig, bmcFS_TX_sig$ensembl_gene_id),])
bmcFS_TX<-cbind(cFS_TX_sig, bmcFS_TX_sig)
write.csv(bmcFS_TX, "ExpressionData/mRNA_out/bmcFS_TX_sig.csv")
#Biomartstuff
hs_id_cFS_TX<-ortho_bmcFS_TX_sig$hsapiens_homolog_ensembl_gene
head(hs_id_cFS_TX)
hs_id1_cFS_TX<-unique(hs_id_cFS_TX)
head(hs_id1_cFS_TX)
hs_ens_geneid_cFS_TX<-getBM(c("hsapiens_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id", "ensembl_transcript_id"), values=hs_id1_cFS_TX, mart = ensembl3)
names(hs_ens_geneid_cFS_TX)[2]<-paste("hs_transcript")
head(hs_ens_geneid_cFS_TX)
write.csv(hs_ens_geneid_cFS_TX, "ExpressionData/targetscan/cFS-TX/transcript_geneid_cFS_TX.csv")

#miR6529
ts_miR6529<-read.csv("ExpressionData/targetscan/cFS-TX/TargetScan7.0__miR-6529b.predicted_targets.txt", header=FALSE, sep = "\t", stringsAsFactors = FALSE)
ts_miR6529_3<-ts_miR6529[which(ts_miR6529$V11 >= -0.3),]
ts_miR6529_3$V2<-sub("\\.\\d$", "", ts_miR6529_3$V2)
colnames(ts_miR6529_3)<-ts_miR6529_3[1,]
ts_miR6529_3<-ts_miR6529_3[-1,]
tail(ts_miR6529_3)
rownames(ts_miR6529_3)<-NULL
ts_miR6529_3<-head(ts_miR6529_3, -1)
write.csv(ts_miR6529_3, "ExpressionData/targetscan/cFS-TX/ts_miR6529_3_mod.csv")
names(ts_miR6529_3)[2]<-paste("hs_transcript")
list_miR6529_3<-list(hs_ens_geneid_cFS_TX, ts_miR6529_3)
merged_list_miR6529_3<-Reduce(function(...) merge(..., by=c("hs_transcript"), all=F), list_miR6529_3)
merged_list_miR6529_3
write.csv(merged_list_miR6529_3, "ExpressionData/targetscan/cFS-TX/DE_targets_DE_miRNAs_miR6529_3.csv")

#miR192
ts_miR192<-read.csv("ExpressionData/targetscan/cFS-TX/TargetScan7.0__miR-192_215.predicted_targets.txt", header=FALSE, sep = "\t", stringsAsFactors = FALSE)
ts_miR192_3<-ts_miR192[which(ts_miR192$V15 >= -0.3),]
ts_miR192_3$V2<-sub("\\.\\d$", "", ts_miR192_3$V2)
colnames(ts_miR192_3)<-ts_miR192_3[1,]
ts_miR192_3<-ts_miR192_3[-1,]
tail(ts_miR192_3, 20)
rownames(ts_miR192_3)<-NULL
ts_miR192_3<-head(ts_miR192_3, -6)
write.csv(ts_miR192_3, "ExpressionData/targetscan/cFS-TX/ts_miR192_3_mod.csv")
names(ts_miR192_3)[2]<-paste("hs_transcript")
list_miR192_3<-list(hs_ens_geneid_cFS_TX, ts_miR192_3)
merged_list_miR192_3<-Reduce(function(...) merge(..., by=c("hs_transcript"), all=F), list_miR192_3)
merged_list_miR192_3
write.csv(merged_list_miR192_3, "ExpressionData/targetscan/cFS-TX/DE_targets_DE_miRNAs_miR192_3.csv")

#miR151
ts_miR151<-read.csv("ExpressionData/targetscan/cFS-TX/TargetScan7.0__miR-151-5p.predicted_targets.txt", header=FALSE, sep = "\t", stringsAsFactors = FALSE)
ts_miR151_3<-ts_miR151[which(ts_miR151$V11 >= -0.3),]
ts_miR151_3$V2<-sub("\\.\\d$", "", ts_miR151_3$V2)
colnames(ts_miR151_3)<-ts_miR151_3[1,]
ts_miR151_3<-ts_miR151_3[-1,]
tail(ts_miR151_3, 40)
rownames(ts_miR151_3)<-NULL
ts_miR151_3<-head(ts_miR151_3, -31)
write.csv(ts_miR151_3, "ExpressionData/targetscan/cFS-TX/ts_miR151_3_mod.csv")
names(ts_miR151_3)[2]<-paste("hs_transcript")
list_miR151_3<-list(hs_ens_geneid_cFS_TX, ts_miR151_3)
merged_list_miR151_3<-Reduce(function(...) merge(..., by=c("hs_transcript"), all=F), list_miR151_3)
merged_list_miR151_3
write.csv(merged_list_miR151_3, "ExpressionData/targetscan/cFS-TX/DE_targets_DE_miRNAs_miR151_3.csv")






#######################################################################################
#The difference between Finnsheep and Texel with flushing diet
########################################################################################
fFS_TX<-results(dds, contrast=c("group", "FS.F", "TX.F"))
fFS_TX_sig<-subset(fFS_TX, padj < 0.05)
nrow(fFS_TX_sig) #600
summary(fFS_TX_sig)
write.csv(fFS_TX_sig, "ExpressionData/mRNA_out/fFS_TX_sig.csv")
rnfFS_TX_sig<-rownames(fFS_TX_sig)
gofFS_TX_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id", "go_id"), values=rnfFS_TX_sig, mart=ensembl)
write.csv(gofFS_TX_sig, "ExpressionData/mRNA_out/gofFS_TX_sig.csv")
bmfFS_TX_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id",
                                                                                      "external_gene_name", "chromosome_name", "description"), values=rnfFS_TX_sig, mart=ensembl)
ortho_bmfFS_TX_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id", "external_gene_name",
                                                                                            "btaurus_homolog_ensembl_gene", "sscrofa_homolog_ensembl_gene", "mmusculus_homolog_ensembl_gene", "hsapiens_homolog_ensembl_gene", "hsapiens_homolog_canonical_transcript_protein"), values=rnfFS_TX_sig, mart=ensembl)
write.csv(ortho_bmfFS_TX_sig, "ExpressionData/mRNA_out/ortho_bmfFS_TX_sig.csv")

sortbmfFS_TX_sig<-cbind(bmfFS_TX_sig[match(rnfFS_TX_sig, bmfFS_TX_sig$ensembl_gene_id),])
bmfFS_TX<-cbind(fFS_TX_sig, bmfFS_TX_sig)
write.csv(bmfFS_TX, "ExpressionData/mRNA_out/bmfFS_TX_sig.csv")

#Biomartstuff
hs_id<-ortho_bmfFS_TX_sig$hsapiens_homolog_ensembl_gene
head(hs_id)
hs_id1<-unique(hs_id)
head(hs_id1)
nrow(hs_id1)
hs_ens_geneid<-getBM(c("hsapiens_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id", "ensembl_transcript_id"), values=hs_id1, mart = ensembl3)
write.csv(hs_ens_geneid, "ExpressionData/targetscan/fFS-TX/transcript_geneid_fFS_TX.csv")
names(hs_ens_geneid)[2]<-paste("hs_transcript")

#miR432
ts_miR432<-read.csv("ExpressionData/targetscan/fFS-TX/TargetScan7.0__miR-432.predicted_targets.csv", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
ts_miR432_3<-ts_miR432[which(ts_miR432$V11 >= -0.3),]
ts_miR432_3$V2<-sub("\\.\\d$", "", ts_miR432_3$V2)
colnames(ts_miR432_3)<-ts_miR432_3[1,]
ts_miR432_3<-ts_miR432_3[-1,]
tail(ts_miR432_3)
rownames(ts_miR432_3)<-NULL
ts_miR432_3<-head(ts_miR432_3, -4)
write.csv(ts_miR432_3, "ExpressionData/targetscan/fFS-TX/ts_miR432_3_mod.csv")
names(ts_miR432_3)[2]<-paste("hs_transcript")
lista<-list(hs_ens_geneid, ts_miR432_3)
merged_lista = Reduce(function(...) merge(..., by=c("hs_transcript"), all=F), lista)
write.csv(merged_lista, "ExpressionData/targetscan/fFS-TX/DE_targets_DE_miRNAs_miR432_3.csv")


#miR-433-3p
ts_miR433<-read.csv("ExpressionData/targetscan/fFS-TX/TargetScan7.0__miR-433-3p.predicted_targets.txt", header=FALSE, sep = "\t", stringsAsFactors = FALSE)
ts_miR433_3<-ts_miR433[which(ts_miR433$V11 >= -0.3),]
ts_miR433_3$V2<-sub("\\.\\d$", "", ts_miR433_3$V2)
colnames(ts_miR433_3)<-ts_miR433_3[1,]
ts_miR433_3<-ts_miR433_3[-1,]
tail(ts_miR433_3)
rownames(ts_miR433_3)<-NULL
ts_miR433_3<-head(ts_miR433_3, -4)
write.csv(ts_miR433_3, "ExpressionData/targetscan/fFS-TX/ts_miR433_3_mod.csv")
names(ts_miR433_3)[2]<-paste("hs_transcript")
list_miR433_3<-list(hs_ens_geneid, ts_miR433_3)
merged_list_miR433_3<-Reduce(function(...) merge(..., by=c("hs_transcript"), all=F), list_miR433_3)
merged_list_miR433_3
write.csv(merged_list_miR433_3, "ExpressionData/targetscan/fFS-TX/DE_targets_DE_miRNAs_miR433_3.csv")

#miR-2285-m-3
ts_miR2285m3<-read.csv("ExpressionData/targetscan/fFS-TX/TargetScan7.0__miR-2285m.predicted_targets.txt", header=FALSE, sep = "\t", stringsAsFactors = FALSE)
ts_miR2285m3_3<-ts_miR2285m3[which(ts_miR2285m3$V11 >= -0.3),]
ts_miR2285m3_3$V2<-sub("\\.\\d$", "", ts_miR2285m3_3$V2)
colnames(ts_miR2285m3_3)<-ts_miR2285m3_3[1,]
ts_miR2285m3_3<-ts_miR2285m3_3[-1,]
tail(ts_miR2285m3_3, n=20)
ts_miR2285m3_3<-head(ts_miR2285m3_3, -7)
rownames(ts_miR2285m3_3)<-NULL
write.csv(ts_miR2285m3_3, "ExpressionData/targetscan/fFS-TX/ts_miR2285m3_3_mod.csv")
names(ts_miR2285m3_3)[2]<-paste("hs_transcript")
list_miR2285m3_3<-list(hs_ens_geneid, ts_miR2285m3_3)
merged_list_miR2285m3_3<-Reduce(function(...) merge(..., by=c("hs_transcript"), all=F), list_miR2285m3_3)
merged_list_miR2285m3_3
write.csv(merged_list_miR2285m3_3, "ExpressionData/targetscan/fFS-TX/DE_targets_DE_miRNAs_miR2285m3_3.csv")

#miR-326
ts_miR326<-read.csv("ExpressionData/targetscan/fFS-TX/TargetScan7.0__miR-326.predicted_targets.txt", header=FALSE, sep = "\t", stringsAsFactors = FALSE)
ts_miR326_3<-ts_miR326[which(ts_miR326$V11 >= -0.3),]
ts_miR326_3$V2<-sub("\\.\\d$", "", ts_miR326_3$V2)
colnames(ts_miR326_3)<-ts_miR326_3[1,]
ts_miR326_3<-ts_miR326_3[-1,]
tail(ts_miR326_3, 20)
ts_miR326_3<-head(ts_miR326_3, -11)
dim(ts_miR326)
dim(ts_miR326_3)
rownames(ts_miR326_3)<-NULL
write.csv(ts_miR326_3, "ExpressionData/targetscan/fFS-TX/ts_miR326_3_mod.csv")
names(ts_miR326_3)[2]<-paste("hs_transcript")
list_miR326_3<-list(hs_ens_geneid, ts_miR326_3)
merged_list_miR326_3<-Reduce(function(...) merge(..., by=c("hs_transcript"), all=F), list_miR326_3)
merged_list_miR326_3
write.csv(merged_list_miR326_3, "ExpressionData/targetscan/fFS-TX/DE_targets_DE_miRNAs_miR326_3.csv")

#miR-130a
ts_miR130a<-read.csv("ExpressionData/targetscan/fFS-TX/TargetScan7.0__miR-130_301_454.predicted_targets.txt", header=FALSE, sep = "\t", stringsAsFactors = FALSE)
ts_miR130a_3<-ts_miR130a[which(ts_miR130a$V15 >= -0.3),]
ts_miR130a_3$V2<-sub("\\.\\d$", "", ts_miR130a_3$V2)
colnames(ts_miR130a_3)<-ts_miR130a_3[1,]
ts_miR130a_3<-ts_miR130a_3[-1,]
tail(ts_miR130a_3, 40)
ts_miR130a_3<-head(ts_miR130a_3, -30)
rownames(ts_miR130a_3)<-NULL
write.csv(ts_miR130a_3, "ExpressionData/targetscan/fFS-TX/ts_miR130a_3_mod.csv")
names(ts_miR130a_3)[2]<-paste("hs_transcript")
list_miR130a_3<-list(hs_ens_geneid, ts_miR130a_3)
merged_list_miR130a_3<-Reduce(function(...) merge(..., by=c("hs_transcript"), all=F), list_miR130a_3)
merged_list_miR130a_3
write.csv(merged_list_miR130a_3, "ExpressionData/targetscan/fFS-TX/DE_targets_DE_miRNAs_miR130a_3.csv")

#miR-2284
ts_miR2284<-read.csv("ExpressionData/targetscan/fFS-TX/TargetScan7.0__miR-2284_2285.predicted_targets.txt", header=FALSE, sep = "\t", stringsAsFactors = FALSE)
ts_miR2284_3<-ts_miR2284[which(ts_miR2284$V11 >= -0.3),]
ts_miR2284_3$V2<-sub("\\.\\d$", "", ts_miR2284_3$V2)
colnames(ts_miR2284_3)<-ts_miR2284_3[1,]
ts_miR2284_3<-ts_miR2284_3[-1,]
tail(ts_miR2284_3, 20)
ts_miR2284_3<-head(ts_miR2284_3, -18)
rownames(ts_miR2284_3)<-NULL
write.csv(ts_miR2284_3, "ExpressionData/targetscan/fFS-TX/ts_miR2284_3_mod.csv")
names(ts_miR2284_3)[2]<-paste("hs_transcript")
list_miR2284_3<-list(hs_ens_geneid, ts_miR2284_3)
merged_list_miR2284_3<-Reduce(function(...) merge(..., by=c("hs_transcript"), all=F), list_miR2284_3)
merged_list_miR2284_3
write.csv(merged_list_miR2284_3, "ExpressionData/targetscan/fFS-TX/DE_targets_DE_miRNAs_miR2284_3.csv")

#miR-1468
ts_miR1468<-read.csv("ExpressionData/targetscan/fFS-TX/TargetScan7.0__miR-1468.predicted_targets.txt", header=FALSE, sep = "\t", stringsAsFactors = FALSE)
ts_miR1468_3<-ts_miR1468[which(ts_miR1468$V11 >= -0.3),]
ts_miR1468_3$V2<-sub("\\.\\d$", "", ts_miR1468_3$V2)
colnames(ts_miR1468_3)<-ts_miR1468_3[1,]
ts_miR1468_3<-ts_miR1468_3[-1,]
tail(ts_miR1468_3, 20)
ts_miR1468_3<-head(ts_miR1468_3, -14)
rownames(ts_miR1468_3)<-NULL
write.csv(ts_miR1468_3, "ExpressionData/targetscan/fFS-TX/ts_miR1468_3_mod.csv")
names(ts_miR1468_3)[2]<-paste("hs_transcript")
list_miR1468_3<-list(hs_ens_geneid, ts_miR1468_3)
merged_list_miR1468_3<-Reduce(function(...) merge(..., by=c("hs_transcript"), all=F), list_miR1468_3)
merged_list_miR1468_3
write.csv(merged_list_miR1468_3, "ExpressionData/targetscan/fFS-TX/DE_targets_DE_miRNAs_miR1468_3.csv")

#miR-361
ts_miR361<-read.csv("ExpressionData/targetscan/fFS-TX/TargetScan7.0__miR-361.predicted_targets.txt", header=FALSE, sep = "\t", stringsAsFactors = FALSE)
ts_miR361_3<-ts_miR361[which(ts_miR361$V11 >= -0.3),]
ts_miR361_3$V2<-sub("\\.\\d$", "", ts_miR361_3$V2)
colnames(ts_miR361_3)<-ts_miR361_3[1,]
ts_miR361_3<-ts_miR361_3[-1,]
tail(ts_miR361_3, 10)
ts_miR361_3<-head(ts_miR361_3, -6)
rownames(ts_miR361_3)<-NULL
write.csv(ts_miR361_3, "ExpressionData/targetscan/fFS-TX/ts_miR361_3_mod.csv")
names(ts_miR361_3)[2]<-paste("hs_transcript")
list_miR361_3<-list(hs_ens_geneid, ts_miR361_3)
merged_list_miR361_3<-Reduce(function(...) merge(..., by=c("hs_transcript"), all=F), list_miR361_3)
merged_list_miR361_3
write.csv(merged_list_miR361_3, "ExpressionData/targetscan/fFS-TX/DE_targets_DE_miRNAs_miR361_3.csv")


##############################################################################################
#The difference between Finnsheep and F1 with control diet
##############################################################################################
cFS_F1 <- results(dds, contrast=c("group", "FS.C", "F1.C"))
cFS_F1_sig<-subset(cFS_F1, padj < 0.05)
nrow(cFS_F1_sig) #2
summary(cFS_F1_sig)
write.csv(cFS_F1_sig, "ExpressionData/mRNA_out/cFS_F1_sig.csv")
rncFS_F1_sig<-rownames(cFS_F1_sig)
gocFS_F1_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id", "go_id"), values=rncFS_F1_sig, mart=ensembl)
write.csv(gocFS_F1_sig, "ExpressionData/mRNA_out/gocFS_F1_sig.csv")
bmcFS_F1_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id",
                                                                                      "external_gene_name", "chromosome_name", "description"), values=rncFS_F1_sig, mart=ensembl)
ortho_bmcFS_F1_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id", "external_gene_name",
                                                                                            "btaurus_homolog_ensembl_gene", "sscrofa_homolog_ensembl_gene", "mmusculus_homolog_ensembl_gene", "hsapiens_homolog_ensembl_gene"), values=rncFS_F1_sig, mart=ensembl)
write.csv(ortho_bmcFS_F1_sig, "ExpressionData/mRNA_out/ortho_bmcFS_F1_sig.csv")
sortbmcFS_F1_sig<-cbind(bmcFS_F1_sig[match(rncFS_F1_sig, bmcFS_F1_sig$ensembl_gene_id),])
bmcFS_F1<-cbind(cFS_F1_sig, bmcFS_F1_sig)
write.csv(bmcFS_F1, "ExpressionData/mRNA_out/bmcFS_F1_sig.csv")

##############################################################################################
#The difference between Finnsheep and F1 with flushing diet
##############################################################################################
fFS_F1 <- results(dds, contrast=c("group", "FS.F", "F1.F"))
fFS_F1_sig<-subset(fFS_F1, padj < 0.05)
nrow(fFS_F1_sig) #47
summary(fFS_F1_sig)
write.csv(fFS_F1_sig, "ExpressionData/mRNA_out/fFS_F1_sig.csv")
rnfFS_F1_sig<-rownames(fFS_F1_sig)
gofFS_F1_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id", "go_id"), values=rnfFS_F1_sig, mart=ensembl)
write.csv(gofFS_F1_sig, "ExpressionData/mRNA_out/gofFS_F1_sig.csv")
bmfFS_F1_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id",
                                                                                      "external_gene_name", "chromosome_name", "description"), values=rnfFS_F1_sig, mart=ensembl)
ortho_bmfFS_F1_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id", "external_gene_name",
                                                                                            "btaurus_homolog_ensembl_gene", "sscrofa_homolog_ensembl_gene", "mmusculus_homolog_ensembl_gene", "hsapiens_homolog_ensembl_gene"), values=rnfFS_F1_sig, mart=ensembl)
write.csv(ortho_bmfFS_F1_sig, "ExpressionData/mRNA_out/ortho_bmfFS_F1_sig.csv")
sortbmfFS_F1_sig<-cbind(bmfFS_F1_sig[match(rnfFS_F1_sig, bmfFS_F1_sig$ensembl_gene_id),])
bmfFS_F1<-cbind(fFS_F1_sig, bmfFS_F1_sig)
write.csv(bmfFS_F1, "ExpressionData/mRNA_out/bmfFS_F1_sig.csv")

################################################################################################
#The difference between Texel and F1 with control diet
################################################################################################
cTX_F1 <- results(dds, contrast=c("group", "TX.C", "F1.C"))
cTX_F1_sig<-subset(cTX_F1, padj < 0.05)
nrow(cTX_F1_sig) #57
summary(cTX_F1_sig)
write.csv(cTX_F1_sig, "ExpressionData/mRNA_out/cTX_F1_sig.csv")
rncTX_F1_sig<-rownames(cTX_F1_sig)
gocTX_F1_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id", "go_id"), values=rncTX_F1_sig, mart=ensembl)
write.csv(gocTX_F1_sig, "ExpressionData/mRNA_out/gocTX_F1_sig.csv")
bmcTX_F1_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id",
                                                                                      "external_gene_name", "chromosome_name", "description"), values=rncTX_F1_sig, mart=ensembl)
ortho_bmcTX_F1_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id", "external_gene_name",
                                                                                            "btaurus_homolog_ensembl_gene", "sscrofa_homolog_ensembl_gene", "mmusculus_homolog_ensembl_gene", "hsapiens_homolog_ensembl_gene"), values=rncTX_F1_sig, mart=ensembl)
write.csv(ortho_bmcTX_F1_sig, "ExpressionData/mRNA_out/ortho_bmcTX_F1_sig.csv")
sortbmcTX_F1_sig<-cbind(bmcTX_F1_sig[match(rncTX_F1_sig, bmcTX_F1_sig$ensembl_gene_id),])
bmcTX_F1<-cbind(cTX_F1_sig, bmcTX_F1_sig)
write.csv(bmcTX_F1, "ExpressionData/mRNA_out/bmcTX_F1_sig.csv")

################################################################################################
#The difference between Texel and F1 with flushing diet
################################################################################################
fTX_F1 <- results(dds, contrast=c("group", "TX.F", "F1.F"))
fTX_F1_sig<-subset(fTX_F1, padj<0.05)
nrow(fTX_F1_sig) #305
summary(fTX_F1_sig)
write.csv(fTX_F1_sig, "ExpressionData/mRNA_out/fTX_F1_sig.csv")
rnfTX_F1_sig<-rownames(fTX_F1_sig)
gofTX_F1_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id", "go_id"), values=rnfTX_F1_sig, mart=ensembl)
write.csv(gofTX_F1_sig, "ExpressionData/mRNA_out/gofTX_F1_sig.csv")
bmfTX_F1_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id",
                                                                                      "external_gene_name", "chromosome_name", "description"), values=rnfTX_F1_sig, mart=ensembl)
ortho_bmfTX_F1_sig<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("ensembl_gene_id", "external_gene_name",
                                                                                            "btaurus_homolog_ensembl_gene", "sscrofa_homolog_ensembl_gene", "mmusculus_homolog_ensembl_gene", "hsapiens_homolog_ensembl_gene"), values=rnfTX_F1_sig, mart=ensembl)
write.csv(ortho_bmfTX_F1_sig, "ExpressionData/mRNA_out/ortho_bmfTX_F1_sig.csv")
sortbmfTX_F1_sig<-cbind(bmfTX_F1_sig[match(rnfTX_F1_sig, bmfTX_F1_sig$ensembl_gene_id),])
bmfTX_F1<-cbind(fTX_F1_sig, bmfTX_F1_sig)
write.csv(bmfTX_F1, "ExpressionData/mRNA_out/bmfTX_F1_sig.csv")


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
write.csv(all_rn, "ExpressionData/mRNA_out/all_rownames.csv")

bm_all_de_entrez<-getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("entrezgene", "ensembl_gene_id",
                                                                                          "ensembl_transcript_id", "ensembl_peptide_id"), values=all_rn, mart=ensembl)
write.csv((bm_all_de_entrez), "ExpressionData/mRNA_out/bm_all_de_entrez.csv")

#It was found out that 
all_cluego<-read.csv("ExpressionData/mRNA_out/cluego_unannotated", header = T)
head(all_cluego)
bm_all_cluego_entrez<-unique(getBM(c("oaries_gene_ensembl"), filters="ensembl_gene_id", attributes=c("entrezgene", "ensembl_gene_id",
                                                                                              "ensembl_transcript_id", "ensembl_peptide_id"), values=all_cluego, mart=ensembl))
bm_all_cluego_entrez
tail(bm_all_cluego_entrez)
write.csv(na.omit(bm_all_cluego_entrez), "ExpressionData/mRNA_out/bm_all_cluego_entrez1.csv")

#Venn diagrams
library(VennDiagram)
pdf("venn_dFS-TX-F1")

vennplot<-venn.diagram(list(df4$rn, df5$rn, df6$rn), NULL, fill=c("red", "green", "yellow"), alpha=c(0.5,0.5, 0.5), cex = 2, cat.fontface=4, category.names=c("FSvsTX", "FSvsF1", "TXvsF1"))
grid.draw(vennplot)
dev.off()

pdf("venn_dFS-TX-F1")

vennplot<-draw.triple.venn(list(df7$rn, df9$rn, df11$rn), NULL, fill=c("red", "green", "yellow"), alpha=c(0.5,0.5, 0.5), cex = 2, cat.fontface=4, category.names=c("cFSvsTX", "cFSvsF1", "cTXvsF1"))
grid.draw(vennplot)
dev.off()

vennplot<-venn.diagram(list(df8$rn, df10$rn, df12$rn), NULL, fill=c("red", "green", "yellow"), alpha=c(0.5,0.5, 0.5), cex = 2, cat.fontface=4, category.names=c("cFSvsTX", "cFSvsF1", "cTXvsF1"))
grid.draw(vennplot)
dev.off()

baseMeanBreed <- sapply( levels(dds$breed), function(lvl) rowMeans( counts(dds,normalized=TRUE)[,dds$breed == lvl] ) )
bm<-as.data.frame(baseMeanBreed)
length(bm$FS[bm$FS >=5]) #16079, total number of genes with basemean value >= 5 in Finnsheep
length(bm$TX[bm$TX >=5]) #16039
length(bm$F1[bm$F1 >=5]) #17345

#Venndiagram with top 500 genes in each breed
bm_fs_50<-rownames(bm.fs)[1:500]
bm.tx<-bm[order(bm$TX, decreasing = TRUE),]
bm_tx_50<-rownames(bm.tx)[1:500]
bm.f1<-bm[order(bm$F1, decreasing = TRUE),]
bm_f1_50<-rownames(bm.f1)[1:500]

vennplot<-venn.diagram(list(bm_fs_50, bm_tx_50, bm_f1_50), NULL, fill=c(1, 2, 3), alpha=c(0.5,0.5, 0.5), resolution = 500, cex = 1.5, cat.fontface="plain", cat.fontfamily="serif", category.names=c("FS", "TX", "F1"))
grid.draw(vennplot)
dev.off()

#bmdFS<-within(bmd, rm(F1, TX)) #optional to extract the subset of data frame that consist ensID for particular breed and base Mean
#bmdF1<-within(bmd, rm(FS,TX))
#bmdTX<-within(bmd, rm(FS,F1))

bmFS<-subset(bm, bm$FS)
a<-subset(bm$FS, bm$FS>5) #16079
nrow(a)
head(a)
length(bmd$F1[bmd$F1>=5]) # 17345

length(bm$TX[bm$TX>=5]) #16039
row.names(bm)[which (bm$FS>=5),]
fs <- which(bm$FS >=5)
nrow(fs)
head(fs)
head(bmd)
head(baseMeanBreed)
nrow(rownames(bmd$FS > 5))
genesFS
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
rld <- rlog(dds)

dds_fFS_TX<-dds[ , dds$group %in% c("FS.F" ,"TX.F")]
#data[!(data$v1 %in% c("b", "d", "e")), ]
rld_dds_fFS_TX<-rlog(dds_fFS_TX)
plotPCA(rld_dds_fFS_TX, intgroup=c("breed", "diet"))
plotPCA(rld, intgroup=c("breed", "diet"))

dds_fTX_F1<-dds[ , dds$group %in% c("TX.F" ,"F1.F")]
#data[!(data$v1 %in% c("b", "d", "e")), ]
rld_dds_fTX_F1<-rlog(dds_fTX_F1)
plotPCA(rld_dds_fTX_F1, intgroup=c("breed", "diet"))

library("DESeq2")
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:20]
nt <- normTransform(dds) # defaults to log2(x+1)
log2.norm.counts <- assay(nt)[select,]
head(log2.norm.counts)
df <- as.data.frame(colData(dds)[,c("breed","diet")])
pheatmap(log2.norm.counts, cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df)
#pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
#         cluster_cols=FALSE, annotation_col=df)

library("genefilter")
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),20)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld)[,c("breed","diet")])
pheatmap(mat, annotation_col=df)

#samplematrix
sampleDists <- dist( t( assay(rld) ) )
sampleDists
sampleDistMatrix <- as.matrix( sampleDists )
#rownames(sampleDistMatrix) <- paste( rld$breed, rld$diet, sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors, annotation_row = df, show_rownames = TRUE, show_colnames = TRUE)

# Principal components analysis
rld_pca <- function (rld_dds_fFS_TX, intgroup = "group", ntop = 500, colors=NULL, legendpos="bottomleft", main="PCA Biplot", textcx=1, ...) {
  require(genefilter)
  require(calibrate)
  require(RColorBrewer)
  rv = rowVars(assay(rld_dds_fFS_TX))
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld_dds_fFS_TX)[select, ]))
  fac = factor(apply(as.data.frame(colData(rld_dds_fFS_TX)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
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
rld_pca(rld_dds_fFS_TX, colors=NULL, intgroup="group", xlim=c(-50, 50))
dev.off()


# Get differential expression results
res <- results(dds)
res
table(res$padj<0.05)
## Order by adjusted p-value
res <- res[order(res$baseMean, decreasing = TRUE), ]
## Merge with normalized count data
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)
## Write results
write.csv(resdata, file="ExpressionData/mRNA_out/diffexpr-results.csv")

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
maplot <- function (res, thresh=0.05, labelsig=FALSE, textcx=1, ...) {
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

plotMA(res, ylim=c(-3, 3))

## Volcano plot with "significant" genes labeled, one example
volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=FALSE, textcx=1, ...) {
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

select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:20]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(counts(dds,normalized=TRUE)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(8,10))

heatmap.2(assay(rld)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(8, 12))
