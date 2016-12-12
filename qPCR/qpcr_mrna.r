library(MCMC.qpcr)
mRNA_qpcr<-read.csv("qPCR/mrna_raw_ed_final.csv")
mRNA_qpcr
mRNA_eff<-read.csv("qPCR/mRNA_efficiency_mod_final.csv")
mRNA_eff
mRNA_qs<-cq2counts(data = mRNA_qpcr, genecols = c(4:9), condcols = c(1:3), effic = mRNA_eff)
mRNA_qs$breed.diet<-as.factor(paste(mRNA_qs$breed, mRNA_qs$diet, sep = "."))
naive<-mcmc.qpcr(data=mRNA_qs, fixed = "breed.diet", controls = "GAPDH")
S1<-HPDsummary(model = naive, data = mRNA_qs)  
S2<-HPDsummary(model = naive, data = mRNA_qs, relative = T)

#mRNA_naive<-mcmc.qpcr(fixed = "breed+diet+breed:diet", random = "gene", data = mRNA_qs, controls = "GAPDH", normalize = T)
#mRNA_smm.naive<-HPDsummary(mRNA_naive, mRNA_qs, relative = T)



miRNA_qpcr<-read.csv("qPCR/miRNA CT values.csv")
miRNA_qpcr
miRNA_eff<-read.csv("qPCR/miRNA_efficiency.csv")
miRNA_eff
miRNA_qs<-cq2counts(data=miRNA_qpcr, genecols = c(4:13), condcols=c(1:3), effic = miRNA_eff)
miRNA_qs$breed.diet<-as.factor(paste(miRNA_qs$breed, miRNA_qs$diet, sep = "."))
miRNA_naive<-mcmc.qpcr(data = miRNA_qs, fixed = "breed.diet", controls = "U6")
miRNA_S1<-HPDsummary(model = miRNA_naive, data = miRNA_qs)

slices<-c(469588, 628394)
lbls<-c("Deleterious", "Tolerated")
pie(slices, labels = lbls)

slices1<-c(163346, 163569, 142673)
lbls1<-c("FS", "TX", "F1")
pie(slices1, labels = lbls1)


slices <- c(10, 12,4, 16, 8)
lbls <- c("US", "UK", "Australia", "Germany", "France")
pie(slices, labels = lbls, main="Pie Chart of Countries")

