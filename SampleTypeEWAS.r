## perform anova across all samples to identify sites with different means between cell types and tissues
## perform levene's test across all samples to identify sites with different variances between cell types and tissues

setwd("/mnt/data1/Eilis/Projects/Asthma/")
load("QC/CombinedQC_WithMarkdown/CellSortedAsthmaERisk_Normalised.rdat")

pheno$Sample.Type<-as.factor(pheno$Sample.Type)
pheno$Sample.Type<-relevel(pheno$Sample.Type, ref = "whole blood") ## take whole blood as reference



out<-matrix(data = NA, nrow = nrow(betas), ncol = 25)
colnames(out)<-c("ANOVA:P", paste("P", levels(pheno$Sample.Type)[-1]), paste("Mean", levels(pheno$Sample.Type)), "LeveneP", paste("SD", levels(pheno$Sample.Type)))
for(i in 1:nrow(betas)){
	model<-lm(betas[i,] ~ pheno$Sample.Type)
	out[i,1]<-anova(model)[1,5]
	out[i,2:8]<-summary(model)$coefficients[-1,4]
	out[i,9:16]<-c(aggregate(betas[i,], by = list(pheno$Sample.Type), FUN = "mean")[,2])
	out[i,17]<-levene.test(betas[i,], pheno$Sample.Type)$p.value
	out[i,18:25]<-c(aggregate(betas[i,], by = list(pheno$Sample.Type), FUN = "sd")[,2])

}

write.csv(out, "CellTypeComparisons/AnovaCompareCellTypes.csv")
