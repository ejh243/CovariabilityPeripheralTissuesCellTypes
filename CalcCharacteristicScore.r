## calculate a variance characteristic score
## between blood cell types compare each cell types
## adjust for mean differences in DNAm first

library(car)
library(parallel)
cl <- makeCluster(16)
clusterExport(cl,"leveneTest")

adjustCT<-function(row, sampleType){
	return(residuals(lm(row ~ sampleType)))
}

varCharScore<-function(row, indVar){
	return(leveneTest(row ~ indVar, center = "median")[1,3])
}


setwd("")
load("CellSortedAsthmaERisk_Normalised.rdat")

pheno$Sample.Type<-as.factor(pheno$Sample.Type)


## filter to autosomal
epicManifest<-read.csv("MethylationEPIC_v-1-0_B2.csv", skip = 7)
epicManifest<-epicManifest[match(rownames(betas), epicManifest$Name),]
betas<-betas[which(epicManifest$CHR != "X" & epicManifest$CHR != "Y"),]

## subset to blood cell types only
include<-c("B-cells","CD4 T-cells","CD8 T-cells","Granulocytes","Monocytes")
betas<-betas[,which(pheno$Sample.Type %in% include)]
pheno<-pheno[match(colnames(betas), pheno$Basename),]

## adjust for differences in means between cell types
betasAdj<-t(apply(betas, 1, adjustCT, pheno$Sample.Type))

pheno$Sample.Type<-as.factor(as.character(pheno$Sample.Type))

out<-matrix(data = NA, nrow = nrow(betasAdj), ncol = 5)
colnames(out)<-levels(pheno$Sample.Type)
rownames(out)<-rownames(betas)
for(each in levels(pheno$Sample.Type)){
	indVar<-rep(0, ncol(betasAdj))
	indVar[which(pheno$Sample.Type == each)]<-1
	indVar<-as.factor(indVar)
	tmp<-parApply(cl, betasAdj, 1, varCharScore, indVar)
	out[,each]<-tmp
	out[,each]<--log10(out[,each])
	## as non-parametric test can test for spcific direction so instead will use sign to indicate direction of effect
	## calculate difference in sd to determine sign
	coef<-sign(parApply(cl, betasAdj[,which(indVar == 1)], 1, sd)-parApply(cl,betasAdj[,which(indVar == 0)], 1, sd))
	out[,each]<-out[,each]*coef
}

write.csv(out, "CellTypeComparisons/VarianceCharacteristicScores.csv")
