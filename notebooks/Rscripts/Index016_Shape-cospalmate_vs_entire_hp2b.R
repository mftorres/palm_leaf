library(ape)
library(phytools)
library(MCMCglmmRAM)
library(dplyr)
completeFun <- function(data, desiredCols) {
	completeVec <- complete.cases(data[, desiredCols])
	return(data[completeVec, ])
}
dataall<-read.csv("./palms_alltraits_curated_20220620.csv",quote="",sep="\t",header=TRUE)
posdis42<-read.tree("./Clean_1_42presampled.trees")
data_pre<-completeFun(dataall,c('CHELSA_ai_stand', 'CHELSA_bio1_stand', 'CHELSA_bio12_stand', 'CHELSA_bio15_stand', 'Max_Rachis_Length_m_stand', 'HeightOverCanopy_stand'))
data_hp2b<-filter(data_pre, cospalmate_binomial == "True" | entire_binomial == "True")
rownames(data_hp2b) <- data_hp2b$tip_name
data_hp2b$cospalmate_binomial<-factor(data_hp2b$cospalmate_binomial)
data_hp2b$CHELSA_ai_stand<-as.numeric(data_hp2b$CHELSA_ai_stand)
data_hp2b$CHELSA_bio1_stand<-as.numeric(data_hp2b$CHELSA_bio1_stand)
data_hp2b$CHELSA_bio12_stand<-as.numeric(data_hp2b$CHELSA_bio12_stand)
data_hp2b$CHELSA_bio15_stand<-as.numeric(data_hp2b$CHELSA_bio15_stand)
data_hp2b$Max_Rachis_Length_m_stand<-as.numeric(data_hp2b$Max_Rachis_Length_m_stand)
data_hp2b$HeightOverCanopy_stand<-as.numeric(data_hp2b$HeightOverCanopy_stand)
tree1<-posdis42[[1]]
missingspp<-setdiff(sort(tree1$tip.label),sort(data_hp2b$tip_name))
data_hp2b$animal <- factor(data_hp2b$tip_name)
Nburn <- 500
Nnitt <- 10000000
Nthin <- 5000
k <- 2
I <- diag(k-1)
J <- matrix(rep(1, (k-1)^2), c(k-1, k-1))
priors<-list(R=list(V=(1/k)*(I+J), fix=1), G=list(G1=list(V=diag(k-1), nu=0.002)))
n_tree=42
hp2b_postdist<-c()
for(tree in 1:n_tree){
	tree2<-drop.tip(posdis42[[tree]],c(missingspp))
	modelhp2b<- MCMCglmm(cospalmate_binomial~CHELSA_ai_stand+CHELSA_bio1_stand+CHELSA_bio12_stand+CHELSA_bio15_stand+Max_Rachis_Length_m_stand+HeightOverCanopy_stand,
			random = ~animal,
			data = data_hp2b,
			reduced = TRUE,
			pedigree = tree2,
			family = "threshold",
			prior = priors,
			verbose = TRUE,
			burnin = Nburn, nitt = Nnitt, thin = Nthin,
			pr = TRUE, pl = TRUE, saveX = TRUE,  saveZ = TRUE)
	hp2b_postdist<-rbind(hp2b_postdist,modelhp2b$Sol)
	write.table(hp2b_postdist,"./Shape-cospalmate_vs_entire_hp2b_postdist.txt",sep="\t")
}
write.table(hp2b_postdist,"./Shape-cospalmate_vs_entire_hp2b_postdist.txt",sep="\t")
save.image("./Shape-cospalmate_vs_entire_hp2b.Rimage")
