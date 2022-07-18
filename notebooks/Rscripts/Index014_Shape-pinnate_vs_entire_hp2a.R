library(ape)
library(phytools)
library(MCMCglmmRAM)
library(dplyr)
library(foreach)
library(doParallel)
numCores <- detectCores()
numCores
registerDoParallel(numCores)
completeFun <- function(data, desiredCols) {
	completeVec <- complete.cases(data[, desiredCols])
	return(data[completeVec, ])
}
dataall<-read.csv("./palms_alltraits_curated_20220620.csv",quote="",sep="\t",header=TRUE)
posdis42<-read.tree("./Clean_1_42presampled.trees")
data_pre<-completeFun(dataall,c('CHELSA_ai_stand', 'CHELSA_bio1_stand', 'CHELSA_bio4_stand', 'CHELSA_bio15_stand', 'Max_Rachis_Length_m_stand', 'HeightOverCanopy_stand'))
data_hp2a<-filter(data_pre, pinnate_binomial == "True" | entire_binomial == "True")
rownames(data_hp2a) <- data_hp2a$tip_name
data_hp2a$pinnate_binomial<-factor(data_hp2a$pinnate_binomial)
data_hp2a$CHELSA_ai_stand<-as.numeric(data_hp2a$CHELSA_ai_stand)
data_hp2a$CHELSA_bio1_stand<-as.numeric(data_hp2a$CHELSA_bio1_stand)
data_hp2a$CHELSA_bio4_stand<-as.numeric(data_hp2a$CHELSA_bio4_stand)
data_hp2a$CHELSA_bio15_stand<-as.numeric(data_hp2a$CHELSA_bio15_stand)
data_hp2a$Max_Rachis_Length_m_stand<-as.numeric(data_hp2a$Max_Rachis_Length_m_stand)
data_hp2a$HeightOverCanopy_stand<-as.numeric(data_hp2a$HeightOverCanopy_stand)
tree1<-posdis42[[1]]
missingspp<-setdiff(sort(tree1$tip.label),sort(data_hp2a$tip_name))
data_hp2a$animal <- factor(data_hp2a$tip_name)
Nburn <- 500
Nnitt <- 10000000
Nthin <- 5000
k <- 2
I <- diag(k-1)
J <- matrix(rep(1, (k-1)^2), c(k-1, k-1))
priors<-list(R=list(V=(1/k)*(I+J), fix=1), G=list(G1=list(V=diag(k-1), nu=0.002)))
n_tree=42
packages=c("ape","phytools","MCMCglmmRAM","dplyr")
hp2a_postdist <-c()
hp2a_postdist <- foreach(i=1:n_tree, .combine=rbind, .packages=packages) %dopar% {
	tree2<-drop.tip(posdis42[[i]],c(missingspp))
	modelhp2a<- MCMCglmm(pinnate_binomial~CHELSA_ai_stand+CHELSA_bio1_stand+CHELSA_bio4_stand+CHELSA_bio15_stand+Max_Rachis_Length_m_stand+HeightOverCanopy_stand,
			random = ~animal,
			data = data_hp2a,
			reduced = TRUE,
			pedigree = tree2,
			family = "threshold",
			prior = priors,
			verbose = TRUE,
			burnin = Nburn, nitt = Nnitt, thin = Nthin,
			pr = TRUE, pl = TRUE, saveX = TRUE,  saveZ = TRUE)
	hp2a_postdist<-rbind(hp2a_postdist,modelhp2a$Sol)
	write.table(hp2a_postdist,"./Shape-pinnate_vs_entire_hp2a_postdist-1.txt",sep="\t")
}
write.table(hp2a_postdist,"./Shape-pinnate_vs_entire_hp2a_postdist-1.txt",sep="\t")
save.image("./Shape-pinnate_vs_entire_hp2a-1.Rimage")
