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
data_pre<-completeFun(dataall,c('CHELSA_ai_stand', 'CHELSA_bio1_stand', 'CHELSA_bio4_stand', 'CHELSA_bio15_stand', 'StemHeightBladeLength_stand'))

data_hp1b<-filter(data_pre, entire_binomial == "True" | dissection == 1)
rownames(data_hp1b) <- data_hp1b$tip_name
data_hp1b$dissection<-factor(data_hp1b$dissection)
data_hp1b$CHELSA_ai_stand<-as.numeric(data_hp1b$CHELSA_ai_stand)
data_hp1b$CHELSA_bio1_stand<-as.numeric(data_hp1b$CHELSA_bio1_stand)
data_hp1b$CHELSA_bio12_stand<-as.numeric(data_hp1b$CHELSA_bio12_stand)
data_hp1b$CHELSA_bio15_stand<-as.numeric(data_hp1b$CHELSA_bio15_stand)
data_hp1b$StemHeightBladeLength_stand<-as.numeric(data_hp1b$StemHeightBladeLength_stand)
tree1<-posdis42[[1]]
missingspp<-setdiff(sort(tree1$tip.label),sort(data_hp1b$tip_name))
data_hp1b$animal <- factor(data_hp1b$tip_name)
Nburn <- 9000
Nnitt <- 20000000
Nthin <- 8000
k <- 2
I <- diag(k-1)
J <- matrix(rep(1, (k-1)^2), c(k-1, k-1))
priors<-list(R=list(V=(1/k)*(I+J), fix=1), G=list(G1=list(V=diag(k-1), nu=0.002)))
n_tree=42
packages=c("ape","phytools","MCMCglmmRAM","dplyr")
hp1b_postdist<-c()
hp1b_postdist <- foreach(i=1:n_tree, .combine=rbind, .packages=packages) %dopar% {
	tree2<-drop.tip(posdis42[[i]],c(missingspp))
	modelhp1b<- MCMCglmm(dissection~CHELSA_ai_stand+CHELSA_bio1_stand+CHELSA_bio12_stand+CHELSA_bio15_stand+StemHeightBladeLength_stand,
			random = ~animal,
			data = data_hp1b,
			reduced = TRUE,
			pedigree = tree2,
			family = "threshold",
			prior = priors,
			verbose = TRUE,
			burnin = Nburn, nitt = Nnitt, thin = Nthin,
			pr = TRUE, pl = TRUE, saveX = TRUE,  saveZ = TRUE)
	hp1b_postdist<-rbind(hp1b_postdist,modelhp1b$Sol)
}
write.table(hp1b_postdist,"./Shape-entire_vs_dissected_hp1b_postdist-2.txt",sep="\t")
save.image("./Shape-entire_vs_dissected_hp1b-2.Rimage")
