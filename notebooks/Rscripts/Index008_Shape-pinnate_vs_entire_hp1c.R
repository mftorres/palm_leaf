library(ape)
library(phytools)
library(MCMCglmmRAM)
completeFun <- function(data, desiredCols) {
	completeVec <- complete.cases(data[, desiredCols])
	return(data[completeVec, ])
}
dataall<-read.csv("./palms_alltraits_curated_20220620.csv",quote="",sep="\t",header=TRUE)
posdis42<-read.tree("./Clean_1_42presampled.trees")
data_pre<-completeFun(dataall,c('CHELSA_ai_stand', 'CHELSA_bio1_stand', 'CHELSA_bio4_stand', 'CHELSA_bio15_stand', 'HeightOverCanopy_stand'))
data_pre<-na.omit(object= dataall, cols=c('CHELSA_ai_stand', 'CHELSA_bio1_stand', 'CHELSA_bio4_stand', 'CHELSA_bio15_stand', 'HeightOverCanopy_stand'))
filter<-c("pinnate","entire")
data_hp1c<-subset(data_pre,shape %in% filter)
rownames(data_hp1c) <- data_hp1c$tip_name
data_hp1c$pinnate_binomial<-factor(data_hp1c$pinnate_binomial)
data_hp1c$CHELSA_ai_stand<-as.numeric(data_hp1c$CHELSA_ai_stand)
data_hp1c$CHELSA_bio1_stand<-as.numeric(data_hp1c$CHELSA_bio1_stand)
data_hp1c$CHELSA_bio4_stand<-as.numeric(data_hp1c$CHELSA_bio4_stand)
data_hp1c$CHELSA_bio15_stand<-as.numeric(data_hp1c$CHELSA_bio15_stand)
data_hp1c$HeightOverCanopy_stand<-as.numeric(data_hp1c$HeightOverCanopy_stand)
tree1<-posdis42[[1]]
missingspp<-setdiff(sort(tree1$tip.label),sort(data_hp1c$tip_name))
data_hp1c$animal <- factor(data_hp1c$tip_name)
Nburn <- 9000
Nnitt <- 2000000
Nthin <- 2000
k <- 2
I <- diag(k-1)
J <- matrix(rep(1, (k-1)^2), c(k-1, k-1))
priors<-list(R=list(V=(1/k)*(I+J), fix=1), G=list(G1=list(V=diag(k-1), nu=0.002)))
n_tree=42
hp1c_postdist<-c()
for(tree in 1:n_tree){
	tree2<-drop.tip(posdis42[[tree]],c(missingspp))
	modelhp1c<- MCMCglmm(pinnate_binomial~CHELSA_ai_stand+CHELSA_bio1_stand+CHELSA_bio4_stand+CHELSA_bio15_stand+HeightOverCanopy_stand,
			random = ~animal,
			data = data_hp1c,
			reduced = TRUE,
			pedigree = tree2,
			family = "threshold",
			prior = priors,
			verbose = TRUE,
			burnin = Nburn, nitt = Nnitt, thin = Nthin,
			pr = TRUE, pl = TRUE, saveX = TRUE,  saveZ = TRUE)
	hp1c_postdist<-rbind(hp1c_postdist,modelhp1c$Sol)
	write.table(hp1c_postdist,"./Shape-pinnate_vs_entire_hp1c_postdist.txt",sep="\t")
}
write.table(hp1c_postdist,"./Shape-pinnate_vs_entire_hp1c_postdist.txt",sep="\t")
save.image("./Shape-pinnate_vs_entire_hp1c.Rimage")
