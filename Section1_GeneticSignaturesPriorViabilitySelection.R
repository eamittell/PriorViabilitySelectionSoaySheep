#### Script for models from section 1 of the manuscript: Quantitative genetic analysis of early life viability resolves a case of the paradox of stasis in body size traits of Soay sheep
### Lizy Mittell
## Section 1: Genetic signatures of prior viability selection
# Cleaned 9th May 2025

# Libraries
library(MCMCglmm)
library(QGglmm)

# File path to data -- edit as appropriate
file_path_data <- "/Users/umer/Documents/Work/Evol_Climate_Change/Analyses/SoaySheep/ProcessedData/"
# File path to output -- edit as appropriate
file_path_output <- ""

# Load pedigree
pedigree <- read.csv(paste(file_path_data,"SheepPedigreeRandomIds_May2025.csv",sep=""))
# Format
pedigree <- pedigree[,-1]
colnames(pedigree)[1] <- "animal"
pedigree$animal <- as.factor(pedigree$animal)
pedigree$dam <- as.factor(pedigree$dam)
pedigree$sire <- as.factor(pedigree$sire)
str(pedigree)

########################
## Genetic signatures ##
########################
############################
## Adult August body mass ##
############################
# Load model data
data <- read.csv(paste(file_path_data,"CovuDataABMFYS_120525_RandomIDs.csv",sep=""))

# Format
data$animal <- as.factor(data$animal)
data$trait <- as.factor(data$trait)
data$CapAgeYears <- as.factor(data$CapAgeYears)
str(data)
summary(data)
summary(data[which(data$trait=="r"),])
summary(data[which(data$trait=="s"),])

# Prune pedigree
pruned <- prunePed(pedigree,keep=data$id)

# Run model
prior1 <- list(R = list(R1=list(V=diag(2), nu=2.002, covu=TRUE, fix=2), R2=list(V=1, nu=0.002)), G=list(G1=list(V=diag(2),nu=2.002), G2=list(V=diag(2),nu=2.002),G3=list(V=diag(1),nu=2.002,alpha.mu=0,alpha.V=1000)))
modelABM <- MCMCglmm(y ~ trait - 1 + at.level(trait,"r"):CapAgeYears + trait:Sex + trait:Twin + trait:PopDensMC + at.level(trait,"r"):jdateMC + at.level(trait,"s"):MumAgeMC, random=~us(trait):animal + us(trait):BirthYear + us(at.level(trait,"s")):MumID + idh(at.level(trait,"r")):id, rcov=~idh(at.level(trait,"s")):id + idh(at.level(trait,"r")):id:time, family = NULL, pedigree=pruned, data=data, prior=prior1, nitt=1505000, thin=1000, burnin=5000)
# Save
save(modelABM,file=paste(file_path_output,"ABMcovuFYS_120525.RData",sep="")) # Date as appropriate

# Process results
load(paste(file_path_output,"ABMcovuFYS_120525.RData",sep="")) # Date as appropriate
model <- modelABM
# Number of mcmc samples
nMCMC <- 1500

# Latent scale to data scale across the full posterior distribution
DataScaleRes <- NULL
for(i in 1:nMCMC){
	# Extraction of the fitted values
	# Survival
	mu.s.Fixed <- mean(model$Sol[i,2])*model$X[,2] + mean(model$Sol[i,18])*model$X[,18] + mean(model$Sol[i,20])*model$X[,20] + mean(model$Sol[i,22])*model$X[,22] + mean(model$Sol[i,24])*model$X[,24]
	# Repeat measure trait
	mu.r.Fixed <- mean(model$Sol[i,1])*model$X[,1] + mean(model$Sol[i,3])*model$X[,3] + mean(model$Sol[i,4])*model$X[,4] + mean(model$Sol[i,5])*model$X[,5] + mean(model$Sol[i,6])*model$X[,6] + mean(model$Sol[i,7])*model$X[,7] + mean(model$Sol[i,8])*model$X[,8] + mean(model$Sol[i,9])*model$X[,9] + mean(model$Sol[i,10])*model$X[,10] + mean(model$Sol[i,11])*model$X[,11] + mean(model$Sol[i,12])*model$X[,12] + mean(model$Sol[i,13])*model$X[,13] + mean(model$Sol[i,14])*model$X[,14] + mean(model$Sol[i,15])*model$X[,15] + mean(model$Sol[i,16])*model$X[,16] + mean(model$Sol[i,17])*model$X[,17] + mean(model$Sol[i,19])*model$X[,19] + mean(model$Sol[i,21])*model$X[,21] + mean(model$Sol[i,23])*model$X[,23]
	
	# Predicted values for transformations
	pred <- cbind(mu.s.Fixed,mu.r.Fixed)

	# Genetic effect
	Gpred <- matrix(c(mean(model$VCV[i,"traits:traits.animal"]), mean(model$VCV[i,"traits:traitr.animal"]), mean(model$VCV[i,"traits:traitr.animal"]), mean(model$VCV[i,"traitr:traitr.animal"])), ncol = 2)
	# Birth year/Cohort effect
	BYpred <- matrix(c(mean(model$VCV[i,"traits:traits.BirthYear"]), mean(model$VCV[i,"traits:traitr.BirthYear"]), mean(model$VCV[i,"traits:traitr.BirthYear"]), mean(model$VCV[i,"traitr:traitr.BirthYear"])), ncol = 2)
	# Maternal effect -- survival only
	Matpred <- matrix(c(mean(model$VCV[i,'at.level(trait, "s"):at.level(trait, "s").MumID']), 0, 0, 0), ncol = 2)
	# Residual/ID level
	RIpred <- matrix(c(mean(model$VCV[i,'at.level(trait, "s").id:at.level(trait, "s").id']), mean(model$VCV[i,'at.level(trait, "s").id:at.level(trait, "r").id']), mean(model$VCV[i,'at.level(trait, "s").id:at.level(trait, "r").id']), mean(model$VCV[i,'at.level(trait, "r").id:at.level(trait, "r").id'])), ncol = 2)
	# Residual effect and measurement error -- trait
	Rpred <- matrix(c(0,0,0, mean(model$VCV[i,'at.level(trait, "r").id:time'])), ncol = 2)
	Ppred <- Gpred + BYpred + Matpred + RIpred + Rpred
	
	# Transform from the latent scale to the data scale
	obs <- QGmvparams(predict=pred,vcv.G=Gpred,vcv.P=Ppred,model=c("binom1.probit","Gaussian"))
	
	# Results to save
	CovaSZ <- obs$vcv.G.obs[1,2] # Genetic covariance
	Saz1 <- obs$vcv.G.obs[1,2]/obs$mean.obs[1] ## Prior viability genetic selection differential trait
	
	res <- c(CovaSZ,Saz1,i)
	
	DataScaleRes <- rbind(DataScaleRes,res)	
}
DataScaleRes <- as.data.frame(DataScaleRes)
names(DataScaleRes) <- c("CovaSZ","Saz1","i")

# Save
save(DataScaleRes, file=paste(file_path_output,"DataScaleRes_ABMcovuFYS_120525.RData",sep="")) # Date as appropriate

#############################
## Adult metatarsal length ##
#############################
# Load model data
data <- read.csv(paste(file_path_data,"CovuDataAMTLFYS_120525_RandomIDs.csv",sep=""))
data$animal <- as.factor(data$animal)
data$trait <- as.factor(data$trait)
data$CapAgeYears <- as.factor(data$CapAgeYears)
str(data)
summary(data)
summary(data[which(data$trait=="r"),])
summary(data[which(data$trait=="s"),])

# Pruned pedigree and age categorical
pruned <- prunePed(pedigree,keep=data$id)

# Bivariate models of first year survival with later-life metatarsal length
prior1 <- list(R = list(R1=list(V=diag(2), nu=2.002, covu=TRUE, fix=2), R2=list(V=1, nu=0.002)), G=list(G1=list(V=diag(2),nu=2.002), G2=list(V=diag(2),nu=2.002), G3=list(V=diag(1),nu=2.002,alpha.mu=0,alpha.V=1000)))
modelAMTL <- MCMCglmm(y ~ trait - 1 + at.level(trait,"r"):CapAgeYears + trait:Sex + trait:Twin + trait:PopDensMC + at.level(trait,"r"):jdateMC + at.level(trait,"s"):MumAgeMC, random=~us(trait):animal + us(trait):BirthYear + us(at.level(trait,"s")):MumID + idh(at.level(trait,"r")):id, rcov=~idh(at.level(trait,"s")):id + idh(at.level(trait,"r")):id:time, family = NULL, pedigree=pedigree, data=data, prior=prior1,nitt=1505000, thin=1000, burnin=5000)
save(modelAMTL,file=paste(file_path_output,"AMTLcovuFYS_120525.RData",sep="")) # Date as appropriate

# Process results
load(paste(file_path_output,file="AMTLcovuFYS_120525.RData",sep="")) # Date as appropriate
model <- modelAMTL
# Number of mcmc samples
nMCMC <- 1500

# Latent scale to data scale across the full posterior distribution
DataScaleRes <- NULL
for(i in 1:nMCMC){
	# Extraction of the fitted values
	# Survival
	mu.s.Fixed <- mean(model$Sol[i,2])*model$X[,2] + mean(model$Sol[i,18])*model$X[,18] + mean(model$Sol[i,20])*model$X[,20] + mean(model$Sol[i,22])*model$X[,22] + mean(model$Sol[i,24])*model$X[,24]
	# Repeat measure trait
	mu.r.Fixed <- mean(model$Sol[i,1])*model$X[,1] + mean(model$Sol[i,3])*model$X[,3] + mean(model$Sol[i,4])*model$X[,4] + mean(model$Sol[i,5])*model$X[,5] + mean(model$Sol[i,6])*model$X[,6] + mean(model$Sol[i,7])*model$X[,7] + mean(model$Sol[i,8])*model$X[,8] + mean(model$Sol[i,9])*model$X[,9] + mean(model$Sol[i,10])*model$X[,10] + mean(model$Sol[i,11])*model$X[,11] + mean(model$Sol[i,12])*model$X[,12] + mean(model$Sol[i,13])*model$X[,13] + mean(model$Sol[i,14])*model$X[,14] + mean(model$Sol[i,15])*model$X[,15] + mean(model$Sol[i,16])*model$X[,16] + mean(model$Sol[i,17])*model$X[,17] + mean(model$Sol[i,19])*model$X[,19] + mean(model$Sol[i,21])*model$X[,21] + mean(model$Sol[i,23])*model$X[,23]
	
	# Predicted values for transformations
	pred <- cbind(mu.s.Fixed,mu.r.Fixed)

	# Genetic effect
	Gpred <- matrix(c(mean(model$VCV[i,"traits:traits.animal"]), mean(model$VCV[i,"traits:traitr.animal"]), mean(model$VCV[i,"traits:traitr.animal"]), mean(model$VCV[i,"traitr:traitr.animal"])), ncol = 2)
	# Birth year/Cohort effect
	BYpred <- matrix(c(mean(model$VCV[i,"traits:traits.BirthYear"]), mean(model$VCV[i,"traits:traitr.BirthYear"]), mean(model$VCV[i,"traits:traitr.BirthYear"]), mean(model$VCV[i,"traitr:traitr.BirthYear"])), ncol = 2)
	# Maternal effect -- survival only
	Matpred <- matrix(c(mean(model$VCV[i,'at.level(trait, "s"):at.level(trait, "s").MumID']), 0, 0, 0), ncol = 2)
	# Residual/ID level
	RIpred <- matrix(c(mean(model$VCV[i,'at.level(trait, "s").id:at.level(trait, "s").id']), mean(model$VCV[i,'at.level(trait, "s").id:at.level(trait, "r").id']), mean(model$VCV[i,'at.level(trait, "s").id:at.level(trait, "r").id']), mean(model$VCV[i,'at.level(trait, "r").id:at.level(trait, "r").id'])), ncol = 2)
	# Residual effect and measurement error -- trait
	Rpred <- matrix(c(0,0,0, mean(model$VCV[i,'at.level(trait, "r").id:time'])), ncol = 2)
	Ppred <- Gpred + BYpred + Matpred + RIpred + Rpred
#	Ppred <- Gpred + BYpred + RIpred + Rpred
	
	# Transform from the latent scale to the data scale
	obs <- QGmvparams(predict=pred,vcv.G=Gpred,vcv.P=Ppred,model=c("binom1.probit","Gaussian"))
	
	# Results to save
	CovaSZ <- obs$vcv.G.obs[1,2] # Genetic covariance
	Saz1 <- obs$vcv.G.obs[1,2]/obs$mean.obs[1] ## Prior viability genetic selection differential trait
	
	res <- c(CovaSZ,Saz1,i)
	
	DataScaleRes <- rbind(DataScaleRes,res)	
}
DataScaleRes <- as.data.frame(DataScaleRes)
names(DataScaleRes) <- c("CovaSZ","Saz1","i")

save(DataScaleRes, file=paste(file_path_output,"DataScaleRes_AMTLcovuFYS_120525.RData"),sep="")) # Date as appropriate

######################
## Male horn length ##
######################
# Load model data
data <- read.csv(paste(file_path_data,"CovuDataMHLFYS_120525_RandomIDs.csv",sep=""))

# Format
data$animal <- as.factor(data$animal)
data$trait <- as.factor(data$trait)
data$CapAgeYears <- as.factor(data$CapAgeYears)
str(data)
summary(data)
summary(data[which(data$trait=="r"),])
summary(data[which(data$trait=="s"),])

# Pruned pedigree and age categorical
pruned <- prunePed(pedigree,keep=data$id)

# Bivariate models of male first year survival with male horn length
prior1 <- list(R = list(R1=list(V=diag(2), nu=2.002, covu=TRUE, fix=2), R2=list(V=1, nu=0.002)), G=list(G1=list(V=diag(2),nu=2.002), G2=list(V=diag(2),nu=2.002), G3=list(V=diag(1),nu=2.002,alpha.mu=0,alpha.V=1000)))
modelMHL_malesonly <- MCMCglmm(y ~ trait - 1 + at.level(trait,"r"):CapAgeYears + trait:Twin + at.level(trait,"s"):PopDensMC + at.level(trait,"s"):MumAgeMC, random=~us(trait):animal + us(trait):BirthYear + us(at.level(trait,"s")):MumID + idh(at.level(trait,"r")):id, rcov=~idh(at.level(trait,"s")):id + idh(at.level(trait,"r")):id:time, family = NULL, pedigree=pruned, data=data, prior=prior1,nitt=2506000,thin=1000,burnin=6000)
save(modelMHL_malesonly,file=paste(file_path_output,"MHLcovuFYS_malesonly_120525.RData",sep="")) # Date as appropriate

# Process results
load(paste(file_path_output,"MHLcovuFYS_malesonly_120525.RData",sep="")) # Date as appropriate
model <- modelMHL_malesonly
# Number of mcmc samples
nMCMC <- 2500

# Latent scale to data scale across the full posterior distribution
DataScaleRes <- NULL
for(i in 1:nMCMC){
	# Extraction of the fitted values
	# Survival
	mu.s.Fixed <- mean(model$Sol[i,2])*model$X[,2] + mean(model$Sol[i,13])*model$X[,13] + mean(model$Sol[i,14])*model$X[,14] + mean(model$Sol[i,15])*model$X[,15]
	# Repeat measure trait
	mu.r.Fixed <- mean(model$Sol[i,1])*model$X[,1] + mean(model$Sol[i,3])*model$X[,3] + mean(model$Sol[i,4])*model$X[,4] + mean(model$Sol[i,5])*model$X[,5] + mean(model$Sol[i,6])*model$X[,6] + mean(model$Sol[i,7])*model$X[,7] + mean(model$Sol[i,8])*model$X[,8] + mean(model$Sol[i,9])*model$X[,9] + mean(model$Sol[i,10])*model$X[,10] + mean(model$Sol[i,11])*model$X[,11] + mean(model$Sol[i,12])*model$X[,12]
		
	# Predicted values for transformations
	pred <- cbind(mu.s.Fixed,mu.r.Fixed)

	# Genetic effect
	Gpred <- matrix(c(mean(model$VCV[i,"traits:traits.animal"]), mean(model$VCV[i,"traits:traitr.animal"]), mean(model$VCV[i,"traits:traitr.animal"]), mean(model$VCV[i,"traitr:traitr.animal"])), ncol = 2)
	# Birth year/Cohort effect
	BYpred <- matrix(c(mean(model$VCV[i,"traits:traits.BirthYear"]), mean(model$VCV[i,"traits:traitr.BirthYear"]), mean(model$VCV[i,"traits:traitr.BirthYear"]), mean(model$VCV[i,"traitr:traitr.BirthYear"])), ncol = 2)
	# Maternal effect -- survival only
	Matpred <- matrix(c(mean(model$VCV[i,'at.level(trait, "s"):at.level(trait, "s").MumID']), 0, 0, 0), ncol = 2)
	# Residual/ID level
	RIpred <- matrix(c(mean(model$VCV[i,'at.level(trait, "s").id:at.level(trait, "s").id']), mean(model$VCV[i,'at.level(trait, "s").id:at.level(trait, "r").id']), mean(model$VCV[i,'at.level(trait, "s").id:at.level(trait, "r").id']), mean(model$VCV[i,'at.level(trait, "r").id:at.level(trait, "r").id'])), ncol = 2)
	# Residual effect and measurement error -- trait
	Rpred <- matrix(c(0,0,0, mean(model$VCV[i,'at.level(trait, "r").id:time'])), ncol = 2)
	Ppred <- Gpred + BYpred + Matpred + RIpred + Rpred
#	Ppred <- Gpred + BYpred + RIpred + Rpred
	
	# Transform from the latent scale to the data scale
	obs <- QGmvparams(predict=pred,vcv.G=Gpred,vcv.P=Ppred,model=c("binom1.probit","Gaussian"))
	
	# Results to save
	CovaSZ <- obs$vcv.G.obs[1,2] # Genetic covariance
	Saz1 <- obs$vcv.G.obs[1,2]/obs$mean.obs[1] ## Prior viability genetic selection differential trait
	
	res <- c(CovaSZ,Saz1,i)
	
	DataScaleRes <- rbind(DataScaleRes,res)	
}
DataScaleRes <- as.data.frame(DataScaleRes)
names(DataScaleRes) <- c("CovaSZ","Saz1","i")

save(DataScaleRes, file=paste(file_path_output,"DataScaleResMHLcovuFYS_malesonly_120525.RData",sep="")) # Date as appropriate

###########################
## Scrotal circumference ##
###########################
# Load model data
data <- read.csv(paste(file_path_data,"CovuDataMSCFYS_120525_RandomIDs.csv",sep=""))

# Format
data$animal <- as.factor(data$animal)
data$trait <- as.factor(data$trait)
data$CapAgeYears <- as.factor(data$CapAgeYears)
str(data)
summary(data)
summary(data[which(data$trait=="r"),])
summary(data[which(data$trait=="s"),])

# Pruned pedigree and age categorical
pruned <- prunePed(pedigree,keep=data$id)

# Bivariate models of male first year survival with scrotal circumference
prior1 <- list(R = list(R1=list(V=diag(2), nu=4.002, covu=TRUE, fix=2), R2=list(V=1, nu=2.002)), G=list(G1=list(V=diag(2),nu=4.002,alpha.mu=c(0,0),alpha.V=diag(c(1000,100))), G2=list(V=diag(2),nu=4.002,alpha.mu=c(0,0),alpha.V=diag(c(1000,100))), G3=list(V=diag(1),nu=2.002,alpha.mu=0,alpha.V=1000)))
modelMSC_malesonly <- MCMCglmm(y ~ trait - 1 + at.level(trait,"r"):CapAgeYears + trait:Twin + trait:PopDensMC + at.level(trait,"r"):jdateMC + at.level(trait,"s"):MumAgeMC, random=~us(trait):animal + us(trait):BirthYear + us(at.level(trait,"s")):MumID + idh(at.level(trait,"r")):id, rcov=~idh(at.level(trait,"s")):id + idh(at.level(trait,"r")):id:time, family = NULL, pedigree=pruned, data=data, prior=prior1,nitt=10500000,thin=500,burnin=20000)
save(modelMSC_malesonly, file=paste(file_path_output,"MSCcovuFYS_malesonly_120525.RData",sep="")) # Date as appropriate

# Process results
load(paste(file_path_output,"MSCcovuFYS_malesonly_120525.RData",sep="")) # Date as appropriate
model <- modelMSC_malesonly
# Number of mcmc samples
nMCMC <- 20960

# Latent scale to data scale across the full posterior distribution
DataScaleRes <- NULL
for(i in 1:nMCMC){
	# Extraction of the fitted values
	# Survival
	mu.s.Fixed <- mean(model$Sol[i,2])*model$X[,2] + mean(model$Sol[i,14])*model$X[,14] + mean(model$Sol[i,16])*model$X[,16] + mean(model$Sol[i,18])*model$X[,18]
	# Repeat measure trait
	mu.r.Fixed <- mean(model$Sol[i,1])*model$X[,1] + mean(model$Sol[i,3])*model$X[,3] + mean(model$Sol[i,4])*model$X[,4] + mean(model$Sol[i,5])*model$X[,5] + mean(model$Sol[i,6])*model$X[,6] + mean(model$Sol[i,7])*model$X[,7] + mean(model$Sol[i,8])*model$X[,8] + mean(model$Sol[i,9])*model$X[,9] + mean(model$Sol[i,10])*model$X[,10] + mean(model$Sol[i,11])*model$X[,11] + mean(model$Sol[i,12])*model$X[,12] + mean(model$Sol[i,13])*model$X[,13] + mean(model$Sol[i,15])*model$X[,15] + mean(model$Sol[i,17])*model$X[,17]
		
	# Predicted values for transformations
	pred <- cbind(mu.s.Fixed,mu.r.Fixed)

	# Genetic effect
	Gpred <- matrix(c(mean(model$VCV[i,"traits:traits.animal"]), mean(model$VCV[i,"traits:traitr.animal"]), mean(model$VCV[i,"traits:traitr.animal"]), mean(model$VCV[i,"traitr:traitr.animal"])), ncol = 2)
	# Birth year/Cohort effect
	BYpred <- matrix(c(mean(model$VCV[i,"traits:traits.BirthYear"]), mean(model$VCV[i,"traits:traitr.BirthYear"]), mean(model$VCV[i,"traits:traitr.BirthYear"]), mean(model$VCV[i,"traitr:traitr.BirthYear"])), ncol = 2)
	# Maternal effect -- survival only
	Matpred <- matrix(c(mean(model$VCV[i,'at.level(trait, "s"):at.level(trait, "s").MumID']), 0, 0, 0), ncol = 2)
	# Residual/ID level
	RIpred <- matrix(c(mean(model$VCV[i,'at.level(trait, "s").id:at.level(trait, "s").id']), mean(model$VCV[i,'at.level(trait, "s").id:at.level(trait, "r").id']), mean(model$VCV[i,'at.level(trait, "s").id:at.level(trait, "r").id']), mean(model$VCV[i,'at.level(trait, "r").id:at.level(trait, "r").id'])), ncol = 2)
	# Residual effect and measurement error -- trait
	Rpred <- matrix(c(0,0,0, mean(model$VCV[i,'at.level(trait, "r").id:time'])), ncol = 2)
	Ppred <- Gpred + BYpred + Matpred + RIpred + Rpred
#	Ppred <- Gpred + BYpred + RIpred + Rpred
	
	# Transform from the latent scale to the data scale
	obs <- QGmvparams(predict=pred,vcv.G=Gpred,vcv.P=Ppred,model=c("binom1.probit","Gaussian"))
	
	# Results to save
	CovaSZ <- obs$vcv.G.obs[1,2] # Genetic covariance
	Saz1 <- obs$vcv.G.obs[1,2]/obs$mean.obs[1] ## Prior viability genetic selection differential trait
	
	res <- c(CovaSZ,Saz1,i)
	
	DataScaleRes <- rbind(DataScaleRes,res)	
}
DataScaleRes <- as.data.frame(DataScaleRes)
names(DataScaleRes) <- c("CovaSZ","Saz1","i")

save(DataScaleRes, file=paste(file_path_output,"DataScaleResMSCcovuFYS_malesonly_120525.RData",sep="")) # Date as appropriate


############################################################################################################################################################
## Adult size traits -- heritability, selection (using only information on adults) and predicted response to selection (using only information on adults) ##
############################################################################################################################################################
############################
## Adult August body mass ##
############################
# Heritability #
# Load model data
data <- read.csv(paste(file_path_data,"UniVarABM_h2_RandomIds_May2025.csv",sep=""))

# Format
data$animal <- as.factor(data$animal)
data$ID <- as.factor(data$ID)
data$CapYear <- as.factor(data$CapYear)
data$BirthYear <- as.factor(data$BirthYear)
data$CapAgeYears <- as.factor(data$CapAgeYears)
str(data)
summary(data)

# Prune pedigree
pruned <- prunePed(pedigree,keep=data$ID)

# Univariate animal model
priorAAW1 <- list(R = list(R1=list(V=1, nu=2.002)), G=list(G1=list(V=1,nu=2.002),G2=list(V=1,nu=2.002),G3=list(V=1,nu=2.002),G4=list(V=1,nu=2.002)))
UniAAW1 <- MCMCglmm(Weight ~ CapAgeYears + Sex + Twin + jdateMC + PopDensMC, random=~animal + CapYear + ID + BirthYear, rcov=~units, family="gaussian", data=data, pedigree=pruned, prior=priorAAW1,nitt=10,thin=1,burnin=1)#, nitt=1505000,thin=500,burnin=5000)
save(UniAAW1,file="UniAAW1_120525.RData")

# Selection #
# Load model data
data <- read.csv(paste(file_path_data,"DataAdultPhenoSel_ABM_RandomIds_May2025.csv",sep=""))

# Format
data$id <- as.factor(data$id)
data$trait <- as.factor(data$trait)
data$CapAgeYears <- as.factor(data$CapAgeYears)
str(data)
summary(data)

# Bivariate phenotypic models
# Both sexes
priorC <- list(R = list(R1=list(V=diag(2), nu=2.002, covu=TRUE), R2=list(V=1, nu=2.002)))
BivModelCovuAAM <- MCMCglmm(y ~ trait - 1 + at.level(trait,"r"):CapAgeYears + trait:Sex + trait:Twin + at.level(trait,"r"):jdateMC + trait:PopDensMC, random=~idh(at.level(trait,"r")):id, rcov=~idh(at.level(trait,"s")):id + idh(at.level(trait,"r")):id:time, family = NULL, data=data, prior=priorC, nitt=155000,thin=50,burnin=5000)
save(BivModelCovuAAM,file="phenSelCovuAAM_130525.RData")

# Males
malesC <- data[which(data$Sex==2),]
priorC <- list(R = list(R1=list(V=diag(2), nu=2.002, covu=TRUE), R2=list(V=1, nu=2.002)))
BivModelCovuAAMM <- MCMCglmm(y ~ trait - 1 + at.level(trait,"r"):CapAgeYears + trait:Twin + at.level(trait,"r"):jdateMC + trait:PopDensMC, random=~idh(at.level(trait,"r")):id, rcov=~idh(at.level(trait,"s")):id + idh(at.level(trait,"r")):id:time, family = NULL, data=malesC, prior=priorC nitt=155000,thin=50,burnin=5000)
save(BivModelCovuAAMM,file="phenSelCovuAAM_males_130525.RData")

# Males Age cat 
femalesC <- data[which(data$Sex==1),]
priorC <- list(R = list(R1=list(V=diag(2), nu=2.002, covu=TRUE), R2=list(V=1, nu=2.002)))
BivModelCovuAAMF <- MCMCglmm(y ~ trait - 1 + at.level(trait,"r"):CapAgeYears + trait:Twin + at.level(trait,"r"):jdateMC + trait:PopDensMC, random=~idh(at.level(trait,"r")):id, rcov=~idh(at.level(trait,"s")):id + idh(at.level(trait,"r")):id:time, family = NULL, data=femalesC, prior=priorC, nitt=155000,thin=50,burnin=5000)
save(BivModelCovuAAMF,file="phenSelCovuAAM_females_130525.RData")

#############################
## Adult metatarsal length ##
#############################
# Load model data
data <- read.csv(paste(file_path_data,"DataAdultPhenoSel_AMTL_RandomIds_May2025.csv",sep=""))

# Format
data$id <- as.factor(data$id)
data$trait <- as.factor(data$trait)
data$CapAgeYears <- as.factor(data$CapAgeYears)
str(data)
summary(data)

# Bivariate phenotypic models
# Both sexes
priorC <- list(R = list(R1=list(V=diag(2), nu=2.002, covu=TRUE), R2=list(V=1, nu=2.002)))
BivModelCovuMTL <- MCMCglmm(y ~ trait - 1 + at.level(trait,"r"):CapAgeYears + trait:Sex + at.level(trait,"r"):Twin + at.level(trait,"r"):jdateMC + trait:PopDensMC, random=~idh(at.level(trait,"r")):id, rcov=~idh(at.level(trait,"s")):id + idh(at.level(trait,"r")):id:time, family = NULL, data=data, prior=priorC,nitt=26000,thin=10,burnin=3000)
save(BivModelCovuMTL,file="phenSelCovuMTL_130525.RData")

# Males
malesC <- data[which(data$Sex==2),]
priorC <- list(R = list(R1=list(V=diag(2), nu=2.002, covu=TRUE), R2=list(V=1, nu=2.002)))
BivModelCovuMTLM <- MCMCglmm(y ~ trait - 1 + at.level(trait,"r"):CapAgeYears + at.level(trait,"r"):Twin + at.level(trait,"r"):jdateMC + trait:PopDensMC, random=~idh(at.level(trait,"r")):id, rcov=~idh(at.level(trait,"s")):id + idh(at.level(trait,"r")):id:time, family = NULL, data=malesC, prior=priorC,nitt=46000,thin=10,burnin=3000)
save(BivModelCovuMTLM,file="phenSelCovuMTL_males_130525.RData")

# Females
femalesC <- data[which(data$Sex==1),]
priorC <- list(R = list(R1=list(V=diag(2), nu=2.002, covu=TRUE), R2=list(V=1, nu=2.002)))
BivModelCovuMTLF <- MCMCglmm(y ~ trait - 1 + at.level(trait,"r"):CapAgeYears + at.level(trait,"r"):Twin + at.level(trait,"r"):jdateMC + trait:PopDensMC, random=~idh(at.level(trait,"r")):id, rcov=~idh(at.level(trait,"s")):id + idh(at.level(trait,"r")):id:time, family = NULL, data=femalesC, prior=priorC,nitt=26000,thin=10,burnin=3000)
save(BivModelCovuMTLF,file="phenSelCovuMTL_females_130525.RData")


######################
## Male horn length ##
######################
# Load model data
data <- read.csv(paste(file_path_data,"DataAdultPhenoSel_MHL_RandomIds_May2025.csv",sep=""))

# Format
data$id <- as.factor(data$id)
data$trait <- as.factor(data$trait)
data$CapAgeYears <- as.factor(data$CapAgeYears)
str(data)
summary(data)

# Bivariate phenotypic models
priorC <- list(R = list(R1=list(V=diag(2), nu=2.002, covu=TRUE), R2=list(V=1, nu=2.002)))
BivModelCovuMHL <- MCMCglmm(y ~ trait - 1 + at.level(trait,"r"):CapAgeYears + trait:Twin + at.level(trait,"r"):jdateMC + trait:PopDensMC, random=~idh(at.level(trait,"r")):id, rcov=~idh(at.level(trait,"s")):id + idh(at.level(trait,"r")):id:time, family = NULL, data=data, prior=priorC,nitt=26000,thin=10,burnin=3000)
save(BivModelCovuMHL,file="phenSelCovuMHL_130525.RData") 


###########################
## Scrotal circumference ##
###########################
# Load model data
data <- read.csv(paste(file_path_data,"DataAdultPhenoSel_MSC_RandomIds_May2025.csv",sep=""))

# Format
data$id <- as.factor(data$id)
data$trait <- as.factor(data$trait)
data$CapAgeYears <- as.factor(data$CapAgeYears)
str(data)
summary(data)

# Bivariate phenotypic models
priorC <- list(R = list(R1=list(V=diag(2), nu=2.002, covu=TRUE), R2=list(V=1, nu=2.002)))
BivModelCovuMSC <- MCMCglmm(y ~ trait - 1 + at.level(trait,"r"):CapAgeYears + trait:Twin + at.level(trait,"r"):jdateMC + trait:PopDensMC, random=~idh(at.level(trait,"r")):id, rcov=~idh(at.level(trait,"s")):id + idh(at.level(trait,"r")):id:time, family = NULL, data=data, prior=priorC,nitt=36000,thin=20,burnin=3000)
save(BivModelCovuMSC,file="phenSelCovuMSC_130525.RData") 

