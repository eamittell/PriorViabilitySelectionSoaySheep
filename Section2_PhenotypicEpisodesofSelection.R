#### Script for models from section 2 of the manuscript: Quantitative genetic analysis of early life viability resolves a case of the paradox of stasis in body size traits of Soay sheep
### Lizy Mittell
## Section 2: Phenotypic episodes of selection
# Cleaned 9th May 2025

# Libraries
library(MCMCglmm)
library(mvtnorm)

# File path to data -- edit as appropriate
#file_path_data <- ""
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

#############################################
## Viability selection differentials lambs ##
#############################################
####################
# August body mass # 
####################
# Load data
data <- read.csv(paste(file_path_data,"DataPhenoSel_LABM_RandomIds_May2025.csv",sep=""))

# Bivariate models of lamb august weight with first year survival
## 2815 lamb august weights (from 2815 lambs as only one obs)
prior1 <- list(R = list(R1=list(V=diag(2), nu=2.002,fix=2)))
modelBivViabLAM <- MCMCglmm(cbind(LambAugWeight,Sur1) ~ trait - 1 + trait:Twin + trait:Sex + trait:PopDensMC + at.level(trait,"LambAugWeight"):jdateMC + at.level(trait,"LambAugWeight"):birthjdateMC, rcov=~us(trait):units, family = c("gaussian","threshold"), data=data, prior=prior1,nitt= 103000,burnin=3000,thin=10)
save(modelBivViabLAM,file=paste(file_path_output,"BivViab_LAM_191224.RData",sep=""))

# Males
maleslaw <- data[which(data$Sex==2),]
## 1369
prior1 <- list(R = list(R1=list(V=diag(2), nu=2.002,fix=2)))
modelBivViabLAMM <- MCMCglmm(cbind(LambAugWeight,Sur1) ~ trait - 1 + trait:Twin + trait:PopDensMC + at.level(trait,"LambAugWeight"):jdateMC + at.level(trait,"LambAugWeight"):birthjdateMC, rcov=~us(trait):units, family = c("gaussian","threshold"), data=maleslaw, prior=prior1,nitt= 103000,burnin=3000,thin=10)
save(modelBivViabLAMM,file=paste(file_path_output,"BivViabM_LAM_191224.RData",sep=""))

# Females
femaleslaw <- data[which(data$Sex==1),]
## 1446
prior1 <- list(R = list(R1=list(V=diag(2), nu=2.002,fix=2)))
modelBivViabLAMF <- MCMCglmm(cbind(LambAugWeight,Sur1) ~ trait - 1 + trait:Twin + trait:PopDensMC + at.level(trait,"LambAugWeight"):jdateMC + at.level(trait,"LambAugWeight"):birthjdateMC, rcov=~us(trait):units, family = c("gaussian","threshold"), data=femaleslaw, prior=prior1,nitt= 103000,burnin=3000,thin=10)
save(modelBivViabLAMF,file=paste(file_path_output,"BivViabF_LAM_191224.RData",sep=""))

### Back-transformation
# Both sexes
load(paste(file_path_output,"BivViab_LAM_191224.RData",sep=""))
model <- modelBivViabLAM
nmc <- 10000

ViabResDataScaleLAM <- NULL
for(i in 1:nmc){
# Extraction of the fitted values
# Survival
mu.s.Fixed <- mean(model$Sol[i,2])*model$X[,2] + mean(model$Sol[i,4])*model$X[,4] + mean(model$Sol[i,6])*model$X[,6] + mean(model$Sol[i,8])*model$X[,8]
# Repeat measure trait
mu.r.Fixed <- mean(model$Sol[i,1])*model$X[,1] + mean(model$Sol[i,3])*model$X[,3] + mean(model$Sol[i,5])*model$X[,5] + mean(model$Sol[i,7])*model$X[,7] + mean(model$Sol[i,9])*model$X[,9] + mean(model$Sol[i,10])*model$X[,10]
		
# Phenotypic co-variance matrix
Ppred <- matrix(c(mean(model$VCV[i,"traitSur1:traitSur1.units"]), mean(model$VCV[i,"traitSur1:traitLambAugWeight.units"]), mean(model$VCV[i,"traitSur1:traitLambAugWeight.units"]), mean(model$VCV[i,"traitLambAugWeight:traitLambAugWeight.units"])), ncol = 2)

# Approximation
nlambs <- 5630
p.nm <- rmvnorm(nlambs,c(0,0),Ppred)
muS <- mean(mu.s.Fixed) # Mean survival on the latent scale
sur.lat <- p.nm[,1] + rnorm(nlambs) + mu.s.Fixed
# Convert to data scale
data.s <- (sur.lat>0)+0
# Calculate relative survival
prior.s <- data.s/mean(data.s)
# Calculate covariance with trait
Sz_approx <- cov(p.nm[,2],prior.s)

res <- c(muS,Sz_approx,i)

ViabResDataScaleLAM <- rbind(ViabResDataScaleLAM,res)
}
ViabResDataScaleLAM <- as.data.frame(ViabResDataScaleLAM)
names(ViabResDataScaleLAM) <- c("muS","Sz_approx","i")
save(ViabResDataScaleLAM,file=paste(file_path_output,"ViabResDataScaleLAM.RData",sep=""))

# Males
load(paste(file_path_output,"BivViabM_LAM_191224.RData",sep=""))
model <- modelBivViabLAMM
nmc <- 10000

ViabResDataScaleLAMM <- NULL
for(i in 1:nmc){
# Extraction of the fitted values
# Survival
mu.s.Fixed <- mean(model$Sol[i,2])*model$X[,2] + mean(model$Sol[i,4])*model$X[,4] + mean(model$Sol[i,6])*model$X[,6]
# Repeat measure trait
mu.r.Fixed <- mean(model$Sol[i,1])*model$X[,1] + mean(model$Sol[i,3])*model$X[,3] + mean(model$Sol[i,5])*model$X[,5] + mean(model$Sol[i,7])*model$X[,7] + mean(model$Sol[i,8])*model$X[,8]
		
# Phenotypic co-variance matrix
Ppred <- matrix(c(mean(model$VCV[i,"traitSur1:traitSur1.units"]), mean(model$VCV[i,"traitSur1:traitLambAugWeight.units"]), mean(model$VCV[i,"traitSur1:traitLambAugWeight.units"]), mean(model$VCV[i,"traitLambAugWeight:traitLambAugWeight.units"])), ncol = 2)

# Approximation
nlambs <- 2738
p.nm <- rmvnorm(nlambs,c(0,0),Ppred)
muS <- mean(mu.s.Fixed) # Mean survival on the latent scale
sur.lat <- p.nm[,1] + rnorm(nlambs) + mu.s.Fixed
# Convert to data scale
data.s <- (sur.lat>0)+0
# Calculate relative survival
prior.s <- data.s/mean(data.s)
# Calculate covariance with trait
Sz_approx <- cov(p.nm[,2],prior.s)

res <- c(muS,Sz_approx,i)

ViabResDataScaleLAMM <- rbind(ViabResDataScaleLAMM,res)
}
ViabResDataScaleLAMM <- as.data.frame(ViabResDataScaleLAMM)
names(ViabResDataScaleLAMM) <- c("muS","Sz_approx","i")
save(ViabResDataScaleLAMM,file=paste(file_path_output,"ViabResDataScaleLAM_Males.RData",sep=""))

# Females
load(paste(file_path_output,"BivViabF_LAM_191224.RData",sep=""))
model <- modelBivViabLAMF
nmc <- 10000

ViabResDataScaleLAMF <- NULL
for(i in 1:nmc){
# Extraction of the fitted values
# Survival
mu.s.Fixed <- mean(model$Sol[i,2])*model$X[,2] + mean(model$Sol[i,4])*model$X[,4] + mean(model$Sol[i,6])*model$X[,6]
# Repeat measure trait
mu.r.Fixed <- mean(model$Sol[i,1])*model$X[,1] + mean(model$Sol[i,3])*model$X[,3] + mean(model$Sol[i,5])*model$X[,5] + mean(model$Sol[i,7])*model$X[,7] + mean(model$Sol[i,8])*model$X[,8]
		
# Phenotypic co-variance matrix
Ppred <- matrix(c(mean(model$VCV[i,"traitSur1:traitSur1.units"]), mean(model$VCV[i,"traitSur1:traitLambAugWeight.units"]), mean(model$VCV[i,"traitSur1:traitLambAugWeight.units"]), mean(model$VCV[i,"traitLambAugWeight:traitLambAugWeight.units"])), ncol = 2)

# Approximation
nlambs <- 2892
p.nm <- rmvnorm(nlambs,c(0,0),Ppred) 
muS <- mean(mu.s.Fixed) # Mean survival on the latent scale
sur.lat <- p.nm[,1] + rnorm(nlambs) + mu.s.Fixed 
# Convert to data scale
data.s <- (sur.lat>0)+0 
# Calculate relative survival
prior.s <- data.s/mean(data.s)
# Calculate covariance with trait
Sz_approx <- cov(p.nm[,2],prior.s) 

res <- c(muS,Sz_approx,i)

ViabResDataScaleLAMF <- rbind(ViabResDataScaleLAMF,res)
}
ViabResDataScaleLAMF <- as.data.frame(ViabResDataScaleLAMF)
names(ViabResDataScaleLAMF) <- c("muS","Sz_approx","i")
save(ViabResDataScaleLAMF,file=paste(file_path_output,"ViabResDataScaleLAM_Females.RData",sep=""))

############################
# August metatarsal length #
############################
# Load data
data <- read.csv(paste(file_path_data,"DataPhenoSel_LAMTL_RandomIds_May2025.csv",sep=""))

# Bivariate models of lamb august metatarsal length with first year survival
## 2711 lamb august weights (from 2711 lambs as only one obs)
prior1 <- list(R = list(R1=list(V=diag(2), nu=2.002,fix=2)))
modelBivViabLAMTL <- MCMCglmm(cbind(LambHindLeg,Sur1) ~ trait - 1 + trait:Twin + trait:Sex + trait:PopDensMC + at.level(trait,"LambHindLeg"):jdateMC + at.level(trait,"LambHindLeg"):birthjdateMC, rcov=~us(trait):units, family = c("gaussian","threshold"), data=data, prior=prior1,nitt= 103000,burnin=3000,thin=10)
save(modelBivViabLAMTL,file=paste(file_path_output,"BivViab_LAMTL_191224.RData",sep=""))

# Males
maleslhll <- data[which(data$Sex==2),]
## 1313
prior1 <- list(R = list(R1=list(V=diag(2), nu=2.002,fix=2)))
modelBivViabLAMTLM <- MCMCglmm(cbind(LambHindLeg,Sur1) ~ trait - 1 + trait:Twin + trait:PopDensMC + at.level(trait,"LambHindLeg"):jdateMC + at.level(trait,"LambHindLeg"):birthjdateMC, rcov=~us(trait):units, family = c("gaussian","threshold"), data= maleslhll, prior=prior1,nitt=10,thin=1,burnin=1)#,nitt= 103000,burnin=3000,thin=10)
save(modelBivViabLAMTLM,file=paste(file_path_output,"BivViabM_LAMTL_191224.RData",sep=""))

# Females
femaleslhll <- data[which(data$Sex==1),]
## 1398
prior1 <- list(R = list(R1=list(V=diag(2), nu=2.002,fix=2)))
modelBivViabLAMTLF <- MCMCglmm(cbind(LambHindLeg,Sur1) ~ trait - 1 + trait:Twin + trait:PopDensMC + at.level(trait,"LambHindLeg"):jdateMC + at.level(trait,"LambHindLeg"):birthjdateMC, rcov=~us(trait):units, family = c("gaussian","threshold"), data= femaleslhll, prior=prior1,nitt=10,thin=1,burnin=1)#,nitt= 103000,burnin=3000,thin=10)
save(modelBivViabLAMTLF,file=paste(file_path_output,"BivViabF_LAMTL_191224.RData",sep=""))

### Back-transformation
# Both sexes
load(paste(file_path_output,"BivViab_LAMTL_191224.RData",sep=""))
model <- modelBivViabLAMTL
nmc <- 10000

ViabResDataScaleLAM <- NULL
for(i in 1:nmc){
# Extraction of the fitted values
# Survival
mu.s.Fixed <- mean(model$Sol[i,2])*model$X[,2] + mean(model$Sol[i,4])*model$X[,4] + mean(model$Sol[i,6])*model$X[,6] + mean(model$Sol[i,8])*model$X[,8]
# Repeat measure trait
mu.r.Fixed <- mean(model$Sol[i,1])*model$X[,1] + mean(model$Sol[i,3])*model$X[,3] + mean(model$Sol[i,5])*model$X[,5] + mean(model$Sol[i,7])*model$X[,7] + mean(model$Sol[i,9])*model$X[,9] + mean(model$Sol[i,10])*model$X[,10]
		
# Phenotypic co-variance matrix
Ppred <- matrix(c(mean(model$VCV[i,"traitSur1:traitSur1.units"]), mean(model$VCV[i,"traitSur1:traitLambHindLeg.units"]), mean(model$VCV[i,"traitSur1:traitLambHindLeg.units"]), mean(model$VCV[i,"traitLambHindLeg:traitLambHindLeg.units"])), ncol = 2)

# Approximation
nlambs <- 5422
p.nm <- rmvnorm(nlambs,c(0,0),Ppred)
muS <- mean(mu.s.Fixed) # Mean survival on the latent scale
sur.lat <- p.nm[,1] + rnorm(nlambs) + mu.s.Fixed
# Convert to data scale
data.s <- (sur.lat>0)+0
# Calculate relative survival
prior.s <- data.s/mean(data.s)
# Calculate covariance with trait
Sz_approx <- cov(p.nm[,2],prior.s)

res <- c(muS,Sz_approx,i)

ViabResDataScaleLAMTL <- rbind(ViabResDataScaleLAMTL,res)
}
ViabResDataScaleLAMTL <- as.data.frame(ViabResDataScaleLAMTL)
names(ViabResDataScaleLAMTL) <- c("muS","Sz_approx","i")
save(ViabResDataScaleLAMTL,file=paste(file_path_output,"ViabResDataScaleLAMTL.RData",sep=""))

# Males
load(paste(file_path_output,"BivViabM_LAMTL_191224.RData",sep=""))
model <- modelBivViabLLAMTLM
nmc <- 10000

ViabResDataScaleLAMTLM <- NULL
for(i in 1:nmc){
# Extraction of the fitted values
# Survival
mu.s.Fixed <- mean(model$Sol[i,2])*model$X[,2] + mean(model$Sol[i,4])*model$X[,4] + mean(model$Sol[i,6])*model$X[,6]
# Repeat measure trait
mu.r.Fixed <- mean(model$Sol[i,1])*model$X[,1] + mean(model$Sol[i,3])*model$X[,3] + mean(model$Sol[i,5])*model$X[,5] + mean(model$Sol[i,7])*model$X[,7] + mean(model$Sol[i,8])*model$X[,8]
		
# Phenotypic co-variance matrix
Ppred <- matrix(c(mean(model$VCV[i,"traitSur1:traitSur1.units"]), mean(model$VCV[i,"traitSur1:traitLambHindLeg.units"]), mean(model$VCV[i,"traitSur1:traitLambHindLeg.units"]), mean(model$VCV[i,"traitLambHindLeg:traitLambHindLeg.units"])), ncol = 2)

# Approximation
nlambs <- 2626
p.nm <- rmvnorm(nlambs,c(0,0),Ppred)
muS <- mean(mu.s.Fixed) # Mean survival on the latent scale
sur.lat <- p.nm[,1] + rnorm(nlambs) + mu.s.Fixed
# Convert to data scale
data.s <- (sur.lat>0)+0
# Calculate relative survival
prior.s <- data.s/mean(data.s)
# Calculate covariance with trait
Sz_approx <- cov(p.nm[,2],prior.s)

res <- c(muS,Sz_approx,i)

ViabResDataScaleLAMTLM <- rbind(ViabResDataScaleLAMTLM,res)
}
ViabResDataScaleLAMTLM <- as.data.frame(ViabResDataScaleLAMTLM)
names(ViabResDataScaleLAMTLM) <- c("muS","Sz_approx","i")
save(ViabResDataScaleLAMTLM,file=paste(file_path_output,"ViabResDataScaleLAMTLM_Males.RData",sep=""))

# Females
load(paste(file_path_output,"BivViabF_LAMTL_191224.RData",sep=""))
model <- modelBivViabLAMTLF
nmc <- 10000

ViabResDataScaleLAMTLF <- NULL
for(i in 1:nmc){
# Extraction of the fitted values
# Survival
mu.s.Fixed <- mean(model$Sol[i,2])*model$X[,2] + mean(model$Sol[i,4])*model$X[,4] + mean(model$Sol[i,6])*model$X[,6]
# Repeat measure trait
mu.r.Fixed <- mean(model$Sol[i,1])*model$X[,1] + mean(model$Sol[i,3])*model$X[,3] + mean(model$Sol[i,5])*model$X[,5] + mean(model$Sol[i,7])*model$X[,7] + mean(model$Sol[i,8])*model$X[,8]
		
# Phenotypic co-variance matrix
Ppred <- matrix(c(mean(model$VCV[i,"traitSur1:traitSur1.units"]), mean(model$VCV[i,"traitSur1:traitLambHindLeg.units"]), mean(model$VCV[i,"traitSur1:traitLambHindLeg.units"]), mean(model$VCV[i,"traitLambHindLeg:traitLambHindLeg.units"])), ncol = 2)

# Approximation
nlambs <- 2796
p.nm <- rmvnorm(nlambs,c(0,0),Ppred) 
muS <- mean(mu.s.Fixed) # Mean survival on the latent scale
sur.lat <- p.nm[,1] + rnorm(nlambs) + mu.s.Fixed 
# Convert to data scale
data.s <- (sur.lat>0)+0 
# Calculate relative survival
prior.s <- data.s/mean(data.s)
# Calculate covariance with trait
Sz_approx <- cov(p.nm[,2],prior.s) 

res <- c(muS,Sz_approx,i)

ViabResDataScaleLAMTLF <- rbind(ViabResDataScaleLAMTLF,res)
}
ViabResDataScaleLAMTLF <- as.data.frame(ViabResDataScaleLAMTLF)
names(ViabResDataScaleLAMTLF) <- c("muS","Sz_approx","i")
save(ViabResDataScaleLAMTLF,file=paste(file_path_output,"ViabResDataScaleLAMTL_Females.RData",sep=""))

####################
# Male horn length #
####################
# Load data
data <- read.csv(paste(file_path_data,"DataPhenoSel_LMHL_RandomIds_May2025.csv",sep=""))

# Bivariate models of lamb male horn length with first year survival
## 1178 lamb  male horn length (from 1178 lambs as only one obs)
prior1 <- list(R = list(R1=list(V=diag(2), nu=2.002,fix=2)))
modelBivViabLMHL <- MCMCglmm(cbind(LambHornLen,Sur1) ~ trait - 1 + trait:Twin + trait:PopDensMC + at.level(trait,"LambHornLen"):jdateMC + at.level(trait,"LambHornLen"):birthjdateMC, rcov=~us(trait):units, family = c("gaussian","threshold"), data= malelambHL2, prior=prior1,nitt=103000,burnin=3000,thin=10)
save(modelBivViabLMHL,file="BivViab_LMHL_191224.RData") 

### Back-transformation
model <- modelBivViabLMHL
nmc <- 10000

ViabResDataScaleLMHL <- NULL
for(i in 1:nmc){
# Extraction of the fitted values
# Survival
mu.s.Fixed <- mean(model$Sol[i,2])*model$X[,2] + mean(model$Sol[i,4])*model$X[,4] + mean(model$Sol[i,6])*model$X[,6]
# Repeat measure trait
mu.r.Fixed <- mean(model$Sol[i,1])*model$X[,1] + mean(model$Sol[i,3])*model$X[,3] + mean(model$Sol[i,5])*model$X[,5] + mean(model$Sol[i,7])*model$X[,7] + mean(model$Sol[i,8])*model$X[,8]
		
# Phenotypic co-variance matrix
Ppred <- matrix(c(mean(model$VCV[i,"traitSur1:traitSur1.units"]), mean(model$VCV[i,"traitSur1:traitLambHornLen.units"]), mean(model$VCV[i,"traitSur1:traitLambHornLen.units"]), mean(model$VCV[i,"traitLambHornLen:traitLambHornLen.units"])), ncol = 2)

# Approximation
nlambs <- 2356
p.nm <- rmvnorm(nlambs,c(0,0),Ppred) 
muS <- mean(mu.s.Fixed) # Mean survival on the latent scale 
sur.lat <- p.nm[,1] + rnorm(nlambs) + mu.s.Fixed
# Convert to data scale
data.s <- (sur.lat>0)+0 
# Calculate relative survival
prior.s <- data.s/mean(data.s)
# Calculate covariance with trait
Sz_approx <- cov(p.nm[,2],prior.s) 

res <- c(muS,Sz_approx,i)

ViabResDataScaleLMHL <- rbind(ViabResDataScaleLMHL,res)
}
ViabResDataScaleLMHL <- as.data.frame(ViabResDataScaleLMHL)
names(ViabResDataScaleLMHL) <- c("muS","Sz_approx","i")
save(ViabResDataScaleLMHL,file=paste(file_path_output,"ViabResDataScaleLMHL.RData",sep=""))


##############################
# Male scrotal circumference #
##############################
# Load data
data <- read.csv(paste(file_path_data,"DataPhenoSel_LMSC_RandomIds_May2025.csv",sep=""))

# Bivariate models of lamb male scrotal circ with first year survival
## 1145 lamb  male scrotal circ (from 1145 lambs as only one obs)
prior1 <- list(R = list(R1=list(V=diag(2), nu=2.002,fix=2)))
modelBivViabLMSC <- MCMCglmm(cbind(LambBolCirc,Sur1) ~ trait - 1 + trait:Twin + trait:PopDensMC + at.level(trait,"LambBolCirc"):jdateMC + at.level(trait,"LambBolCirc"):birthjdateMC, rcov=~us(trait):units, family = c("gaussian","threshold"), data= lambMSC2, prior=prior1,nitt=103000,burnin=3000,thin=10)
save(modelBivViabLMSC,file="BivViab_LMSC_191224.RData") #3214

### Back-transformation
model <- modelBivViabLMSC
nmc <- 10000

ViabResDataScaleLMSC <- NULL
for(i in 1:nmc){
# Extraction of the fitted values
# Survival
mu.s.Fixed <- mean(model$Sol[i,2])*model$X[,2] + mean(model$Sol[i,4])*model$X[,4] + mean(model$Sol[i,6])*model$X[,6]
# Repeat measure trait
mu.r.Fixed <- mean(model$Sol[i,1])*model$X[,1] + mean(model$Sol[i,3])*model$X[,3] + mean(model$Sol[i,5])*model$X[,5] + mean(model$Sol[i,7])*model$X[,7] + mean(model$Sol[i,8])*model$X[,8]
		
# Phenotypic co-variance matrix
Ppred <- matrix(c(mean(model$VCV[i,"traitSur1:traitSur1.units"]), mean(model$VCV[i,"traitSur1:traitLambBolCirc.units"]), mean(model$VCV[i,"traitSur1:traitLambBolCirc.units"]), mean(model$VCV[i,"traitLambBolCirc:traitLambBolCirc.units"])), ncol = 2)

# Approximation
nlambs <- 2290
p.nm <- rmvnorm(nlambs,c(0,0),Ppred) 
muS <- mean(mu.s.Fixed) # Mean survival on the latent scale
sur.lat <- p.nm[,1] + rnorm(nlambs) + mu.s.Fixed 
# Convert to data scale
data.s <- (sur.lat>0)+0 
# Calculate relative survival
prior.s <- data.s/mean(data.s)
# Calculate covariance with trait
Sz_approx <- cov(p.nm[,2],prior.s) 

res <- c(muS,Sz_approx,i)

ViabResDataScaleLMSC <- rbind(ViabResDataScaleLMSC,res)
}
ViabResDataScaleLMSC <- as.data.frame(ViabResDataScaleLMSC)
names(ViabResDataScaleLMSC) <- c("muS","Sz_approx","i")
save(ViabResDataScaleLMSC,file=paste(file_path_output,"ViabResDataScaleLMSC.RData",sep=""))


##############################################################################
### Phenotypic covariance between lamb august weight and adult body weight ###
##############################################################################
####################
# August body mass # 
####################
# Load data
data <- read.csv(paste(file_path_data,"DataPhenoCovu_ABM_LABM_RandomIds_May2025.csv",sep=""))

# Format
data$CapAgeYears <- as.factor(data$CapAgeYears)
data$trait <- as.factor(data$trait)
data$id <- as.factor(data$id)
summary(data)

# Bivariate models of lamb august weight with later-life body weight
# Both sexes
prior2 <- list(R = list(R1=list(V=diag(2), nu=2.002, covu=TRUE), R2=list(V=1, nu=2.002)))
modelCovuLAW <- MCMCglmm(y ~ trait - 1 + trait:Twin + trait:PopDensMC + trait:jdateMC + at.level(trait,"s"):birthjdateMC  + trait:Sex + at.level(trait,"r"):CapAgeYears + at.level(trait,"s"):MumAge + at.level(trait,"r"):time, random=~idh(at.level(trait,"r")):id, rcov=~idh(at.level(trait,"s")):id + idh(at.level(trait,"r")):id:time, family = NULL, data=data, prior=prior2,nitt=66000,thin=50,burnin=3000)
save(modelCovuLAW,file="phenCovP_LAW_ABW_130125.RData")

## Males
malesC <- data[which(data$Sex==2),]
prior2 <- list(R = list(R1=list(V=diag(2), nu=2.002, covu=TRUE), R2=list(V=1, nu=2.002)))
modelCovuLAWmales <- MCMCglmm(y ~ trait - 1 + trait:Twin + trait:PopDensMC + trait:jdateMC + at.level(trait,"s"):birthjdateMC + at.level(trait,"r"):CapAgeYears + at.level(trait,"s"):MumAge + at.level(trait,"r"):time, random=~idh(at.level(trait,"r")):id, rcov=~idh(at.level(trait,"s")):id + idh(at.level(trait,"r")):id:time, family = NULL, data=malesC, prior=prior2, nitt=206000,thin=100,burnin=3000)
save(modelCovuLAWmales,file="phenCovP_LAW_ABW_males_130125.RData") # 

## Females
femalesC <- data[which(data$Sex==1),]
prior2 <- list(R = list(R1=list(V=diag(2), nu=2.002, covu=TRUE), R2=list(V=1, nu=2.002)))
modelCovuLAWfemales <- MCMCglmm(y ~ trait - 1 + trait:Twin + trait:PopDensMC + trait:jdateMC + at.level(trait,"s"):birthjdateMC + at.level(trait,"r"):CapAgeYears + at.level(trait,"s"):MumAge + at.level(trait,"r"):time, random=~idh(at.level(trait,"r")):id, rcov=~idh(at.level(trait,"s")):id + idh(at.level(trait,"r")):id:time, family = NULL, data=femalesC, prior=prior2, nitt=206000,thin=100,burnin=3000)
save(modelCovuLAWfemales,file="phenCovP_LAW_ABW_females_130125.RData") # 

############################
# August metatarsal length #
############################
# Load data
data <- read.csv(paste(file_path_data,"DataPhenoCovu_AMTL_LMTL_RandomIds_May2025.csv",sep=""))

# Format
data$CapAgeYears <- as.factor(data$CapAgeYears)
data$trait <- as.factor(data$trait)
data$id <- as.factor(data$id)
summary(data)

# Bivariate models of lamb august metatarsal length with later-life metatarsal length
prior2 <- list(R = list(R1=list(V=diag(2), nu=2.002, covu=TRUE), R2=list(V=1, nu=2.002)))
modelCovuLHLL <- MCMCglmm(y ~ trait - 1 + trait:Twin + trait:PopDensMC + trait:jdateMC + at.level(trait,"s"):birthjdateMC  + trait:Sex + at.level(trait,"r"):CapAgeYears + at.level(trait,"s"):MumAge + at.level(trait,"r"):time, random=~idh(at.level(trait,"r")):id, rcov=~idh(at.level(trait,"s")):id + idh(at.level(trait,"r")):id:time, family = NULL, data=data, prior=prior2,nitt=66000,thin=50,burnin=3000)
save(modelCovuLHLL,file="phenCovP_LHLL_AHLL_130125.RData")

## Males
malesC <- data[which(data$Sex==2),]
length(unique(malesC$id[which(malesC$trait=="s")])) ## 1273 lamb august weights (from 1273 lambs as only one obs)
dim(malesC[which(malesC$trait=="r"),])
length(unique(malesC$id[which(malesC$trait=="r")])) ## 1134 obs adult august weight from 597 adult sheep (age 2+ -- 15)
prior2 <- list(R = list(R1=list(V=diag(2), nu=2.002, covu=TRUE), R2=list(V=1, nu=2.002)))
modelCovuLHLLmales <- MCMCglmm(y ~ trait - 1 + trait:Twin + trait:PopDensMC + trait:jdateMC + at.level(trait,"s"):birthjdateMC + at.level(trait,"r"):CapAgeYears + at.level(trait,"s"):MumAge + at.level(trait,"r"):time, random=~idh(at.level(trait,"r")):id, rcov=~idh(at.level(trait,"s")):id + idh(at.level(trait,"r")):id:time, family = NULL, data=malesC, prior=prior2,nitt=1003000,thin=500,burnin=3000)
save(modelCovuLHLLmales,file="phenCovP_LHLL_AHLL_males_130125.RData")

## Females
femalesC <- data[which(data$Sex==1),]
length(unique(femalesC$id[which(femalesC$trait=="s")])) ## 1358 lamb august weights (from 1358 lambs as only one obs)
dim(femalesC[which(femalesC$trait=="r"),])
length(unique(femalesC$id[which(femalesC$trait=="r")])) ## 3639 obs adult august weight from 1097 adult sheep (age 2+ -- 15)
prior2 <- list(R = list(R1=list(V=diag(2), nu=2.002, covu=TRUE), R2=list(V=1, nu=2.002)))
modelCovuLHLLfemales <- MCMCglmm(y ~ trait - 1 + trait:Twin + trait:PopDensMC + trait:jdateMC + at.level(trait,"s"):birthjdateMC + at.level(trait,"r"):CapAgeYears + at.level(trait,"s"):MumAge + at.level(trait,"r"):time, random=~idh(at.level(trait,"r")):id, rcov=~idh(at.level(trait,"s")):id + idh(at.level(trait,"r")):id:time, family = NULL, data=femalesC, prior=prior2,nitt=66000,thin=50,burnin=3000)
save(modelCovuLHLLfemales,file="phenCovP_LHLL_AHLL_females_130125.RData")

####################
# Male horn length #
####################
# Load data
data <- read.csv(paste(file_path_data,"DataPhenoCovu_AMHL_LMHL_RandomIds_May2025.csv",sep=""))

# Format
data$CapAgeYears <- as.factor(data$CapAgeYears)
data$trait <- as.factor(data$trait)
data$id <- as.factor(data$id)
summary(data)

# Bivariate models of lamb male horn length with later-life male horn length
prior2 <- list(R = list(R1=list(V=diag(2), nu=2.002, covu=TRUE), R2=list(V=1, nu=2.002)))
modelCovuLMHL <- MCMCglmm(y ~ trait - 1 + trait:Twin + trait:PopDensMC + trait:jdateMC + at.level(trait,"s"):birthjdateMC + at.level(trait,"r"):CapAgeYears + at.level(trait,"s"):MumAge + at.level(trait,"r"):time, random=~idh(at.level(trait,"r")):id, rcov=~idh(at.level(trait,"s")):id + idh(at.level(trait,"r")):id:time, family = NULL, data=data, prior=prior2,nitt=36000,thin=20,burnin=3000)
save(modelCovuLMHL,file="phenCovP_LMHL_AMHL_130125.RData")

##############################
# Male scrotal circumference #
##############################
# Load data
data <- read.csv(paste(file_path_data,"DataPhenoCovu_AMSC_LMSC_RandomIds_May2025.csv",sep=""))

# Format
data$CapAgeYears <- as.factor(data$CapAgeYears)
data$trait <- as.factor(data$trait)
data$id <- as.factor(data$id)
summary(data)

# Bivariate models of lamb male scrotal circumference with later-life male scrotal circumference
prior2 <- list(R = list(R1=list(V=diag(2), nu=2.002, covu=TRUE), R2=list(V=1, nu=2.002)))
modelCovuLMSC <- MCMCglmm(y ~ trait - 1 + trait:Twin + trait:PopDensMC + trait:jdateMC + at.level(trait,"s"):birthjdateMC  + at.level(trait,"r"):CapAgeYears + at.level(trait,"s"):MumAge + at.level(trait,"r"):time, random=~idh(at.level(trait,"r")):id, rcov=~idh(at.level(trait,"s")):id + idh(at.level(trait,"r")):id:time, family = NULL, data=data, prior=prior2,nitt=56000,thin=50,burnin=3000)
save(modelCovuLMSC,file="phenCovP_LMSC_AMSC_130125.RData")


###########################################################
### Phenotypic covariance among lamb august size traits ###
###########################################################
# Load data
data <- read.csv(paste(file_path_data,"DataPhenoCovarLambTraits_RandomIds_May2025.csv",sep=""))

# Run multi-response model
# Both sexes
lambModel <- MCMCglmm(c(LambAugWeight, LambHindLeg, LambHornLen, LambBolCirc) ~ at.level(trait,c("LambAugWeight","LambHindLeg")):Sex + trait:PopDensBirth + trait:MumAge + trait:Twin + at.level(trait,c("LambAugWeight","LambHindLeg","LambBolCirc")):CapJDate, rcov=~us(trait):units, family=c("gaussian","gaussian","gaussian","gaussian"), data=data, nitt=103000, thin=10, burnin=3000)
save(lambModel,file="PmatrixLambTraits_130525.RData")

# Get selection gradients
# P matrix lamb traits
load(file=paste(file_path_output,"PmatrixLambTraits_130525.RData",sep=""))
# Selection differentials viability lamb traits
load(file=paste(file_path_output,"ViabResDataScaleLAW.RData",sep=""))
load(file=paste(file_path_output,"ViabResDataScaleLHL.RData",sep=""))
load(file=paste(file_path_output,"ViabResDataScaleLMHL.RData",sep=""))
load(file=paste(file_path_output,"ViabResDataScaleLMSC.RData",sep=""))

nMCMC <- 10000
lambSelGrads <- NULL
for(i in 1:nMCMC){
	lambP <-matrix(c(lambModel$VCV[i,"traitLambAugWeight:traitLambAugWeight.units"],lambModel$VCV[i,"traitLambHindLeg:traitLambAugWeight.units"],lambModel$VCV[i,"traitLambHornLen:traitLambAugWeight.units"],lambModel$VCV[i,"traitLambBolCirc:traitLambAugWeight.units"],lambModel$VCV[i,"traitLambAugWeight:traitLambHindLeg.units"],lambModel$VCV[i,"traitLambHindLeg:traitLambHindLeg.units"],lambModel$VCV[i,"traitLambHornLen:traitLambHindLeg.units"],lambModel$VCV[i,"traitLambBolCirc:traitLambHindLeg.units"],lambModel$VCV[i,"traitLambAugWeight:traitLambHornLen.units"],lambModel$VCV[i,"traitLambHindLeg:traitLambHornLen.units"],lambModel$VCV[i,"traitLambHornLen:traitLambHornLen.units"],lambModel$VCV[i,"traitLambBolCirc:traitLambHornLen.units"],lambModel$VCV[i,"traitLambAugWeight:traitLambBolCirc.units"],lambModel$VCV[i,"traitLambHindLeg:traitLambBolCirc.units"],lambModel$VCV[i,"traitLambHornLen:traitLambBolCirc.units"],lambModel$VCV[i,"traitLambBolCirc:traitLambBolCirc.units"]),ncol=4)
	lambP

	lambS <- c(ViabResDataScaleLAW$Sz_approx[i]/sqrt(lambModel$VCV[i,"traitLambAugWeight:traitLambAugWeight.units"]), ViabResDataScaleLHL$Sz_approx[i]/sqrt(lambModel$VCV[i,"traitLambHindLeg:traitLambHindLeg.units"]), ViabResDataScaleLMHL$Sz_approx[i]/sqrt(lambModel$VCV[i,"traitLambHornLen:traitLambHornLen.units"]), ViabResDataScaleLMSC$Sz_approx[i]/sqrt(lambModel$VCV[i,"traitLambBolCirc:traitLambBolCirc.units"]))

	beta <- solve(lambP) %*% lambS
	
	lambSelGrads <- cbind(lambSelGrads,beta)
}
lambSelGrads <- as.data.frame(t(lambSelGrads))
names(lambSelGrads) <- c("LambAugustWeight","LambHindlegLength","MaleLambHornLength","LambScrotalCirc")
summary(lambSelGrads)
save(lambSelGrads,file=paste(file_path_models,"lambSelGrads_130525.RData",sep=""))

# Get standardised selection gradients
nMCMC <- 10000
lambSelGrads_stand <- NULL
for(i in 1:nMCMC){
	lambP <-matrix(c(lambModel$VCV[i,"traitLambAugWeight:traitLambAugWeight.units"],lambModel$VCV[i,"traitLambHindLeg:traitLambAugWeight.units"],lambModel$VCV[i,"traitLambHornLen:traitLambAugWeight.units"],lambModel$VCV[i,"traitLambBolCirc:traitLambAugWeight.units"],lambModel$VCV[i,"traitLambAugWeight:traitLambHindLeg.units"],lambModel$VCV[i,"traitLambHindLeg:traitLambHindLeg.units"],lambModel$VCV[i,"traitLambHornLen:traitLambHindLeg.units"],lambModel$VCV[i,"traitLambBolCirc:traitLambHindLeg.units"],lambModel$VCV[i,"traitLambAugWeight:traitLambHornLen.units"],lambModel$VCV[i,"traitLambHindLeg:traitLambHornLen.units"],lambModel$VCV[i,"traitLambHornLen:traitLambHornLen.units"],lambModel$VCV[i,"traitLambBolCirc:traitLambHornLen.units"],lambModel$VCV[i,"traitLambAugWeight:traitLambBolCirc.units"],lambModel$VCV[i,"traitLambHindLeg:traitLambBolCirc.units"],lambModel$VCV[i,"traitLambHornLen:traitLambBolCirc.units"],lambModel$VCV[i,"traitLambBolCirc:traitLambBolCirc.units"]),ncol=4)
	lambP

	lambS <- c(ViabResDataScaleLAW$Sz_approx[i]/sqrt(lambModel$VCV[i,"traitLambAugWeight:traitLambAugWeight.units"]), ViabResDataScaleLHL$Sz_approx[i]/sqrt(lambModel$VCV[i,"traitLambHindLeg:traitLambHindLeg.units"]), ViabResDataScaleLMHL$Sz_approx[i]/sqrt(lambModel$VCV[i,"traitLambHornLen:traitLambHornLen.units"]), ViabResDataScaleLMSC$Sz_approx[i]/sqrt(lambModel$VCV[i,"traitLambBolCirc:traitLambBolCirc.units"]))
	
	beta_stand <- solve(cov2cor(lambP))%*% lambS

	lambSelGrads_stand <- cbind(lambSelGrads_stand,beta_stand)
}
lambSelGrads_stand <- as.data.frame(t(lambSelGrads_stand))
names(lambSelGrads_stand) <- c("LambAugustWeight","LambHindlegLength","MaleLambHornLength","LambScrotalCirc")
summary(lambSelGrads_stand)
save(lambSelGrads_stand,file=paste(file_path_models,"lambSelGrads_stand_130525.RData",sep=""))

# Males
maleLambs <- data[which(data$Sex=="2"),]
# P-matrix for lamb traits
maleLambModel <- MCMCglmm(c(LambAugWeight, LambHindLeg, LambHornLen, LambBolCirc) ~ trait:PopDensBirth + trait:MumAge + trait:Twin + at.level(trait,c("LambAugWeight","LambHindLeg","LambBolCirc")):CapJDate, rcov=~us(trait):units, family=c("gaussian","gaussian","gaussian","gaussian"), data= maleLambs, nitt=1003000, thin=100, burnin=3000)
save(maleLambModel,file="PmatrixMaleLambTraits_270125.RData") 

# Get selection gradients
# P matrix male lamb traits
load(file=paste(file_path_models,"PmatrixMaleLambTraits_270125.RData",sep="")) #10000
# Selection differentials viability lamb traits
load(file=paste(file_path_models,"ViabResDataScaleLAW_Males.RData",sep="")) #10000
load(file=paste(file_path_models,"ViabResDataScaleLHL_Males.RData",sep="")) #10000
load(file=paste(file_path_models,"ViabResDataScaleLMHL.RData",sep="")) #10000
load(file=paste(file_path_models,"ViabResDataScaleLMSC.RData",sep="")) #10000

lambModel <- maleLambModel_notMC
nMCMC <- 10000
malelambSelGrads <- NULL
for(i in 1:nMCMC){
	lambP <-matrix(c(lambModel$VCV[i,"traitLambAugWeight:traitLambAugWeight.units"],lambModel$VCV[i,"traitLambHindLeg:traitLambAugWeight.units"],lambModel$VCV[i,"traitLambHornLen:traitLambAugWeight.units"],lambModel$VCV[i,"traitLambBolCirc:traitLambAugWeight.units"],lambModel$VCV[i,"traitLambAugWeight:traitLambHindLeg.units"],lambModel$VCV[i,"traitLambHindLeg:traitLambHindLeg.units"],lambModel$VCV[i,"traitLambHornLen:traitLambHindLeg.units"],lambModel$VCV[i,"traitLambBolCirc:traitLambHindLeg.units"],lambModel$VCV[i,"traitLambAugWeight:traitLambHornLen.units"],lambModel$VCV[i,"traitLambHindLeg:traitLambHornLen.units"],lambModel$VCV[i,"traitLambHornLen:traitLambHornLen.units"],lambModel$VCV[i,"traitLambBolCirc:traitLambHornLen.units"],lambModel$VCV[i,"traitLambAugWeight:traitLambBolCirc.units"],lambModel$VCV[i,"traitLambHindLeg:traitLambBolCirc.units"],lambModel$VCV[i,"traitLambHornLen:traitLambBolCirc.units"],lambModel$VCV[i,"traitLambBolCirc:traitLambBolCirc.units"]),ncol=4)
	lambP

	lambS <- c(ViabResDataScaleLAWM$Sz_approx[i]/sqrt(lambModel$VCV[i,"traitLambAugWeight:traitLambAugWeight.units"]), ViabResDataScaleLHLM$Sz_approx[i]/sqrt(lambModel$VCV[i,"traitLambHindLeg:traitLambHindLeg.units"]), ViabResDataScaleLMHL$Sz_approx[i]/sqrt(lambModel$VCV[i,"traitLambHornLen:traitLambHornLen.units"]), ViabResDataScaleLMSC$Sz_approx[i]/sqrt(lambModel$VCV[i,"traitLambBolCirc:traitLambBolCirc.units"]))

	beta <- solve(lambP) %*% lambS
	
	malelambSelGrads <- cbind(malelambSelGrads,beta)
}
malelambSelGrads <- as.data.frame(t(malelambSelGrads))
names(malelambSelGrads) <- c("LambAugustWeight","LambHindlegLength","MaleLambHornLength","LambScrotalCirc")
summary(malelambSelGrads)
malelambSelGrad_notMC <- malelambSelGrads
save(malelambSelGrad_notMC,file=paste(file_path_models,"MalelambSelGrads_030225.RData",sep=""))

# Get standardised selection gradients
lambModel <- maleLambModel_notMC
nMCMC <- 10000
malelambSelGrads_stand <- NULL
for(i in 1:nMCMC){
	lambP <-matrix(c(lambModel$VCV[i,"traitLambAugWeight:traitLambAugWeight.units"],lambModel$VCV[i,"traitLambHindLeg:traitLambAugWeight.units"],lambModel$VCV[i,"traitLambHornLen:traitLambAugWeight.units"],lambModel$VCV[i,"traitLambBolCirc:traitLambAugWeight.units"],lambModel$VCV[i,"traitLambAugWeight:traitLambHindLeg.units"],lambModel$VCV[i,"traitLambHindLeg:traitLambHindLeg.units"],lambModel$VCV[i,"traitLambHornLen:traitLambHindLeg.units"],lambModel$VCV[i,"traitLambBolCirc:traitLambHindLeg.units"],lambModel$VCV[i,"traitLambAugWeight:traitLambHornLen.units"],lambModel$VCV[i,"traitLambHindLeg:traitLambHornLen.units"],lambModel$VCV[i,"traitLambHornLen:traitLambHornLen.units"],lambModel$VCV[i,"traitLambBolCirc:traitLambHornLen.units"],lambModel$VCV[i,"traitLambAugWeight:traitLambBolCirc.units"],lambModel$VCV[i,"traitLambHindLeg:traitLambBolCirc.units"],lambModel$VCV[i,"traitLambHornLen:traitLambBolCirc.units"],lambModel$VCV[i,"traitLambBolCirc:traitLambBolCirc.units"]),ncol=4)
	lambP

	lambS <- c(ViabResDataScaleLAWM$Sz_approx[i]/sqrt(lambModel$VCV[i,"traitLambAugWeight:traitLambAugWeight.units"]), ViabResDataScaleLHLM$Sz_approx[i]/sqrt(lambModel$VCV[i,"traitLambHindLeg:traitLambHindLeg.units"]), ViabResDataScaleLMHL$Sz_approx[i]/sqrt(lambModel$VCV[i,"traitLambHornLen:traitLambHornLen.units"]), ViabResDataScaleLMSC$Sz_approx[i]/sqrt(lambModel$VCV[i,"traitLambBolCirc:traitLambBolCirc.units"]))

	beta_stand <- solve(cov2cor(lambP))%*% lambS

	malelambSelGrads_stand <- cbind(malelambSelGrads_stand,beta_stand)
}
malelambSelGrads_stand <- as.data.frame(t(malelambSelGrads_stand))
names(malelambSelGrads_stand) <- c("LambAugustWeight","LambHindlegLength","MaleLambHornLength","LambScrotalCirc")
summary(malelambSelGrads_stand)
save(malelambSelGrads_stand,file=paste(file_path_models,"MalelambSelGrads_stand_270125.RData",sep=""))

# Females
femaleLambs <- data[which(data$Sex=="1"),]
# P-matrix for lamb traits
femaleLambModel <- MCMCglmm(c(LambAugWeight, LambHindLeg) ~ trait:PopDensBirth + trait:MumAge + trait:Twin + trait:CapJDate, rcov=~us(trait):units, family=c("gaussian","gaussian"), data= femaleLambs, nitt=1003000, thin=100, burnin=3000)
save(femaleLambModel,file="PmatrixFemaleLambTraits_270125.RData")

# Get selection gradients
# P matrix male lamb traits
load(file=paste(file_path_models,"PmatrixFemaleLambTraits_270125.RData",sep=""))
# Selection differentials viability lamb traits
load(file=paste(file_path_models,"ViabResDataScaleLAW_Females.RData",sep="")) 
load(file=paste(file_path_models,"ViabResDataScaleLHL_Females.RData",sep="")) 

nMCMC <- 10000
femalelambSelGrads_stand <- NULL
for(i in 1:nMCMC){
	lambP <-matrix(c(femaleLambModel$VCV[i,"traitLambAugWeight:traitLambAugWeight.units"], femaleLambModel$VCV[i,"traitLambHindLeg:traitLambAugWeight.units"], femaleLambModel$VCV[i,"traitLambAugWeight:traitLambHindLeg.units"], femaleLambModel$VCV[i,"traitLambHindLeg:traitLambHindLeg.units"]),ncol=2)
	lambP

	lambS <- c(ViabResDataScaleLAWF$Sz_approx[i]/sqrt(femaleLambModel$VCV[i,"traitLambAugWeight:traitLambAugWeight.units"]), ViabResDataScaleLHLF$Sz_approx[i]/sqrt(femaleLambModel$VCV[i,"traitLambHindLeg:traitLambHindLeg.units"]))
	
	beta_stand <- solve(cov2cor(lambP))%*% lambS

	femalelambSelGrads_stand <- cbind(femalelambSelGrads_stand,beta_stand)
}
femalelambSelGrads_stand <- as.data.frame(t(femalelambSelGrads_stand))
names(femalelambSelGrads_stand) <- c("LambAugustWeight","LambHindlegLength")
summary(femalelambSelGrads_stand)
save(femalelambSelGrads_stand,file=paste(file_path_models,"FemalelambSelGrads_270125.RData",sep=""))

# Get standardised selection gradients
nMCMC <- 10000
femalelambSelGrads <- NULL
for(i in 1:nMCMC){
	lambP <-matrix(c(femaleLambModel$VCV[i,"traitLambAugWeight:traitLambAugWeight.units"], femaleLambModel$VCV[i,"traitLambHindLeg:traitLambAugWeight.units"], femaleLambModel$VCV[i,"traitLambAugWeight:traitLambHindLeg.units"], femaleLambModel$VCV[i,"traitLambHindLeg:traitLambHindLeg.units"]),ncol=2)
	lambP

	lambS <- c(ViabResDataScaleLAWF$Sz_approx[i]/sqrt(femaleLambModel$VCV[i,"traitLambAugWeight:traitLambAugWeight.units"]), ViabResDataScaleLHLF$Sz_approx[i]/sqrt(femaleLambModel$VCV[i,"traitLambHindLeg:traitLambHindLeg.units"]))

	beta <- solve(lambP) %*% lambS
	
	femalelambSelGrads <- cbind(femalelambSelGrads,beta)
}
femalelambSelGrads <- as.data.frame(t(femalelambSelGrads))
names(femalelambSelGrads) <- c("LambAugustWeight","LambHindlegLength")
summary(femalelambSelGrads)
save(femalelambSelGrads,file=paste(file_path_models,"FemalelambSelGrads_030225.RData",sep=""))









