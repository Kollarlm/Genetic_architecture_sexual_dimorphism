library(MCMCglmm) 
library(tidyr) 
library(plyr) 
library(dplyr)
library(MuMIn)
library(ggplot2)
library(car)
library(gdata)
library(matrixcalc)
library(coda)
library(ellipse)
library(MasterBayes)
library(corrplot)
library(parallel)
library(coda)
library(ggpubr)
library(wesanderson)


##############################################################
# Creating comparable variable across sexes for reproduction #
##############################################################

# Dataframes
#data <- read.csv("~/PATH/Data_for_models/LH.traits.data.NONA.csv")
data <- read.csv("~/PATH/data.homemade.traits.csv")

###### For the growth and developmental dataset ######

# Formatting data for growth and development 
data$X <- NULL
data$famid <- as.factor(data$famid)
data$sampid <- as.factor(data$sampid) 
data$ssex <- as.factor(data$ssex)
data$Plate <- as.factor(data$Plate)
data$Days_gam <- as.numeric(data$Days_gam)
data$Total_gam <- as.numeric(data$Total_gam)

# Scale growth and development traits
data$area_wk3 <- scale(data$area_wk3)
data$perim_wk3 <- scale(data$perim_wk3)
data$circ_wk3 <- scale(data$circ_wk3)
data$Total_gam <- scale(data$Total_gam)

# Formatting data for morpholohy and physiology
data$X <- NULL
data$famid <- as.factor(data$famid)
data$sampid <- as.factor(data$sampid) 
data$ssex <- as.factor(data$ssex)
data$Date.PTR <- as.factor(data$Date.PTR)

###### For the morphology and physiology dataset ######

# Log transform VOC traits
data$Total.conc <- log10(data$Total.conc + 1e-10)
data$Avg.conc.raw <- log10(data$Avg.conc.raw + 1e-10)
data$Total.conc.stand.mean <- log10(data$Total.conc.stand.mean + 1e-10)
data$Avg.conc.divmean <- log10(data$Avg.conc.divmean + 1e-10)

# Scaling the reproductive traits seperately for males and females
data.f <- subset(data, data$ssex == "f")
data.f <- droplevels(data.f)
data.m <- subset(data, data$ssex == "m")
data.m <- droplevels(data.m)
data.f$rel_repro <- data.f$data.repro / mean(data.f$data.repro) # mean(f$rel_repro) = 1
data.m$rel_repro <- data.m$data.repro / mean(data.m$data.repro) # mean(m$rel_repro) = 1
data <- rbind(data.f, data.m)

# Scale morphology and physiology traits
data$Avg.conc.raw <- scale(data$Avg.conc.raw)
data$Total.conc <- scale(data$Total.conc)
data$Total.masses <- scale(data$Total.masses)
data$data.Leaf_Length_Average <- scale(data$data.Leaf_Length_Average)
data$rel_repro <- scale(data$rel_repro)
######################################################################################################################################

#####################################
## MCMC simulation of the G matrix ##
#####################################

# # Adjust "gaus", "pois", and "prior" accordingly

n=4 # dim of Gm and Gf
Blength=n^2 # length of the B matrix
m=2 # Number of matrices in an array
Gn=n*m #
Glength=Gn^2 # length of the G matrix 
Gnames=c("female","male")

## Error structure vector
gaus <- rep("gaussian", n) 

prior2 <- list(R=list(V=diag(0.01, Gn), nu=7.001),
               G = list(G1 = list(V=diag(0.001, Gn), nu=7.001),
                        G2 = list(V = 0.001, nu = 7.001),
                        G3 = list(V = 0.001, nu = 7.001)))


parexp <- list(R = list(V = diag(0.01, Gn), nu = 0.002),
               G = list(G1 = list(V = diag(0.001, Gn), nu = Gn, alpha.mu=rep(0, Gn), alpha.V = diag(1000,Gn)),
                        G2 = list(V = 0.001, nu = 0.002, alpha.mu=0, alpha.v = 1000),
                        G3 = list(V = 0.001, nu = 0.002, alpha.mu=0, alpha.v = 1000)))


## Running the models for 1.) Growth and development
##                        2.) Morphology and physiology

## Adding in the relative reproduction traits
Final.VOC.4.repro <- MCMCglmm(c(Total.masses, Total.conc.stand.mean, data.Leaf_Length_Average, rel_repro)  ~
                                trait-1 + trait:ssex, # Fixed effects
                              random = ~us(trait:ssex):famid + sampid + Date.PTR, # Random effects (us() estimates covariances)
                              rcov = ~idh(trait:ssex):units, # residual matrix, see Haddfield's course notes, top of pg. 71
                              data = data,
                              family = gaus,
                              prior = prior2,
                              nitt=1005000,
                              burnin=5000,
                              thin=1000,
                              slice=TRUE,
                              verbose = T,
                              pr=T)

plot(Final.VOC.4.repro)
save(Final.VOC.4.repro, file = "Final.VOC.4.repro.obj")

#LH trait model with only traits with additive genetic variance
scaled.ThreeLHtraits.parexp <- MCMCglmm(cbind(perim_wk3, circ_wk3, Total_gam) ~
                                          trait-1 + trait:ssex, # Fixed effects
                                        random = ~us(trait:ssex):famid + # Random effects (us() estimates covariances)
                                          sampid +
                                          Plate,
                                        rcov = ~idh(trait:ssex):units, # no residual covariance between sexes because those measure stem from different individuals, see Haddfield's course notes, top of pg. 71
                                        data = data,
                                        family = gaus,
                                        prior = parexp, # parameter expanded
                                        nitt=1005000,
                                        burnin=5000,
                                        thin=1000,
                                        slice=TRUE,
                                        verbose = T,
                                        pr=T)


plot(scaled.ThreeLHtraits.parexp)
save(scaled.ThreeLHtraits.parexp, file = "scaled.ThreeLHtraits.parexp.obj") 

## Trace plots and posterior distributions
par(mar=c(1,1,1,1))
plot(scaled.TwoLHtraits.parexp$VCV)

## Autocorrlation plots
autocorr.plot(scaled.TwoLHtraits.parexp$VCV)
autocorr(scaled.TwoLHtraits.parexp$VCV)

######################################################################################################################################

############################
## Gelman Rubin Criterion ##
############################

## Replace model for 
set.seed(1)
scaled.threeLHtraits.parexp_GRC <- mclapply(1:4, function(i) {
  MCMCglmm(cbind(perim_wk3, circ_wk3, Total_gam) ~
             trait-1 + trait:ssex, # Fixed effects
           random = ~us(trait:ssex):famid + # Random effects (us() estimates covariances)
             sampid,
           rcov = ~idh(trait:ssex):units, # no residual covariance between sexes because those measure stem from different individuals, see Haddfield's course notes, top of pg. 71
           data = data,
           family = gaus,
           prior = parexp, # parameter expanded
           nitt=1005,
           burnin=5,
           thin=1,
           slice=TRUE,
           verbose = T,
           pr=T)
}, mc.cores = 4)


scaled.threeLHtraits.parexp_GRC <- lapply(scaled.threeLHtraits.parexp_GRC, function(m) m$Sol)
scaled.threeLHtraits.parexp_GRC <- do.call(mcmc.list, scaled.threeLHtraits.parexp_GRC)

## Plot 
par(mfrow=c(4,2), mar=c(2,2,1,2))
pdf("scaled.threeLHtraits.parexp_GRC.pdf")
gelman.plot(scaled.threeLHtraits.parexp_GRC, auto.layout=F)
dev.off()

gelman.diag(scaled.threeLHtraits.parexp_GRC)

## Plot 
par(mfrow=c(8,2), mar=c(2, 1, 1, 1))
pdf("scaled.threeLHtraits.parexp_GRC2.pdf")
plot(scaled.threeLHtraits.parexp_GRC, ask=F, auto.layout=F)
dev.off()

####################
## Loading models ##
####################

load(file ="~/PATH/scaled.threeLHtraits.parexp.obj") # LH Traits
summary(scaled.ThreeLHtraits.parexp)

load(file ="~/PATH/Final.VOC.4.repro.obj") # VOC traits
summary(Final.VOC.4.repro)

##################################################### 
## Setting all paramters for submatrix comparisons ##
##################################################### 

### Setting parameters for the code ###
## Change suffix and model for LH traits and VOCs
n=4 # dim of Gm and Gf. Use #3 for Growth and development model
Blength=n^2 # Length of the B matrix
m=2 # Number of matrices in an array
Gn=n*m 
Glength=Gn^2 # Length of the G matrix
model <- Final.VOC.4.repro # Morphology and physiology
suffix <- "Final.VOC.4.repro" # Morphology and physiology
# model <- scaled.ThreeLHtraits.parexp # Growth and Development
# suffix <- "scaled.ThreeLHtraits.parexp" # Growth and Development
o <- 1 #Starting value for extracting matrix
p <- n+1 # Starting value for extracting matrix (n + 1)
q <- Gn # Ending value for extracting matrix
#traitnames = c("Juvenile growth (P)", "Juvenille growth form", "Mature Tissue") # Growth and Development
traitnames = c("Total masses", "Total conc", "Leaf Length", "Relative repro") # Morphology and physiology

##############################################################################################
# Get genetic variances (diag(G.cov)), rMFs (diag(B.cor)), and resample the correlation of rMFs~SDs #
##############################################################################################

# Get covariance matrix (G.cov) and HPD intervals
G.cov <- matrix(colMeans(model$VCV[,1:Glength]),nrow=Gn,ncol=Gn,byrow=F)

# Get posterior 95% HPD intervals for G.cov
# By applying this function to Gmcmc 
Gmat_HPD <-  function (Gmat) {
  corG2 <- lapply(Gmat, as.vector)
  corG2 <- do.call(rbind,corG2)
  corG2 <- as.mcmc(corG2)
  fhpd <- function (x) { x <- as.mcmc(x); HPDinterval(x, prob=0.95) }
  hpd <- round(apply(corG2,2, fhpd ),3)
  int.res <- paste(hpd[1,],hpd[2,],sep=" , ")
  mat.int <- matrix(int.res,nrow=Gn,ncol=Gn, byrow=F)
  return(mat.int)
}

# Extract 1000 stored G matrices from model. 
f <- function(x){list(matrix(x,nrow=Gn,ncol=Gn,byrow=F))}
Gmcmc <- lapply(apply(model$VCV[,1:Glength], 1, f), "[[", 1) 

# Apply Gmat_HPD function to Gmcmc to get the HPD intervals for G.cov
CI_G.cov <- Gmat_HPD(Gmcmc) 

# Genetic variances (Vg, which I think = Va in this case (they're basically "isofemale lines"))
fem.Va <- diag(G.cov)[o:n]
male.Va <- diag(G.cov)[p:q]

# and CIs
fem.CIs <- diag(CI_G.cov)[o:n]
male.CIs <- diag(CI_G.cov)[p:q]

# Get correlation matrix (G.cor)
G.cor <- matrix(posterior.mode(posterior.cor(model$VCV[,1:Glength])),Gn,Gn)

# With HPDintervals for the correlation coefficients
CI_G.cor <- HPDinterval(posterior.cor(model$VCV[,1:Glength]))

# Put those CIs into a matrix
G.cor_CIs <- 1:Glength
for (i in 1:Glength) {
  G.cor_CIs[i] <- paste(CI_G.cor[i,1], CI_G.cor[i,2], sep=" , ")
}
CI_G.cor <- matrix(G.cor_CIs[1:Glength], nrow=Gn, ncol=Gn, byrow=F)

# Get correlation matrix for B (B.cor)
B.cor <- G.cor[o:n,p:q] 

# Get HPD intervals for B.cor
CI_B.cor <- CI_G.cor[o:n,p:q] 

##################################################################################
# Intersexual genetic correlations (rMFs) and their HPD intervals for each trait #
##################################################################################

rMF <- diag(B.cor)
CIs_rMF <- diag(CI_B.cor)

# Sexual dimorphism estimates 
sex.dim <- colMeans(model$Sol[,p:q])


# A summary of the male and female genetic variances, sexual dimorphism and rmf.
Table_one <- as.data.frame(cbind(traitnames, fem.Va, fem.CIs, male.Va, male.CIs, rMF, CIs_rMF, sex.dim))

# Write table 1
write.csv(Table_one, paste0("~/Desktop/Quantitative_Genetics_Moss2/Figures_Resubmission/Table1_", suffix, ".csv"))


#################################################################################################################################

# Comparing Gm and Gf using multiple methods #
# 1.) Hansens difference (d)
# 2.) Eigenstensor analysis
# 3.) Simplified eigenstensor analysis
# Code has been adapted from Puentes et al 2016 and Aguirre et al 2014

##################################
#  1.)  Hanson's difference (d)  #
##################################
# Code adapted from Puentes et al 2016

# Median estimate of G.cov
Gmcmc.median <- matrix(apply(model$VCV[,1:Glength],2,median),nrow=Gn,ncol=Gn,byrow=F)
Gf <- Gmcmc.median[o:n,o:n] 
Gm <- Gmcmc.median[p:q,p:q] 

# Point estimate of Hanson's d between male and female sub-matrices Gf and Gm 
d_M.F <- sqrt(mean(eigen(Gf-Gm)$values^2)) * (1 - (var(eigen(Gf-Gm)$values^2)/mean(eigen(Gf-Gm)$values^2)^2) / (4*(length(eigen(Gf-Gm)$values)+2) ) )#r
d_table <- data.frame("Sex"= c("Gf-Gm"), 
                      "Response d"= round(c(d_M.F), 8), "lo.ci"=NA, "up.ci"=NA)


## Check uncertainty of differences and test if d != 0

# # Extract 1000 stored  sub-matrices for... 
f.sub <- function(x){list(matrix(x,nrow=n,ncol=n,byrow=F))} # n*n sub-matrix function

# Gf
Gf.mcmc <- list()
for (i in 1:1000) {
  Gf.mcmc[[i]] <- Gmcmc[[i]][o:n, o:n]
}

# Gm
Gm.mcmc <- list()
for (i in 1:1000) {
  Gm.mcmc[[i]] <- Gmcmc[[i]][p:q, p:q]
}

sex.Gs <- list(Gf.mcmc, Gm.mcmc)
MCMCsamp=length(Gf.mcmc) 
Garray.sex <- array(,c(n,n,m,MCMCsamp))
Gnames=c("female","male")
dimnames(Garray.sex) <- list(traitnames,traitnames,Gnames)
for (k in 1:m) {
  for (i in 1:MCMCsamp) {
    Garray.sex[,,k,i]<-sex.Gs[[k]][[i]]
  }
}


for (k in 1:2) {
  Garray <- Garray.sex
  if(k  == 1) {sex1 <- 1; sex2 <- 2} 
  
  yy <- numeric(1000)  # this is the empty vector for how many times you want to resample d
  for (i in 1:1000) {
    nr1   <- sample(1:1000,1) # this picks a random number between 1:1000 
    nr2 	<- sample(1:1000,1) # same
    G11      <- Garray[,, sex1,nr1] # the random number is used to sample a random one of the 1000 matrices in Garray.sex that correspond to pop1, i.e. k == 1, i.e. female matrices
    G12      <- Garray[,, sex1,nr2] # and then again.
    G21      <- Garray[,, sex2,nr1] # Then same, but for pop2, i.e. k == 2, i.e. male matrices
    G22      <- Garray[,, sex2,nr2] # and again. So that below the do something like (d+d)-(d+d); see annotations at end of those lines
    
    yy[i] <- ( (sqrt(mean(eigen(G11-G12)$values^2))*(1 - (var(eigen(G11-G12)$values^2)/mean(eigen(G11-G12)$values^2)^2) / (4*(length(eigen(G11-G12)$values)+2) ) )) + # Hanson's d between two randomly sampled female posterior matrix
                 (sqrt(mean(eigen(G21-G22)$values^2))*(1 - (var(eigen(G21-G22)$values^2)/mean(eigen(G21-G22)$values^2)^2) / (4*(length(eigen(G21-G22)$values)+2) ) ))) - # summed with the d between two randomly sampled male matrix
      ((sqrt(mean(eigen(G11-G22)$values^2))*(1 - (var(eigen(G11-G22)$values^2)/mean(eigen(G11-G22)$values^2)^2) / (4*(length(eigen(G11-G22)$values)+2) ) )) + # minus d between a random male and random female matrix
         (sqrt(mean(eigen(G12-G21)$values^2))*(1 - (var(eigen(G12-G21)$values^2)/mean(eigen(G12-G21)$values^2)^2) / (4*(length(eigen(G12-G21)$values)+2) ) ))) # plus another calculation of d between a random male and random female matrix.
    # The difference between the sum of the d estimates in the first two lines minus the sum of the d estimates in the last two lines is estimated many times in order to get an HPD interval around that difference.
    
  }
  d_table[ k, c("lo.ci","up.ci")] <- round(HPDinterval(as.mcmc(yy), prob = 0.95)[ c(2, 1)],8)*-1 
  # *-1 to get the signs correct since d is positive  
}

# Print Table 
print(d_table) 

################################################
# 2.) Eigaen tensor analysis between Gf and Gm #
################################################
# Code adapted from Aguirre et al 2014

#START
covtensor <- function(Gs){
  if (dim(Gs)[[1]] != dim(Gs)[[2]]){
    stop("G array must be of order n x n x m x MCMCsamp")
  }
  if (is.na(dim(Gs)[4])) {
    stop("There are no MCMCsamples")
  }
  neigten <- n*(n+1)/2 
  #Number of eigentensors
  MCMC.S <- array(,c(neigten, neigten, MCMCsamp))
  dimnames(MCMC.S) <- list(paste("e", 1:neigten, sep=""), paste("e", 1:neigten, sep=""))
  for (k in 1:MCMCsamp){
    MCMCG <- Gs[,,,k] 
    MCMCvarmat <- t(apply(MCMCG, 3, diag)) 
    #find the variances of the kth G and store them 
    MCMCcovmat <- t(apply(MCMCG, 3, lowerTriangle)) 
    #find the covariances of the kth G and store them
    MCMC.S[1:n,1:n, k] <- cov(MCMCvarmat, MCMCvarmat) 
    #fill the upper left quadrant of the kth S
    MCMC.S[(n+1):neigten,(n+1):neigten, k] <- 2*cov(MCMCcovmat, MCMCcovmat)
    #fill the lower right quadrant of the kth S
    MCMC.S[1:n,(n+1):neigten, k] <- sqrt(2)*cov(MCMCvarmat, MCMCcovmat)
    #fill the upper right quadrant of the kth S
    MCMC.S[(n+1):neigten,1:n, k] <- sqrt(2)*cov(MCMCcovmat, MCMCvarmat)
    #fill the lower left quadrant of the kthS
  }  
  av.S <- apply(MCMC.S, 1:2, mean)
  #posterior mean S
  av.S.val <- eigen(av.S)$values
  #eigenvalues of posterior mean S 
  av.S.vec <- eigen(av.S)$vectors
  #eigenvalues of posterior mean S
  eTmat <- array(, c(n, n, neigten))
  dimnames(eTmat) <- list(traitnames, traitnames, paste("E", 1:neigten, sep=""))  
  for (i in 1:neigten){
    emat <- matrix(0, n, n) 
    lowerTriangle(emat) <- 1/sqrt(2)*av.S.vec[(n+1):neigten,i]
    emat <- emat + t(emat)
    diag(emat) <- av.S.vec[1:n,i]
    eTmat[,,i] <- emat 
  }
  #construct the second-order eigentensors of posterior mean S
  eT.eigen <- array(, c(n+1, n, neigten))
  for (i in 1:neigten){
    eT.eigen[1,,i] <- t(eigen(eTmat[,,i])$values) 
    #Eigenvalues of the ith eigentensor
    eT.eigen[2:(n+1),,i] <- eigen(eTmat[,,i])$vectors 
    #Eigenvectors of the ith eigentensor
    eT.eigen[,,i] <- eT.eigen[,order(abs(eT.eigen[1,,i]), decreasing = T), i]
  }
  MCMC.S.val <- matrix(, MCMCsamp, neigten)
  colnames(MCMC.S.val) <- paste("E", 1:neigten, sep="")
  for (i in 1:MCMCsamp){
    for(j in 1:neigten){
      MCMC.S.val[i,j] <- t(av.S.vec[,j]) %*% MCMC.S[,,i] %*% av.S.vec[,j]
    }
  }
  #posterior distribution of the genetic variance for the eigenvectors of posterior mean S
  av.G.coord <- array(, c(m, neigten, 1))
  dimnames(av.G.coord) <- list(Gnames, paste("E", 1:neigten, sep=""))
  for (i in 1:neigten){
    av.G.coord[,i,] <- apply((apply(Gs, 1:3, mean)) , 3, frobenius.prod, y = eTmat[,,i])
  }
  #Coordinates of the jth avG for the eigentensors of posterior mean S
  MCMC.G.coord <- array(, c(m, neigten, MCMCsamp))
  dimnames(MCMC.G.coord) <- list(Gnames, paste("E", 1:neigten, sep=""))
  for (i in 1:neigten){
    MCMC.G.coord[,i,] <- apply(Gs, 3:4, frobenius.prod, y = eTmat[,,i])
  }
  #Coordinates of the kth MCMC sample of the jth G for the eigentensors of posterior mean S
  tensor.summary <- data.frame(rep(av.S.val,each=n), t(data.frame(eT.eigen)))
  colnames(tensor.summary) <- c("S.eigval", "eT.val", traitnames)
  rownames(tensor.summary)<- paste(paste("e", rep(1:neigten, each=n), sep=""), rep(1:n,neigten), sep=".")
  list(tensor.summary = tensor.summary, av.S = av.S, eTmat = eTmat, av.G.coord = av.G.coord, MCMC.S = MCMC.S, MCMC.S.val = MCMC.S.val, MCMC.G.coord = MCMC.G.coord)
}


# Applying "covtensor" function to the observed G array and storing the data as MCMC.covtensor
MCMC.covtensor <- covtensor(Garray.sex) 
#summary(MCMC.covtensor$tensor.summary) #Summary see Aguirre doc for components

# Calcualte the max number of nonzero-eigentensor
nnonzero <- min(n*(n+1)/2,m-1) 

# Generate a null distribution of variance among G
MCMC.covtensor.rand <- covtensor(rand.Garray) 

# Examining posterior means and 95% HPD intervals for nonzero-egentensors for the observed and randomized G arrays.
HPD.eT.val <- cbind(HPDinterval(as.mcmc(MCMC.covtensor$MCMC.S.val[,1:nnonzero]), prob=0.95), 
                    HPDinterval(as.mcmc(MCMC.covtensor.rand$MCMC.S.val[,1:nnonzero]), prob=0.95))
round(HPD.eT.val, 3)

## Plot Eigenvalues of the non-zero eignetensors for posterior mean, S and 95% HPD intervals of the eignvalues of the non-zero eigntensors of S for each MCMC samples of the observed randomized G arrays
pdf(paste0("~/Desktop/Quantitative_Genetics_Moss2/Figures_Resubmission/Eigentensor_", suffix, ".pdf"))
plot((1:1)-0.1,colMeans(MCMC.covtensor$MCMC.S.val)[1:1],type="p",xlab="Genetic Covariance Tensor",ylab="alpha",
     pch=16,cex=1,xaxt="n",frame.plot=F, ylim=c(-.01,.5),xlim=c(0.5,3.5),cex.lab=1.2,las=1)
axis(1,at=1:1,labels=c(paste("E",rep(1:1),sep="")),cex=1.1)
points((1:1)+0.1,colMeans(MCMC.covtensor.rand$MCMC.S.val)[1:1],type="p",pch=1,cex=1)
arrows((1:1)-0.1,colMeans(MCMC.covtensor$MCMC.S.val)[1:1],(1:1)-0.1,HPD.eT.val[1:1,1],length=0.1,angle=90)
arrows((1:1)-0.1,colMeans(MCMC.covtensor$MCMC.S.val)[1:1],(1:1)-0.1,HPD.eT.val[1:1,2],length=0.1,angle=90)
arrows((1:1)+0.1,colMeans(MCMC.covtensor.rand$MCMC.S.val)[1:1],(1:1)+0.1,HPD.eT.val[1:1,3],length=0.1,angle=90,lty=5)
arrows((1:1)+0.1,colMeans(MCMC.covtensor.rand$MCMC.S.val)[1:1],(1:1)+0.1,HPD.eT.val[1:1,4],length=0.1,angle=90,lty=5)
legend("topright",legend=c("observed","randomised"),lty=c(4,1),pch=c(16,1),cex=1.2,bty="n")
dev.off()

trait.combo <- round(MCMC.covtensor$tensor.summary[1:(n*2),2:dim(MCMC.covtensor$tensor.summary)[2]], 3) #Examining first first 10 rows of E1 and E2 to see trait combinations

#############################################
## Simplified difference between Gf and Gm ##
#############################################

# 1000 differneces between Gf and Gm

Ds <- 1:1000

for (i in 1:MCMCsamp) {
  Ds[i] <- sum(eigen(sex.Gs[[1]][[i]] - sex.Gs[[2]][[i]])$values) 
}


# Result
posterior.mode(as.mcmc(Ds), na.rm=T) 
HPDinterval(as.mcmc(Ds)) # 95% HPD intervals

# The CIs of the trace of the difference matrix overlap zero, but we can test this relative to the null
# Sample the 1000 Gfs and Gms so as to randomly swap the sex labels to get new list of 2000 random n*n Gs (1000 for each sex)
# set up the empty lists (2) of blank vectors (to be coerced into matrices in the loop)
rand.Gf <- vector("list", MCMCsamp)
rand.Gm <- vector("list", MCMCsamp)
rand.GfGm <- list(rand.Gf, rand.Gm)

for (k in 1:m) { # for each of the two lists in rand.GfGm
  for (i in 1:MCMCsamp) { # for each of the 1000 empty vectors (to be coerced into a matrix) in each list
    nr1 <- sample(1:2,1) # draw, randomly, a 1 or a 2 (for indexing, randomly, which list to draw from in sex.Gs)
    rand.GfGm[[k]][[i]] <- matrix(sex.Gs[[nr1]][[i]],nrow=n,ncol=n) # take the ith matrix of the randomly male-or-female list of matrices and put it into the kth list
  }
} 

null.Ds <- 1:1000

for (i in 1:MCMCsamp) {
  null.Ds[i] <- sum(eigen(rand.GfGm[[1]][[i]] - rand.GfGm[[2]][[i]])$values) 
}


# Now test whether the distribution of Ds is different from that of the null.Ds
posterior.mode(as.mcmc(null.Ds), na.rm=T) 
HPDinterval(as.mcmc(null.Ds)) # 95% HPD intervals

P.Ds_v_nullDs <- ifelse(null.Ds < posterior.mode(as.mcmc(Ds), na.rm=T), 0, 1) # Knowing that the difference is greater than the null, we call those null.Ds that are < the Ds 0, else 1...

# Bootstrapped P value
(1-sum(P.Ds_v_nullDs)/1000)*2 # ...and then calculate the 2-tailed P value.

##########################################
# Plotting differences between Gf and Gm #
#           OR                           #
#     Ellipse plots                      #
##########################################

# Covariances
fem.G.cov <- G.cov[o:n,o:n] 
male.G.cov <- G.cov[p:q,p:q] 
up.fem.G <- fem.G.cov[upper.tri(fem.G.cov)]
low.fem.G <- fem.G.cov[lower.tri(fem.G.cov)]
up.male.G <- male.G.cov[upper.tri(male.G.cov)]
low.male.G <- male.G.cov[lower.tri(male.G.cov)]

# CIs for covariances
fem.CI_G.cov <- CI_G.cov[o:n,o:n] 
rownames(fem.CI_G.cov) <- traitnames
colnames(fem.CI_G.cov) <- traitnames
male.CI_G.cov <- CI_G.cov[p:q,p:q] 
rownames(male.CI_G.cov) <- traitnames
colnames(male.CI_G.cov) <- traitnames
CI_up.fem.G <- fem.CI_G.cov[upper.tri(fem.CI_G.cov)]
CI_low.fem.G <- fem.CI_G.cov[lower.tri(fem.CI_G.cov)]
CI_up.male.G <- male.CI_G.cov[upper.tri(male.CI_G.cov)]
CI_low.male.G <- male.CI_G.cov[lower.tri(male.CI_G.cov)]

#Figure 1
# Plotting elipses (Code adapted from Puentes et al 2016)
# Run the "myplotcor" function at the end of the script first.
png("corr_plot_VOC.png",type="cairo", units="in", 
    width=5, height=5, pointsize=8, res=300)
sign.mat <- matrix(0,ncol=Gn,nrow=Gn)
par(col = "tomato2",lty=1,lwd=2)
#sign.mat[4,1] <- 1
mat <- cov2cor(fem.G.cov)
colnames(mat)<- traitnames
rownames(mat)<- traitnames
diag(mat) <- NA
myplotcor(mat,diag=T, type="lower",col="tomato2",new=1,sign.mat=sign.mat)
#sign.mat[4,1] <- 0
par(col = "turquoise4", lty=6,lwd=2)
mat <- cov2cor(male.G.cov)
colnames(mat)<- traitnames
rownames(mat)<- traitnames
diag(mat) <- NA
myplotcor(mat,diag=T, type="lower",col="turquoise4",new=0,sign.mat=sign.mat)
legend("topleft", c("Males", "Female"), col=c("turquoise4", "tomato2"), lty =c(6,1), bty="n")
dev.off()

# Elipses with entire G matrix
png("corr_plot_Gmf_Growth.png",type="cairo", units="in", 
    width=8, height=9, pointsize=18, res=372)
sign.mat <- matrix(0,ncol=Gn,nrow=Gn)
par(col = "black",lty=1,lwd=2)
#sign.mat[4,1] <- 1
mat <- cov2cor(G.cov)
colnames(mat)<- traitnames2
rownames(mat)<- traitnames2
diag(mat) <- NA
myplotcor(mat,diag=T, type="lower",col="black",new=1,sign.mat=sign.mat)
sign.mat <- matrix(0,ncol=Gn,nrow=Gn)
par(col = "black",lty=1,lwd=2)
diag(mat) <- NA
myplotcor(mat,diag=T, type="upper",col="black",new=0,sign.mat=sign.mat)
#legend("topleft", c("Males", "Female"), col=c("black", "black"), lty =c(6,1), bty="n")
dev.off()
##################################################################################################################



###############################
# Symmetry of B matrix
###############################
# Following Gosden and Chenoweth (2014), and Sztepanacz and Houle (2019)

# Quick reference:
# From: https://urldefense.proofpoint.com/v2/url?u=https-3A__www.sciencedirect.com_topics_mathematics_skew-2Dsymmetric-2Dmatrix&d=DwIGAw&c=sJ6xIWYx-zLMB3EPkvcnVg&r=XBGf66T2ZPyLIsd9UZy7Vw&m=nlWDcqbye4Bd8ZGNFYE99PG6SAcbDQfrES6_NsbBJFY&s=k-oPT7RnHKFDDqFhM-805b1NsfRXxYJVcrSXaHh72XU&e= 
#
# Theorem 1.14
# 
# If A is a square matrix, then
# (1)
# A + A^T is symmetric, and
# 
# (2)
# A − A^T is skew-symmetric.

# The matrix we want is B, which we'll relable here as A:

# median estimate of G.cov 
G.median <- matrix(apply(model$VCV[,1:Glength],2,median),nrow=Gn,ncol=Gn,byrow=F)
G.cov <- matrix(colMeans(model$VCV[,1:Glength]),nrow=Gn,ncol=Gn,byrow=F)

# B
A <- G.cov[1:n, (n+1):Gn]

# B^t 
# could use this or could use t(A), they should be the same.
At <- G.cov[(n+1):Gn, 1:n]

# Symetric component of B
S <- 0.5*(A+t(A))

# Skew-symetric component of B
N <- 0.5*(A-t(A))

# # These should add up to == A, thus,  let's check
A - (S+N) # should == 0; yes
# also
S - (t(S)) # should == 0; yes
# and
N - (-t(N)) # should == 0; yes

# Mean sum of squares for S, N and A
SS.S <- mean((Matrix::rowSums(S ** 2) + Matrix::colSums(S ** 2)) / 2)
SS.N <- mean((Matrix::rowSums(N ** 2) + Matrix::colSums(N ** 2)) / 2)
SS.A <- mean((Matrix::rowSums(A ** 2) + Matrix::colSums(A ** 2)) / 2)

# proportion of A that is symetric?
SS.S / SS.A 

# proportion of A that is skew symetric (or assymmetric)?
SS.N / SS.A 


# Get HPD intervals for that point estimate 

# Extract 1000 stored G matrices from model (same as above; ok to overwrite). 
f <- function(x){list(matrix(x,nrow=Gn,ncol=Gn,byrow=F))}
Gmcmc <- lapply(apply(model$VCV[,1:Glength], 1, f), "[[", 1) # Full G matrices

# 1000 B matrices
A <- list()
for (i in 1:1000) {
  A[[i]] <- Gmcmc[[i]][1:n, (n+1):Gn]
}

# 1000 symetric components
S <- list()
for (i in 1:1000) {
  S[[i]] <- 0.5*(A[[i]]+t(A[[i]]))
}

# 1000 skew symetric components
N <- list()
for (i in 1:1000) {
  N[[i]] <- 0.5*(A[[i]]-t(A[[i]]))
}

# set up empty vectors
res.SS.S <- 1:1000
res.SS.N <- 1:1000
res.SS.A <- 1:1000

# for the mean sums of squares 
for (i in 1:1000) {
  res.SS.S[i] <- mean((Matrix::rowSums(S[[i]] ** 2) + Matrix::colSums(S[[i]] ** 2)) / 2)
  res.SS.N[i] <- mean((Matrix::rowSums(N[[i]] ** 2) + Matrix::colSums(N[[i]] ** 2)) / 2)
  res.SS.A[i] <- mean((Matrix::rowSums(A[[i]] ** 2) + Matrix::colSums(A[[i]] ** 2)) / 2)
  
}

# Proportion of A that is N (or S if you wanted to look at that; they're complimentary)
prop.S <- 1:1000
prop.N <- 1:1000
for (i in 1:1000) {
  prop.S[i] <- res.SS.S[i] / res.SS.A[i] 
  prop.N[i] <- res.SS.N[i] / res.SS.A[i] 
  
}


# A more accurate point estimate for the proportion of B that is skew symetric is given by the resampled posterior mode...


# skew sym. or asym.
posterior.mode(as.mcmc(prop.N), na.rm=T) # (note any NAs, so that denominator can be adjusted below)
HPDinterval(as.mcmc(prop.N)) # 95% HPD intervals
# symmetric
posterior.mode(as.mcmc(prop.S), na.rm=T) # (note any NAs, so that denominator can be adjusted below)
HPDinterval(as.mcmc(prop.S)) # 95% HPD intervals


####################################################
## Antagonistic v. concordant subspace analysis ###
####################################################
# Following Sztepanacz and Houle (2019)

# n-dimentional identity matrix
I <- diag(n)

# n eigen vectors that span the space of I 
e.I <- eigen(I)$vectors

# Make that anything but "1's" (Sztepanacz and Houle, 2019) and call it Em 
Em <- e.I/sqrt(2)

# Sm
Em.neg <- Em*(-1) # negative
Sm.top <- cbind(Em, Em) # top part
Sm.bot <- cbind(Em, Em.neg) # bottom part
Sm <- rbind(Sm.top, Sm.bot) # compiled

# Project G onto this space
Gca <- t(Sm) %*% G.cov %*% (Sm)

# index the antagonistic and concordant subspaces
Gc <- Gca[1:n,1:n]
Ga <- Gca[(n+1):Gn,(n+1):Gn]

# Get "trace" of each (trace = the sum of the eigen values)
con.trace <- sum(eigen(Gc)$values)
ant.trace <- sum(eigen(Ga)$values)
tot.trace.Gca <- sum(eigen(Gca)$values)
tot.trace.G.cov <- sum(eigen(G.cov)$values)

# Concordant and antagonistic genetic variances
tot.trace.Gca # this should ==
tot.trace.G.cov # ...this.
con.trace # concordant variance
ant.trace # antagonistic variance

# Percent variance
con.trace/tot.trace.G.cov
ant.trace/tot.trace.G.cov

# Genetic variances for each of the eigen vectors of G.cov, and their concordant and antagonistic genetic variances
G.cov_eigenvalues <- eigen(G.cov)$values
G.c_eigenvalues <- eigen(Gc)$values
G.a_eigenvalues <- eigen(Ga)$values
G.cov_eigenvalues 
G.c_eigenvalues 
G.a_eigenvalues 

# Get HPD mode and HPD intervals for those estimates

# set up empty vectors and df
prop.con <- 1:1000
prop.ant <- 1:1000
ev1 <- 1:1000
ev2 <- 1:1000
ev3 <- 1:1000
ev4 <- 1:1000
ev5 <- 1:1000
ev6 <- 1:1000
# ev7 <- 1:1000
# ev8 <- 1:1000

c.a.Vg <- cbind(ev1, ev2, ev3, ev4, ev5, ev6)#, ev7, ev8) #These change based on the number of traits
tot.Vg <- cbind(ev1, ev2, ev3, ev4, ev5, ev6)#, ev7, ev8) #These change based on the number of traits

# repeat setps above but looping through the 1000 stored estimates 
for (i in 1:MCMCsamp) {
  Gca.temp <- t(Sm) %*% Gmcmc[[i]] %*% (Sm) # Project G onto Sm
  Gc.temp <- Gca.temp[1:n,1:n]# index the antagonistic...
  Ga.temp <- Gca.temp[(n+1):Gn,(n+1):Gn] # ...and concordant subspaces
  con <- sum(eigen(Gc.temp)$values)
  ant <- sum(eigen(Ga.temp)$values)
  tot <- sum(eigen(Gmcmc[[i]])$values)
  prop.con[i] <- con/tot # proportion concordant
  prop.ant[i] <- ant/tot # proportion antagonistic
  c.a.Vg[i,] <- c(eigen(Gc.temp)$values, eigen(Ga.temp)$values)
  tot.Vg[i,] <- eigen(Gmcmc[[i]])$values
}


# overal proportion of Gmf that sexually concordant
posterior.mode(as.mcmc(prop.con), na.rm=T) # 
HPDinterval(as.mcmc(prop.con)) # 95% HPD intervals

# overal proportion of Gmf that sexually antagonistic
posterior.mode(as.mcmc(prop.ant), na.rm=T) # 
HPDinterval(as.mcmc(prop.ant)) # 95% HPD intervals

# set up the figure file
Vg.con.ant.vec <- 1:Gn
con.ant.CI <- 1:Gn
Vg.tot.vec <- 1:Gn
tot.CI <- 1:Gn
plot.a.c.tot <- cbind(Vg.con.ant.vec, con.ant.CI, Vg.tot.vec, tot.CI)


for (i in 1:Gn) {
  plot.a.c.tot[i,1] <- posterior.mode(as.mcmc(c.a.Vg[,i]), na.rm=T) # 
  plot.a.c.tot[i,2] <- paste(HPDinterval(as.mcmc(c.a.Vg[,i]))[1], HPDinterval(as.mcmc(c.a.Vg[,i]))[2])
  plot.a.c.tot[i,3] <- posterior.mode(as.mcmc(tot.Vg[,i]), na.rm=T) # 
  plot.a.c.tot[i,4] <- paste(HPDinterval(as.mcmc(tot.Vg[,i]))[1], HPDinterval(as.mcmc(tot.Vg[,i]))[2])
}

### Plotting antagonistic vs concordant selection.
### Will need to reformat data use the melt function or by hand in excel.
## Below, reformated (long-format) data reread in.

plot_ant_con_LH <- read.csv("~/PATH/plot.a.c.tot_scaled.ThreeLHtraits.parexp_adjusted.csv")
plot_ant_con_LH$Vector <- as.factor(plot_ant_con_LH$Vector)
plot_ant_con_LH$Vg <- as.factor(plot_ant_con_LH$Vg)

plot_ant_con_VOC <- read.csv("~/PATH/plot.a.c.tot_Final.VOC.4.repro_adjusted.csv")
plot_ant_con_VOC$Vector <- as.factor(plot_ant_con_VOC$Vector)
plot_ant_con_VOC$Vg <- as.factor(plot_ant_con_VOC$Vg)


LH_plot <- ggplot(plot_ant_con_LH, aes(xmin = as.numeric(Vector) - .4,xmax = as.numeric(Vector) + .4, 
                                       x= Vector, ymin= 10E-5, ymax=Vg.vect)) +
  geom_rect(position=position_dodge(.7), aes(fill = Vg)) +
  scale_y_log10()+
  geom_errorbar(aes(ymin=  lower.CI, ymax=  upper.CI, group = Vg), 
                colour="black",  size=.2, width = .2,position = position_dodge(width = 0.70))+
  ylab("log10(Eigenvalue)")+
  xlab("Vector") +
  scale_fill_manual(values=wes_palette(n=3, name="Darjeeling1"),
                    name= "Genetic variance",
                    labels = c("GCA Antagonistic", "GCA Concordant", "Total GMF"))+
  theme_bw()

# Figure 2
VOC_plot <- ggplot(plot_ant_con_VOC, aes(xmin = as.numeric(Vector) - .4,xmax = as.numeric(Vector) + .4, 
                                         x= Vector, ymin= 10E-5, ymax=Vg.vect)) +
  geom_rect(position=position_dodge(.7), aes(fill = Vg)) +
  scale_y_log10()+
  geom_errorbar(aes(ymin=  lower.CI, ymax=  upper.CI, group = Vg), 
                colour="black",  size=.2, width = .2,position = position_dodge(width = 0.70))+
  ylab("log10(Eigenvalue)")+
  xlab("Vector") +
  scale_fill_manual(values=wes_palette(n=3, name="Darjeeling1"),
                    name= "Genetic variance",
                    labels = c("GCA Antagonistic", "GCA Concordant", "Total GMF"))+
  theme_bw()


png("GCA_plot.png",type="cairo", units="in", 
    width=7, height=5, pointsize=8, res=300)
ggarrange(LH_plot, VOC_plot, ncol = 2, nrow = 1, common.legend = TRUE, legend = "bottom", labels= c("A", "B") )
dev.off()

####################### Function for elipse plot #################################
## Puentes et al 2016

#START
myplotcor <- function (corr, outline = TRUE, col = "grey", numbers = FALSE, 
                       type = c("full", "lower", "upper"), diag = (type == "full"), 
                       bty = "n", axes = FALSE, xlab = "", ylab = "", asp = 1, cex.lab = par("cex.lab"), 
                       cex = 0.75 * par("cex"), mar = 0.1 + c(2, 2, 4, 2), new, sign.mat, ...) 
{
  savepar <- par(pty = "s", mar = mar)
  on.exit(par(savepar))
  if (is.null(corr)) 
    return(invisible())
  if ((!is.matrix(corr)) || (round(min(corr, na.rm = TRUE), 
                                   6) < -1) || (round(max(corr, na.rm = TRUE), 6) > 1)) 
    stop("Need a correlation matrix")
  if(new=="1") {plot.new()}
  else {par(new = TRUE)}
  rowdim <- dim(corr)[1]
  coldim <- dim(corr)[2]
  rowlabs <- dimnames(corr)[[1]]
  collabs <- dimnames(corr)[[2]]
  if (is.null(rowlabs)) 
    rowlabs <- 1:rowdim
  if (is.null(collabs)) 
    collabs <- 1:coldim
  rowlabs <- as.character(rowlabs)
  collabs <- as.character(collabs)
  col <- rep(col, length = length(corr))
  dim(col) <- dim(corr)
  type <- match.arg(type)
  cols <- 1:coldim
  rows <- 1:rowdim
  xshift <- 0
  yshift <- 0
  if (!diag) {
    if (type == "upper") {
      cols <- 2:coldim
      rows <- 1:(rowdim - 1)
      xshift <- 1
    }
    else if (type == "lower") {
      cols <- 1:(coldim - 1)
      rows <- 2:rowdim
      yshift <- -1
    }
  }
  maxdim <- max(length(rows), length(cols))
  plt <- par("plt")
  xlabwidth <- max(strwidth(rowlabs[rows], units = "figure", 
                            cex = cex.lab))/(plt[2] - plt[1])
  xlabwidth <- xlabwidth * maxdim/(1 - xlabwidth)
  ylabwidth <- max(strwidth(collabs[cols], units = "figure", 
                            cex = cex.lab))/(plt[4] - plt[3])
  ylabwidth <- ylabwidth * maxdim/(1 - ylabwidth)
  plot(c(-xlabwidth - 0.5, maxdim + 0.5), c(0.5, maxdim + 1 + 
                                              ylabwidth), type = "n", bty = bty, axes = axes, xlab = "", 
       ylab = "", asp = asp, cex.lab = cex.lab, ...)
  text(rep(0, length(rows)), length(rows):1, labels = rowlabs[rows], 
       adj = 1, cex = cex.lab)
  text(cols - xshift, rep(length(rows) + 1, length(cols)), 
       labels = collabs[cols], srt = 90, adj = 0, cex = cex.lab)
  mtext(xlab, 1, 0)
  mtext(ylab, 2, 0)
  mat <- diag(c(1, 1))
  plotcorrInternal <- function() {
    if (i == j && !diag) 
      return()
    if (!numbers) {
      mat[1, 2] <- corr[i, j]
      mat[2, 1] <- mat[1, 2]
      ell <- ellipse(mat, t = 0.43)
      ell[, 1] <- ell[, 1] + j - xshift
      ell[, 2] <- ell[, 2] + length(rows) + 1 - i - yshift
      #polygon(ell, col = col[i, j])	
      if (outline)
        if(sign.mat[i, j]==1) {lines(ell, lwd=2)}
      else {lines(ell)}
    }
    else {
      text(j + 0.3 - xshift, length(rows) + 1 - i - yshift, 
           round(10 * corr[i, j], 0), adj = 1, cex = cex)
    }
  }
  for (i in 1:dim(corr)[1]) {
    for (j in 1:dim(corr)[2]) {
      if (type == "full") {
        plotcorrInternal()
      }
      else if (type == "lower" && (i >= j)) {
        plotcorrInternal()
      }
      else if (type == "upper" && (i <= j)) {
        plotcorrInternal()
      }
    }
  }
  invisible()
}
# END