library(MCMCglmm) 
library(tidyr) 
library(plyr) 
library(dplyr)
library(MuMIn)
library(ggplot2)
library(MasterBayes)
library(car)
library(matrixcalc)
library(MASS)

###########
# Data
###########

#This file is what I'm using...
#data <- read.csv("~/PATH/LH.traits.data.NONA.csv") #Growth and developmental traits
data <- read.csv("~/PATH/data.homemade.traits.csv") #Morphology and physiology traits

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

####################
## Define terms  ##
####################

n=1 # number of traits
m=2 # number of sexes
Blength=n^2 # size of Gm, Gf, or B
Gn=n*m # total length of diag of Gmf
Glength=Gn^2 # size of Gmf

# Error structure 
gaus <- rep("gaussian", n) # n traits, all with Gaussian error structure

# Setting priors
full.s.pu <- list(R=list(V=diag(0.01, Gn), nu=1.001),
                      G = list(G1 = list(V=diag(0.001, Gn), nu=1.001),
                               G2 = list(V = 0.001, nu = 1.001),
                               G3 = list(V = 0.001, nu = 1.001)))

full.u.pu <- list(R=list(V=0.01, nu=1.001),
                      G = list(G1 = list(V=0.001, nu=1.001),
                               G2 = list(V = 0.001, nu = 1.001),
                               G3 = list(V = 0.001, nu = 1.001)))

red.pu <- list(R=list(V=0.01, nu=1.001),
                      G = list(G1 = list(V=0.001, nu=1.001),
                            G2 = list(V = 0.001, nu = 1.001)))

red.pu.2 <- list(R=list(V=0.01, nu=1.001),
                      G = list(G1 = list(V=0.001, nu=1.001),
                              G2 = list(V = 0.001, nu = 1.001),
                              G3 = list(V = 0.001, nu = 1.001)))


#############
## Models  ##
#############

# Full sexed univariate model
full.model.rel_repro <- MCMCglmm((rel_repro) ~ ssex, 
                                 random = ~us(ssex):famid + 
                                   sampid+
                                   Date.PTR, 
                                 rcov = ~idh(ssex):units, 
                                 data = data,
                                 family = gaus,
                                 prior = full.s.pu.VOC, 
                                 nitt=1005000,
                                 burnin=5000,
                                 thin=1000,
                                 slice=TRUE,
                                 verbose = T,
                                 pr=T)

par(mar=c(1,1,1,1))
plot(full.model.rel_repro)
full.model.rel_repro$DIC

# Unsexed full univariate model
full.model.unsexed.rel_repro <- MCMCglmm((rel_repro) ~ ssex, 
                                         random = ~us(1):famid + 
                                         sampid +
                                         Date.PTR,
                                         #rcov = ~idh(1):units, # no residual covariance between sexes because those measure stem from different individuals, see Haddfield's course notes, top of pg. 71
                                         data = data,
                                         family = gaus,
                                         prior = full.u.pu.VOC, 
                                         nitt=1005000,
                                         burnin=5000,
                                         thin=1000,
                                         slice=TRUE,
                                         verbose = T,
                                         pr=T)


full.model.unsexed.rel_repro$DIC

## Reduced univariate model dropping just family
reduced.sex.rel_repro <- MCMCglmm((rel_repro) ~ ssex, 
                                  random = ~sampid + 
                                  Date.PTR, 
                                  #rcov = ~idh(1):units, # no residual covariance between sexes because those measure stem from different individuals, see Haddfield's course notes, top of pg. 71
                                  data = data,
                                  family = gaus,
                                  prior = red.pu.VOC, 
                                  nitt=1005000,
                                  burnin=5000,
                                  thin=1000,
                                  slice=TRUE,
                                  verbose = T,
                                  pr=T)

reduced.sex.rel_repro$DIC

############################
## Full model without sex ##
############################

# Motivation for doing this is that mean diff between sex could be a major feature in data that
# makes it hard to fit other suttle differences.
reduced.unsex.rel_repro <- MCMCglmm((rel_repro) ~ 1, 
                                    random = ~famid +
                                    sampid + 
                                    Date.PTR, 
                                    rcov = ~idh(1):units, 
                                    data = data,
                                    family = gaus,
                                    prior = red.pu.2.VOC, 
                                    nitt=1005000,
                                    burnin=5000,
                                    thin=1000,
                                    slice=TRUE,
                                    verbose = T,
                                    pr=T)

reduced.unsex.rel_repro$DIC 

reduced.unsex.nofam.rel_repro <- MCMCglmm((rel_repro) ~ 1, 
                                          random = ~sampid + 
                                            Date.PTR, 
                                          rcov = ~idh(1):units, 
                                          data = data,
                                          family = gaus,
                                          prior = red.pu.VOC, 
                                          nitt=1005000,
                                          burnin=5000,
                                          thin=1000,
                                          slice=TRUE,
                                          verbose = T,
                                          pr=T)

reduced.unsex.nofam.rel_repro$DIC

