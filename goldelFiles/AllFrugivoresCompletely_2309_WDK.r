# Bastian Göldel
# Project: Global Palm and Mammals Functional traits distributions and relation  (category fruits = all frugivores)
# added: extinct frugivores & replaced current ranges of present frugivores by natural ranges


# libraries
#library("lavaan", lib="P:/Documents/Uni Aarhus/PhD/data from Daniel/Data for Bastian/NewWorld")

#Workspace
rm(list=ls())
#setwd("C:/Analysis_Bastian/Mammals and Palms/Office")
setwd("D:/Data/wkissli1/My Documents/Projects/Global/06_Palms/02_Bastian/Analysis/Global_Frugivores")

library(lavaan)


######## CATEGORY FRUIT = ALL (also extinct and natural) Frugivores #########


#Load data
Frugivores<-read.table("ALL_MAMMALS_COMPLETELY_230915.txt", header=T,sep="\t")
str(Frugivores)
names(Frugivores)
Frugivores <- na.omit(Frugivores) # ignores NAs
is.na(Frugivores)                      # check
str(Frugivores)
names(Frugivores)
summary(Frugivores)



# PCA of current climate variables
all_climate <- Frugivores[,12:17]
str(all_climate)
all_climate
PCA_clim<-prcomp(all_climate, scale=TRUE)
summary(PCA_clim)    # first 3 axes explain ~93.8%
biplot(PCA_clim, cex=0.9, xlim=c(-0.23,0.23))

PCA_clim$x[,1:3]
Frugivores<-cbind(Frugivores, PCA_clim$x[,1:3])
is.na(Frugivores)
str(Frugivores)
summary(Frugivores)



# Histograms
hist(log(Frugivores$FruitSize))         # palm fruit size
hist(log(Frugivores$BodySize))          # mammal body size
hist(Frugivores$NPP_mean)               # NPP (net prim. production)
hist(Frugivores$ensLGM_Tmean)           # mean LGM ANOM Temp
hist(log(Frugivores$ensLGM_Pmean))      # mean LGM ANOM Prec
hist(Frugivores$Soil)                   # % of sandy soil (depth 45cm)
hist(Frugivores$PC1)
hist(Frugivores$PC2)
hist(Frugivores$PC3)


# plot single relationships to check for non-linear response
plot(log(Frugivores$FruitSize) ~ log(Frugivores$BodySize))
plot(log(Frugivores$FruitSize) ~ Frugivores$PC1)
plot(log(Frugivores$FruitSize) ~ Frugivores$PC2)
plot(log(Frugivores$FruitSize) ~ Frugivores$PC3)
plot(log(Frugivores$FruitSize) ~ Frugivores$NPP_mean)
plot(log(Frugivores$FruitSize) ~ Frugivores$ensLGM_Tmean)
plot(log(Frugivores$FruitSize) ~ log(Frugivores$ensLGM_Pmean))
plot(log(Frugivores$FruitSize) ~ Frugivores$Soil)

plot(log(Frugivores$BodySize) ~ Frugivores$PC1)
plot(log(Frugivores$BodySize) ~ Frugivores$PC2)
plot(log(Frugivores$BodySize) ~ Frugivores$PC3)
plot(log(Frugivores$BodySize) ~ Frugivores$NPP_mean)
plot(log(Frugivores$BodySize) ~ Frugivores$ensLGM_Tmean)
plot(log(Frugivores$BodySize) ~ log(Frugivores$ensLGM_Pmean))
plot(log(Frugivores$BodySize) ~ Frugivores$Soil)



#------------------ SEM analyses --------------------------------
# Criteria for fit measures:
# Model Chi-Square with its df and p-value: prefer p-value greater than 0.05
# Root Mean Square Error of Approximation (RMSEA): prefer lower 90%CI to be < 0.05
# Comparative Fit Index (CFI): prefer value greater than 0.90


#The dataset
colnames(Frugivores)
str(Frugivores)
summary(Frugivores)

## Recode vars to roughly same scale (recommended!)
FruitSize <- log(Frugivores$FruitSize)
BodySize <- log(Frugivores$BodySize)
NPP_mean <- Frugivores$NPP_mean/10000
ensLGM_Tmean <- Frugivores$ensLGM_Tmean/100
ensLGM_Pmean <- log(Frugivores$ensLGM_Pmean)
Soil <- Frugivores$Soil/10
PC1 <- Frugivores$PC1
PC2 <- Frugivores$PC2
PC3 <- Frugivores$PC3

Frugivores.dat <- data.frame(FruitSize, BodySize, NPP_mean, ensLGM_Tmean, ensLGM_Pmean, Soil, PC1, PC2, PC3)
summary(Frugivores.dat)
varTable(Frugivores.dat)



##### ALL biogeographical regions together ######
# SEM approach
Frugivores_Model_1 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + ensLGM_Pmean + PC1 + PC2 + PC3
BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean
NPP_mean ~ Soil + PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean
'
sem.fit.1 <- sem(Frugivores_Model_1, data=Frugivores.dat)
summary(sem.fit.1, stand=T, rsq=T, fit.measures=T)    

modindices(sem.fit.1)
resid(sem.fit.1, type="standardized")


# delete PC2 from the fruit size model
Frugivores_Model_2 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + ensLGM_Pmean + PC1 + PC3
BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean
NPP_mean ~ Soil + PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean
'
sem.fit.2 <- sem(Frugivores_Model_2, data=Frugivores.dat)
summary(sem.fit.2, stand=T, rsq=T, fit.measures=T)  


# delete ensLGM_Tmean from the body size model
Frugivores_Model_3 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + ensLGM_Pmean + PC1 + PC3
BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Pmean
NPP_mean ~ Soil + PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean
'
sem.fit.3 <- sem(Frugivores_Model_3, data=Frugivores.dat)
summary(sem.fit.3, stand=T, rsq=T, fit.measures=T)  


# delete Soil from the NPP model
Frugivores_Model_4 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + ensLGM_Pmean + PC1 + PC3
BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Pmean
NPP_mean ~ PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean
'
sem.fit.4 <- sem(Frugivores_Model_4, data=Frugivores.dat)
summary(sem.fit.4, stand=T, rsq=T, fit.measures=T)


# delete NPP from the body size model
Frugivores_Model_5 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + ensLGM_Pmean + PC1 + PC3
BodySize ~ PC1 + PC2 + PC3 + ensLGM_Pmean
NPP_mean ~ PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean
'
sem.fit.5 <- sem(Frugivores_Model_5, data=Frugivores.dat)
summary(sem.fit.5, stand=T, rsq=T, fit.measures=T)


# delete PC1 from the body size model
Frugivores_Model_6 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + ensLGM_Pmean + PC1 + PC3
BodySize ~ PC2 + PC3 + ensLGM_Pmean
NPP_mean ~ PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean
'
sem.fit.6 <- sem(Frugivores_Model_6, data=Frugivores.dat)
summary(sem.fit.6, stand=T, rsq=T, fit.measures=T)


# delete ensLGM_Pmean from the body size model
Frugivores_Model_7 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + ensLGM_Pmean + PC1 + PC3
BodySize ~ PC2 + PC3
NPP_mean ~ PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean
'
sem.fit.7 <- sem(Frugivores_Model_7, data=Frugivores.dat)
summary(sem.fit.7, stand=T, rsq=T, fit.measures=T)


# delete PC3 from the fruit size model
Frugivores_Model_8 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + ensLGM_Pmean + PC1
BodySize ~ PC2 + PC3
NPP_mean ~ PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean
'
sem.fit.8 <- sem(Frugivores_Model_8, data=Frugivores.dat)
summary(sem.fit.8, stand=T, rsq=T, fit.measures=T)


# delete PC1 from the fruit size model
Frugivores_Model_9 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + ensLGM_Pmean
BodySize ~ PC2 + PC3
NPP_mean ~ PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean
'
sem.fit.9 <- sem(Frugivores_Model_9, data=Frugivores.dat)
summary(sem.fit.9, stand=T, rsq=T, fit.measures=T)


# delete ensLGM_Pmean from the fruit size model
Frugivores_Model_10 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean
BodySize ~ PC2 + PC3
NPP_mean ~ PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean
'
sem.fit.10 <- sem(Frugivores_Model_10, data=Frugivores.dat)
summary(sem.fit.10, stand=T, rsq=T, fit.measures=T)


# delete PC1 from the NPP model
Frugivores_Model_11 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean
BodySize ~ PC2 + PC3
NPP_mean ~ PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean
'
sem.fit.11 <- sem(Frugivores_Model_11, data=Frugivores.dat)
summary(sem.fit.11, stand=T, rsq=T, fit.measures=T)

#RESULTS (by Daniel): 
# present and past climate affects NPP
# mean frugivore body sizes is driven by climate
# palm mean fruit size is driven by 
#   (1) LGM TEMP fluctuations (larger fruit sizes with higher TEMP anomalies)
#   (2) NPP (negative effect, i.e. mean fruit size decreases with increasing NPP)
#   (3) Frugivore body size (positive effect, i.e. larger fruit sizes with larger mean animal body size) 
#   (4) Soil (positive effect, i.e. larger fruits with more soil types???)


##### NEW WORLD (by Daniel) ######
str(Frugivores)
Frugivores_NW <- Frugivores [ which(Frugivores$THREEREALM =='NewWorld'), ]
str(Frugivores_NW)

# PCA of current climate variables
all_climate <- Frugivores_NW[,12:17]
str(all_climate)
all_climate
PCA_clim<-prcomp(all_climate, scale=TRUE)
summary(PCA_clim)    # first 3 axes explain ~92%
biplot(PCA_clim, cex=0.9, xlim=c(-0.23,0.23))

PCA_clim$x[,1:3]                           
Frugivores_NW<-cbind(Frugivores_NW, PCA_clim$x[,1:3])
is.na(Frugivores_NW)
str(Frugivores_NW)
summary(Frugivores_NW)

hist(log(Frugivores_NW$FruitSize))         # palm fruit size
hist(log(Frugivores_NW$BodySize))          # mammal body size
# Histograms (predictor variables)
hist(Frugivores_NW$PC1)                    # mean annual precipitation
hist(Frugivores_NW$PC2)                    # precipitation seasonlity
hist(Frugivores_NW$PC3)                    # mean annual temperature
hist(Frugivores_NW$NPP_mean)               # NPP (net prim. production)
hist(Frugivores_NW$ensLGM_Tmean)           # mean LGM ANOM Temp
hist(log(Frugivores_NW$ensLGM_Pmean))      # mean LGM ANOM Prec
hist(Frugivores_NW$Soil)

# adapt data
FruitSize <- log(Frugivores_NW$FruitSize)
BodySize <- log(Frugivores_NW$BodySize)
NPP_mean <- Frugivores_NW$NPP_mean/10000
ensLGM_Tmean <- Frugivores_NW$ensLGM_Tmean/100
ensLGM_Pmean <- log(Frugivores_NW$ensLGM_Pmean)
Soil <- Frugivores_NW$Soil/10
PC1 <- Frugivores_NW$PC1
PC2 <- Frugivores_NW$PC2
PC3 <- Frugivores_NW$PC3

Frugivores_NW.dat <- data.frame(FruitSize, BodySize, NPP_mean, ensLGM_Tmean, ensLGM_Pmean, Soil, PC1, PC2, PC3)
summary(Frugivores_NW.dat)
varTable(Frugivores_NW.dat)

# SEM approach
Frugivores_NW_Model_1 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + ensLGM_Pmean + PC1 + PC2 + PC3
BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean
NPP_mean ~ Soil + PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean
'
sem.fit.1 <- sem(Frugivores_NW_Model_1, data=Frugivores_NW.dat)
summary(sem.fit.1, stand=T, rsq=T, fit.measures=T)

modindices(sem.fit.1)
resid(sem.fit.1, type="standardized")

# delete PC1 from the FruitSize model
Frugivores_NW_Model_2 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + ensLGM_Pmean + PC2 + PC3
BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean
NPP_mean ~ Soil + PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean
'
sem.fit.2 <- sem(Frugivores_NW_Model_2, data=Frugivores_NW.dat)
summary(sem.fit.2, stand=T, rsq=T, fit.measures=T)

# delete ensLGM_Tmean from the FruitSize model
Frugivores_NW_Model_3 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Pmean + PC2 + PC3
BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean
NPP_mean ~ Soil + PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean
'
sem.fit.3 <- sem(Frugivores_NW_Model_3, data=Frugivores_NW.dat)
summary(sem.fit.3, stand=T, rsq=T, fit.measures=T)

# delete ensLGM_Tmean from the BodySize model
Frugivores_NW_Model_4 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Pmean + PC2 + PC3
BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Pmean
NPP_mean ~ Soil + PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean
'
sem.fit.4 <- sem(Frugivores_NW_Model_4, data=Frugivores_NW.dat)
summary(sem.fit.4, stand=T, rsq=T, fit.measures=T)


# delete PC3 from the FruitSize model
Frugivores_NW_Model_5 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Pmean + PC2
BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Pmean
NPP_mean ~ Soil + PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean
'
sem.fit.5 <- sem(Frugivores_NW_Model_5, data=Frugivores_NW.dat)
summary(sem.fit.5, stand=T, rsq=T, fit.measures=T)

# delete PC3 from the NPP_mean model
Frugivores_NW_Model_6 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Pmean + PC2
BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Pmean
NPP_mean ~ Soil + PC1 + PC2 + ensLGM_Pmean + ensLGM_Tmean
'
sem.fit.6 <- sem(Frugivores_NW_Model_6, data=Frugivores_NW.dat)
summary(sem.fit.6, stand=T, rsq=T, fit.measures=T)


# delete ensLGM_Pmean from the NPP_mean model
Frugivores_NW_Model_7 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Pmean + PC2
BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Pmean
NPP_mean ~ Soil + PC1 + PC2 + ensLGM_Tmean
'
sem.fit.7 <- sem(Frugivores_NW_Model_7, data=Frugivores_NW.dat)
summary(sem.fit.7, stand=T, rsq=T, fit.measures=T)

# delete PC1 from the NPP_mean model
Frugivores_NW_Model_8 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Pmean + PC2
BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Pmean
NPP_mean ~ Soil + PC2 + ensLGM_Tmean
'
sem.fit.8 <- sem(Frugivores_NW_Model_8, data=Frugivores_NW.dat)
summary(sem.fit.8, stand=T, rsq=T, fit.measures=T)

# delete ensLGM_Pmean from the FruitSize model
Frugivores_NW_Model_9 <- 'FruitSize ~ BodySize + Soil + NPP_mean + PC2
BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Pmean
NPP_mean ~ Soil + PC2 + ensLGM_Tmean
'
sem.fit.9 <- sem(Frugivores_NW_Model_9, data=Frugivores_NW.dat)
summary(sem.fit.9, stand=T, rsq=T, fit.measures=T)


# delete Soil from the FruitSize model
Frugivores_NW_Model_10 <- 'FruitSize ~ BodySize + NPP_mean + PC2
BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Pmean
NPP_mean ~ Soil + PC2 + ensLGM_Tmean
'
sem.fit.10 <- sem(Frugivores_NW_Model_10, data=Frugivores_NW.dat)
summary(sem.fit.10, stand=T, rsq=T, fit.measures=T)


#RESULTS: 
# present and past climate + soil affects NPP
# mean frugivore body sizes is driven by climate, NPP, and LGM PREC anomaly
# palm mean fruit size is mainly driven by 
#   -Frugivore body size (positive effect, i.e. larger fruit sizes with larger mean animal body size) 
# but also by
#  -NPP
#  -PC2 





################################################################################
#PREVIOUS SCRIPT FROM BASTIAN
##### Neotropics ######

Frugivores_Neo <- Frugivores [ which(Frugivores$REALM_LONG =='Neotropics'), ]
Frugivores_Neo

# PCA of current climate variables
all_climate <- Frugivores_Neo[,12:17]
str(all_climate)
all_climate
PCA_clim<-prcomp(all_climate, scale=TRUE)
summary(PCA_clim)    # first 3 axes explain ~92%
biplot(PCA_clim, cex=0.9, xlim=c(-0.23,0.23))

PCA_clim$x[,1:3]                           
Frugivores_Neo<-cbind(Frugivores_Neo, PCA_clim$x[,1:3])
is.na(Frugivores_Neo)
str(Frugivores_Neo)
summary(Frugivores_Neo)


hist(log(Frugivores_Neo$FruitSize))         # palm fruit size
hist(log(Frugivores_Neo$BodySize))          # mammal body size
# Histograms (predictor variables)
hist(Frugivores_Neo$PC1)                    # mean annual precipitation
hist(Frugivores_Neo$PC2)                    # precipitation seasonlity
hist(Frugivores_Neo$PC3)                    # mean annual temperature
hist(Frugivores_Neo$NPP_mean)               # NPP (net prim. production)
hist(Frugivores_Neo$ensLGM_Tmean)           # mean LGM ANOM Temp
hist(log(Frugivores_Neo$ensLGM_Pmean))      # mean LGM ANOM Prec
hist(Frugivores_Neo$Soil)


colnames(Frugivores_Neo)
str(Frugivores_Neo)
summary(Frugivores_Neo)

# adapt data
FruitSize <- log(Frugivores_Neo$FruitSize)
BodySize <- log(Frugivores_Neo$BodySize)
NPP_mean <- Frugivores_Neo$NPP_mean/10000
ensLGM_Tmean <- Frugivores_Neo$ensLGM_Tmean/100
ensLGM_Pmean <- log(Frugivores_Neo$ensLGM_Pmean)
Soil <- Frugivores_Neo$Soil/10
PC1 <- Frugivores_Neo$PC1
PC2 <- Frugivores_Neo$PC2
PC3 <- Frugivores_Neo$PC3

Frugivores_Neo.dat <- data.frame(FruitSize, BodySize, NPP_mean, ensLGM_Tmean, ensLGM_Pmean, Soil, PC1, PC2, PC3)
summary(Frugivores_Neo.dat)
varTable(Frugivores_Neo.dat)


# SEM approach
Frugivores_Neo_Model_1 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + ensLGM_Pmean + PC1 + PC2 + PC3
BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean
NPP_mean ~ Soil + PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean
'
sem.fit.1 <- sem(Frugivores_Neo_Model_1, data=Frugivores_Neo.dat)
summary(sem.fit.1, stand=T, rsq=T, fit.measures=T)

modindices(sem.fit.1)
resid(sem.fit.1, type="standardized")


# delete PC2 from the NPP model
Frugivores_Neo_Model_2 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + ensLGM_Pmean + PC1 + PC2 + PC3
BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean
NPP_mean ~ Soil + PC1 + PC3 + ensLGM_Pmean + ensLGM_Tmean
'
sem.fit.2 <- sem(Frugivores_Neo_Model_2, data=Frugivores_Neo.dat)
summary(sem.fit.2, stand=T, rsq=T, fit.measures=T)


# delete ensLGM_Pmean from the body size model
Frugivores_Neo_Model_3 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + ensLGM_Pmean + PC1 + PC2 + PC3
BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tmean
NPP_mean ~ Soil + PC1 + PC3 + ensLGM_Pmean + ensLGM_Tmean
'
sem.fit.3 <- sem(Frugivores_Neo_Model_3, data=Frugivores_Neo.dat)
summary(sem.fit.3, stand=T, rsq=T, fit.measures=T)


# delete PC3 from the NPP model
Frugivores_Neo_Model_4 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + ensLGM_Pmean + PC1 + PC2 + PC3
BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tmean
NPP_mean ~ Soil + PC1 + ensLGM_Pmean + ensLGM_Tmean
'
sem.fit.4 <- sem(Frugivores_Neo_Model_4, data=Frugivores_Neo.dat)
summary(sem.fit.4, stand=T, rsq=T, fit.measures=T)


# delete ensLGM_Tmean from the NPP model
Frugivores_Neo_Model_5 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + ensLGM_Pmean + PC1 + PC2 + PC3
BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tmean
NPP_mean ~ Soil + PC1 + ensLGM_Pmean
'
sem.fit.5 <- sem(Frugivores_Neo_Model_5, data=Frugivores_Neo.dat)
summary(sem.fit.5, stand=T, rsq=T, fit.measures=T)


# delete Soil from the fruit size model
Frugivores_Neo_Model_6 <- 'FruitSize ~ BodySize + NPP_mean + ensLGM_Tmean + ensLGM_Pmean + PC1 + PC2 + PC3
BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tmean
NPP_mean ~ Soil + PC1 + ensLGM_Pmean
'
sem.fit.6 <- sem(Frugivores_Neo_Model_6, data=Frugivores_Neo.dat)
summary(sem.fit.6, stand=T, rsq=T, fit.measures=T)


# delete ensLGM_Pmean from the NPP model
Frugivores_Neo_Model_7 <- 'FruitSize ~ BodySize + NPP_mean + ensLGM_Tmean + ensLGM_Pmean + PC1 + PC2 + PC3
BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tmean
NPP_mean ~ Soil + PC1
'
sem.fit.7 <- sem(Frugivores_Neo_Model_7, data=Frugivores_Neo.dat)
summary(sem.fit.7, stand=T, rsq=T, fit.measures=T)


# delete PC1 from the NPP model
Frugivores_Neo_Model_8 <- 'FruitSize ~ BodySize + NPP_mean + ensLGM_Tmean + ensLGM_Pmean + PC1 + PC2 + PC3
BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tmean
NPP_mean ~ Soil
'
sem.fit.8 <- sem(Frugivores_Neo_Model_8, data=Frugivores_Neo.dat)
summary(sem.fit.8, stand=T, rsq=T, fit.measures=T)


# delete NPP from the fruit size model
Frugivores_Neo_Model_9 <- 'FruitSize ~ BodySize + ensLGM_Tmean + ensLGM_Pmean + PC1 + PC2 + PC3
BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tmean
NPP_mean ~ Soil
'
sem.fit.9 <- sem(Frugivores_Neo_Model_9, data=Frugivores_Neo.dat)
summary(sem.fit.9, stand=T, rsq=T, fit.measures=T)


# delete ensLGM_Pmean from the fruit size model
Frugivores_Neo_Model_10 <- 'FruitSize ~ BodySize + ensLGM_Tmean + PC1 + PC2 + PC3
BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tmean
NPP_mean ~ Soil
'
sem.fit.10 <- sem(Frugivores_Neo_Model_10, data=Frugivores_Neo.dat)
summary(sem.fit.10, stand=T, rsq=T, fit.measures=T)


# delete PC3 from the fuit size model
Frugivores_Neo_Model_11 <- 'FruitSize ~ BodySize + ensLGM_Tmean + PC1 + PC2
BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tmean
NPP_mean ~ Soil
'
sem.fit.11 <- sem(Frugivores_Neo_Model_11, data=Frugivores_Neo.dat)
summary(sem.fit.11, stand=T, rsq=T, fit.measures=T)


# delete ensLGM_Tmean from the fruit size model
Frugivores_Neo_Model_12 <- 'FruitSize ~ BodySize + PC1 + PC2
BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tmean
NPP_mean ~ Soil
'
sem.fit.12 <- sem(Frugivores_Neo_Model_12, data=Frugivores_Neo.dat)
summary(sem.fit.12, stand=T, rsq=T, fit.measures=T)


# delete PC2 from the fruit size model
Frugivores_Neo_Model_13<- 'FruitSize ~ BodySize + PC1
BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tmean
NPP_mean ~ Soil
'
sem.fit.13 <- sem(Frugivores_Neo_Model_13, data=Frugivores_Neo.dat)
summary(sem.fit.13, stand=T, rsq=T, fit.measures=T)






##### Afrotropics ######

Frugivores_Afro <- Frugivores [ which(Frugivores$REALM_LONG =='Afrotropics'), ]
Frugivores_Afro


# PCA of current climate variables
all_climate <- Frugivores_Afro[,12:17]
str(all_climate)
all_climate
PCA_clim<-prcomp(all_climate, scale=TRUE)
summary(PCA_clim)    # first 3 axes explain ~92%
biplot(PCA_clim, cex=0.9, xlim=c(-0.23,0.23))

PCA_clim$x[,1:3]
Frugivores_Neo<-cbind(Frugivores_Afro, PCA_clim$x[,1:3])
is.na(Frugivores_Afro)
str(Frugivores_Afro)
summary(Frugivores_Afro)


hist(log(Frugivores_Afro$FruitSize))    # palm fruit size
hist(log(Frugivores_Afro$BodySize))                # mammal body size
# Histograms (predictor variables)
hist(Frugivores_Afro$PC1)               # mean annual precipitation
hist(Frugivores_Afro$PC2)                # precipitation seasonlity
hist(Frugivores_Afro$PC3)             # mean annual temperature
hist(Frugivores_Afro$NPP_mean)               # NPP (net prim. production)
hist(Frugivores_Afro$ensLGM_Tmean)           # mean LGM ANOM Temp
hist(log(Frugivores_Afro$ensLGM_Pmean))      # mean LGM ANOM Prec
hist(Frugivores_Afro$Soil)


colnames(Frugivores_Afro)
str(Frugivores_Afro)
summary(Frugivores_Afro)

# adapt data
FruitSize <- log(Frugivores_Afro$FruitSize)
BodySize <- log(Frugivores_Afro$BodySize)
NPP_mean <- Frugivores_Afro$NPP_mean/10000
ensLGM_Tmean <- Frugivores_Afro$ensLGM_Tmean/100
ensLGM_Pmean <- log(Frugivores_Afro$ensLGM_Pmean)
Soil <- Frugivores_Afro$Soil/10
PC1 <- Frugivores_Afro$PC1
PC2 <- Frugivores_Afro$PC2
PC3 <- Frugivores_Afro$PC3

Frugivores_Afro.dat <- data.frame(FruitSize, BodySize, NPP_mean, ensLGM_Tmean, ensLGM_Pmean, Soil, PC1, PC2, PC3)
summary(Frugivores_Afro.dat)
varTable(Frugivores_Afro.dat)


# SEM approach
Frugivores_Afro_Model_1 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + ensLGM_Pmean + PC1 + PC2 + PC3
BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean
NPP_mean ~ Soil + PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean
'
sem.fit.1 <- sem(Frugivores_Afro_Model_1, data=Frugivores_Afro.dat)
summary(sem.fit.1, stand=T, rsq=T, fit.measures=T)

modindices(sem.fit.1)
resid(sem.fit.1, type="standardized")


# delete ensLGM_Pmean from the NPP model
Frugivores_Afro_Model_2 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + ensLGM_Pmean + PC1 + PC2 + PC3
BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean
NPP_mean ~ Soil + PC1 + PC2 + PC3 + ensLGM_Tmean
'
sem.fit.2 <- sem(Frugivores_Afro_Model_2, data=Frugivores_Afro.dat)
summary(sem.fit.2, stand=T, rsq=T, fit.measures=T)


# delete PC1 from the NPP model
Frugivores_Afro_Model_3 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + ensLGM_Pmean + PC1 + PC2 + PC3
BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean
NPP_mean ~ Soil + PC2 + PC3 + ensLGM_Tmean
'
sem.fit.3 <- sem(Frugivores_Afro_Model_3, data=Frugivores_Afro.dat)
summary(sem.fit.3, stand=T, rsq=T, fit.measures=T)


# delete NPP from the fruit size model
Frugivores_Afro_Model_4 <- 'FruitSize ~ BodySize + Soil + ensLGM_Tmean + ensLGM_Pmean + PC1 + PC2 + PC3
BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean
NPP_mean ~ Soil + PC2 + PC3 + ensLGM_Tmean
'
sem.fit.4 <- sem(Frugivores_Afro_Model_4, data=Frugivores_Afro.dat)
summary(sem.fit.4, stand=T, rsq=T, fit.measures=T)


# delete ensLGM_Pmean from the body size model
Frugivores_Afro_Model_5 <- 'FruitSize ~ BodySize + Soil + ensLGM_Tmean + ensLGM_Pmean + PC1 + PC2 + PC3
BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tmean
NPP_mean ~ Soil + PC2 + PC3 + ensLGM_Tmean
'
sem.fit.5 <- sem(Frugivores_Afro_Model_5, data=Frugivores_Afro.dat)
summary(sem.fit.5, stand=T, rsq=T, fit.measures=T)


# delete PC3 from the body size model
Frugivores_Afro_Model_6 <- 'FruitSize ~ BodySize + Soil + ensLGM_Tmean + ensLGM_Pmean + PC1 + PC2 + PC3
BodySize ~ NPP_mean + PC1 + PC2 + ensLGM_Tmean
NPP_mean ~ Soil + PC2 + PC3 + ensLGM_Tmean
'
sem.fit.6 <- sem(Frugivores_Afro_Model_6, data=Frugivores_Afro.dat)
summary(sem.fit.6, stand=T, rsq=T, fit.measures=T)


# delete Soil from the NPP model
Frugivores_Afro_Model_7 <- 'FruitSize ~ BodySize + Soil + ensLGM_Tmean + ensLGM_Pmean + PC1 + PC2 + PC3
BodySize ~ NPP_mean + PC1 + PC2 + ensLGM_Tmean
NPP_mean ~ PC2 + PC3 + ensLGM_Tmean
'
sem.fit.7 <- sem(Frugivores_Afro_Model_7, data=Frugivores_Afro.dat)
summary(sem.fit.7, stand=T, rsq=T, fit.measures=T)


# delete ensLGM_Tmean from the fruit size model
Frugivores_Afro_Model_8 <- 'FruitSize ~ BodySize + Soil + ensLGM_Pmean + PC1 + PC2 + PC3
BodySize ~ NPP_mean + PC1 + PC2 + ensLGM_Tmean
NPP_mean ~ PC2 + PC3 + ensLGM_Tmean
'
sem.fit.8 <- sem(Frugivores_Afro_Model_8, data=Frugivores_Afro.dat)
summary(sem.fit.8, stand=T, rsq=T, fit.measures=T)


# delete PC3 from the fruit size model
Frugivores_Afro_Model_9 <- 'FruitSize ~ BodySize + Soil + ensLGM_Pmean + PC1 + PC2
BodySize ~ NPP_mean + PC1 + PC2 + ensLGM_Tmean
NPP_mean ~ PC2 + PC3 + ensLGM_Tmean
'
sem.fit.9 <- sem(Frugivores_Afro_Model_9, data=Frugivores_Afro.dat)
summary(sem.fit.9, stand=T, rsq=T, fit.measures=T)


# delete ensLGM_Pmean from the fruit size model
Frugivores_Afro_Model_10 <- 'FruitSize ~ BodySize + Soil + PC1 + PC2
BodySize ~ NPP_mean + PC1 + PC2 + ensLGM_Tmean
NPP_mean ~ PC2 + PC3 + ensLGM_Tmean
'
sem.fit.10 <- sem(Frugivores_Afro_Model_10, data=Frugivores_Afro.dat)
summary(sem.fit.10, stand=T, rsq=T, fit.measures=T)


# delete Soil from the fruit size model
Frugivores_Afro_Model_11 <- 'FruitSize ~ BodySize + PC1 + PC2
BodySize ~ NPP_mean + PC1 + PC2 + ensLGM_Tmean
NPP_mean ~ PC2 + PC3 + ensLGM_Tmean
'
sem.fit.11 <- sem(Frugivores_Afro_Model_11, data=Frugivores_Afro.dat)
summary(sem.fit.11, stand=T, rsq=T, fit.measures=T)


# delete PC1 from the body size model
Frugivores_Afro_Model_12 <- 'FruitSize ~ BodySize + PC2
BodySize ~ NPP_mean + PC1 + PC2 + ensLGM_Tmean
NPP_mean ~ PC2 + PC3 + ensLGM_Tmean
'
sem.fit.12 <- sem(Frugivores_Afro_Model_12, data=Frugivores_Afro.dat)
summary(sem.fit.12, stand=T, rsq=T, fit.measures=T)


# delete NPP from the body size model
Frugivores_Afro_Model_13 <- 'FruitSize ~ BodySize + PC2
BodySize ~ PC1 + PC2 + ensLGM_Tmean
NPP_mean ~ PC2 + PC3 + ensLGM_Tmean
'
sem.fit.13 <- sem(Frugivores_Afro_Model_13, data=Frugivores_Afro.dat)
summary(sem.fit.13, stand=T, rsq=T, fit.measures=T)


# delete PC1 from the body size model
Frugivores_Afro_Model_14 <- 'FruitSize ~ BodySize + PC2
BodySize ~ PC2 + ensLGM_Tmean
NPP_mean ~ PC2 + PC3 + ensLGM_Tmean
'
sem.fit.14 <- sem(Frugivores_Afro_Model_14, data=Frugivores_Afro.dat)
summary(sem.fit.14, stand=T, rsq=T, fit.measures=T)


# delete ensLGM_Tmean from the body size model
Frugivores_Afro_Model_15 <- 'FruitSize ~ BodySize + PC2
BodySize ~ PC2
NPP_mean ~ PC2 + PC3 + ensLGM_Tmean
'
sem.fit.15 <- sem(Frugivores_Afro_Model_15, data=Frugivores_Afro.dat)
summary(sem.fit.15, stand=T, rsq=T, fit.measures=T)






# SE-ASIA
Frugivores_Austral <- Frugivores [ which(Frugivores$REALM_LONG =='Australasia'), ]
Frugivores_Austral

Frugivores_Indo <- Frugivores [ which(Frugivores$REALM_LONG =='IndoMalay'), ]
Frugivores_Indo

Frugivores_SEAsia <- rbind (Frugivores_Indo, Frugivores_Austral)
Frugivores_SEAsia


# PCA of current climate variables
all_climate <- Frugivores_SEAsia[,12:17]
names(all_climate)
str(all_climate)
all_climate
PCA_clim<-prcomp(all_climate, scale=TRUE)
summary(PCA_clim)    # first 3 axes explain ~96%
biplot(PCA_clim, cex=0.9, xlim=c(-0.23,0.23))

PCA_clim$x[,1:3]
Frugivores_SEAsia<-cbind(Frugivores_SEAsia, PCA_clim$x[,1:3])
is.na(Frugivores_SEAsia)
str(Frugivores_SEAsia)
summary(Frugivores_SEAsia)



colnames(Frugivores_SEAsia)
str(Frugivores_SEAsia)
summary(Frugivores_SEAsia)

# adapt data
FruitSize <- log(Frugivores_SEAsia$FruitSize)
BodySize <- log(Frugivores_SEAsia$BodySize)
NPP_mean <- Frugivores_SEAsia$NPP_mean/10000
ensLGM_Tmean <- Frugivores_SEAsia$ensLGM_Tmean/100
ensLGM_Pmean <- log(Frugivores_SEAsia$ensLGM_Pmean)
Soil <- Frugivores_SEAsia$Soil/10
PC1 <- Frugivores_SEAsia$PC1
PC2 <- Frugivores_SEAsia$PC2
PC3 <- Frugivores_SEAsia$PC3

Frugivores_SEAsia.dat <- data.frame(FruitSize, BodySize, NPP_mean, ensLGM_Tmean, ensLGM_Pmean, Soil, PC1, PC2, PC3)
summary(Frugivores_SEAsia.dat)
varTable(Frugivores_SEAsia.dat)


# SEM approach
Frugivores_SEA_Model_1 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + ensLGM_Pmean + PC1 + PC2 + PC3
BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean
NPP_mean ~ Soil + PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean
'
sem.fit.1 <- sem(Frugivores_SEA_Model_1, data=Frugivores_SEAsia.dat)
summary(sem.fit.1, stand=T, rsq=T, fit.measures=T)


# delete PC2 from the fruit size model
Frugivores_SEA_Model_2 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + ensLGM_Pmean + PC1 + PC3
BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean
NPP_mean ~ Soil + PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean
'
sem.fit.2 <- sem(Frugivores_SEA_Model_2, data=Frugivores_SEAsia.dat)
summary(sem.fit.2, stand=T, rsq=T, fit.measures=T)


# delete NPP from the body size model
Frugivores_SEA_Model_3 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + ensLGM_Pmean + PC1 + PC3
BodySize ~ PC1 + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean
NPP_mean ~ Soil + PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean
'
sem.fit.3 <- sem(Frugivores_SEA_Model_3, data=Frugivores_SEAsia.dat)
summary(sem.fit.3, stand=T, rsq=T, fit.measures=T)


# delete PC2 from the NPP model
Frugivores_SEA_Model_4 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + ensLGM_Pmean + PC1 + PC3
BodySize ~ PC1 + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean
NPP_mean ~ Soil + PC1 + PC3 + ensLGM_Pmean + ensLGM_Tmean
'
sem.fit.4 <- sem(Frugivores_SEA_Model_4, data=Frugivores_SEAsia.dat)
summary(sem.fit.4, stand=T, rsq=T, fit.measures=T)


# delete PC3 from the NPP model
Frugivores_SEA_Model_5 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + ensLGM_Pmean + PC1 + PC3
BodySize ~ PC1 + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean
NPP_mean ~ Soil + PC1 + ensLGM_Pmean + ensLGM_Tmean
'
sem.fit.5 <- sem(Frugivores_SEA_Model_5, data=Frugivores_SEAsia.dat)
summary(sem.fit.5, stand=T, rsq=T, fit.measures=T)


# delete ensLGM_Tmean from the fruit size model
Frugivores_SEA_Model_6 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Pmean + PC1 + PC3
BodySize ~ PC1 + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean
NPP_mean ~ Soil + PC1 + ensLGM_Pmean + ensLGM_Tmean
'
sem.fit.6 <- sem(Frugivores_SEA_Model_6, data=Frugivores_SEAsia.dat)
summary(sem.fit.6, stand=T, rsq=T, fit.measures=T)


# delete ensLGM_Tmean from the NPP model
Frugivores_SEA_Model_7 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Pmean + PC1 + PC3
BodySize ~ PC1 + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean
NPP_mean ~ Soil + PC1 + ensLGM_Pmean
'
sem.fit.7 <- sem(Frugivores_SEA_Model_7, data=Frugivores_SEAsia.dat)
summary(sem.fit.7, stand=T, rsq=T, fit.measures=T)


# delete PC1 from the fruit size model
Frugivores_SEA_Model_8 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Pmean + PC3
BodySize ~ PC1 + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean
NPP_mean ~ Soil + PC1 + ensLGM_Pmean
'
sem.fit.8 <- sem(Frugivores_SEA_Model_8, data=Frugivores_SEAsia.dat)
summary(sem.fit.8, stand=T, rsq=T, fit.measures=T)


# delete ensLGM_Pmean from the body size model
Frugivores_SEA_Model_9 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Pmean + PC3
BodySize ~ PC1 + PC2 + PC3 + ensLGM_Tmean
NPP_mean ~ Soil + PC1 + ensLGM_Pmean
'
sem.fit.9 <- sem(Frugivores_SEA_Model_9, data=Frugivores_SEAsia.dat)
summary(sem.fit.9, stand=T, rsq=T, fit.measures=T)


# delete Soil from the NPP model
Frugivores_SEA_Model_10 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Pmean + PC3
BodySize ~ PC1 + PC2 + PC3 + ensLGM_Tmean
NPP_mean ~ PC1 + ensLGM_Pmean
'
sem.fit.10 <- sem(Frugivores_SEA_Model_10, data=Frugivores_SEAsia.dat)
summary(sem.fit.10, stand=T, rsq=T, fit.measures=T)


# delete Soil from the fruit size model
Frugivores_SEA_Model_11 <- 'FruitSize ~ BodySize + NPP_mean + ensLGM_Pmean + PC3
BodySize ~ PC1 + PC2 + PC3 + ensLGM_Tmean
NPP_mean ~ PC1 + ensLGM_Pmean
'
sem.fit.11 <- sem(Frugivores_SEA_Model_11, data=Frugivores_SEAsia.dat)
summary(sem.fit.11, stand=T, rsq=T, fit.measures=T)


# delete ensLGM_Pmean from the fruit size model
Frugivores_SEA_Model_12 <- 'FruitSize ~ BodySize + NPP_mean + PC3
BodySize ~ PC1 + PC2 + PC3 + ensLGM_Tmean
NPP_mean ~ PC1 + ensLGM_Pmean
'
sem.fit.12 <- sem(Frugivores_SEA_Model_12, data=Frugivores_SEAsia.dat)
summary(sem.fit.12, stand=T, rsq=T, fit.measures=T)


# delete NPP from the fruit size model
Frugivores_SEA_Model_13 <- 'FruitSize ~ BodySize + PC3
BodySize ~ PC1 + PC2 + PC3 + ensLGM_Tmean
NPP_mean ~ PC1 + ensLGM_Pmean
'
sem.fit.13 <- sem(Frugivores_SEA_Model_13, data=Frugivores_SEAsia.dat)
summary(sem.fit.13, stand=T, rsq=T, fit.measures=T)    # almost finished here: more kicked out just because of p (PC2) = 0.052


# delete PC2 from the body size model 
Frugivores_SEA_Model_14 <- 'FruitSize ~ BodySize + PC3
BodySize ~ PC1 + PC3 + ensLGM_Tmean
NPP_mean ~ PC1 + ensLGM_Pmean
'
sem.fit.14 <- sem(Frugivores_SEA_Model_14, data=Frugivores_SEAsia.dat)
summary(sem.fit.14, stand=T, rsq=T, fit.measures=T)


# delete PC1 from the body size model
Frugivores_SEA_Model_14 <- 'FruitSize ~ BodySize + PC3
BodySize ~ PC3 + ensLGM_Tmean
NPP_mean ~ PC1 + ensLGM_Pmean
'
sem.fit.14 <- sem(Frugivores_SEA_Model_14, data=Frugivores_SEAsia.dat)
summary(sem.fit.14, stand=T, rsq=T, fit.measures=T)


# delete ensLGM_Tmean from the body size model
Frugivores_SEA_Model_15 <- 'FruitSize ~ BodySize + PC3
BodySize ~ PC3
NPP_mean ~ PC1 + ensLGM_Pmean
'
sem.fit.15 <- sem(Frugivores_SEA_Model_15, data=Frugivores_SEAsia.dat)
summary(sem.fit.15, stand=T, rsq=T, fit.measures=T)
