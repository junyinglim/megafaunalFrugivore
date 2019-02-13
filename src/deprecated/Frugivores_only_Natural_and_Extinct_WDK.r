# Bastian Göldel
# Project: Global Palm and Mammals Functional traits distributions and relation  (category fruits = all frugivores)
# only extinct frugivores & present frugivores with natural ranges from Søren (without rest of current species)

# additions by WDK 18 Sept 2015

# libraries
#library("lavaan", lib="P:/Documents/Uni Aarhus/PhD/data from Daniel/Data for Bastian/NewWorld")
library ("lavaan")

#Workspace
rm(list=ls())
#setwd("C:/Analysis_Bastian/Mammals and Palms/Office")
setwd("D:/Data/wkissli1/My Documents/Projects/Global/06_Palms/02_Bastian/Data/SEM")


######## CATEGORY FRUIT = ALL Frugivores #########


#Load data
Frugivores<-read.table("ALLINFO_PALMS_MAMMALS_ENVIRON.txt", header=T,sep="\t")
str(Frugivores)
names(Frugivores)
Frugivores <- na.omit(Frugivores) # ignores NAs
is.na(Frugivores)                      # check
str(Frugivores)
names(Frugivores)
summary(Frugivores)


# Histograms (trait variables)
hist(log(Frugivores$FruitSize))    # palm fruit size

hist(Frugivores$BodySize)                # mammal body size

# Histograms (predictor variables)
hist(Frugivores$PREC_Sum)               # mean annual precipitation
hist(Frugivores$PREC_CV)                # precipitation seasonlity
hist(Frugivores$Tmean_mean)             # mean annual temperature
hist(log(Frugivores$Temp_SD))           # temperature seasonality
hist(Frugivores$NPP_mean)               # NPP (net prim. production)
hist(Frugivores$ensLGM_Tmean)           # mean LGM ANOM Temp
hist(log(Frugivores$ensLGM_Pmean))      # mean LGM ANOM Prec
hist(Frugivores$Soil)                   # % of sandy soil (depth 45cm)


# PCA of current climate variables
all_climate <- Frugivores[,15:21]
names (all_climate)
str(all_climate)
all_climate
PCA_clim<-prcomp(all_climate, scale=TRUE)
summary(PCA_clim)    # first 3 axes explain ~93.8%
biplot(PCA_clim, cex=0.9, xlim=c(-0.23,0.23))

###WDK: the thrid axis doesn't add too much, but we can leave it in for now
#third axis
biplot(PCA_clim, choices = 2:3, cex=0.9, xlim=c(-0.23,0.23))
biplot(PCA_clim, choices = c(1,3), cex=0.9, xlim=c(-0.23,0.23))

PCA_clim$x[,1:3]
Frugivores<-cbind(Frugivores, PCA_clim$x[,1:3])
is.na(Frugivores)
str(Frugivores)
summary(Frugivores)

head(PCA_clim$rotation)
# PC1: strong positive with T_cold_quart (0.46), with Tmean_mean (0.41) and PREC_Sum (0.36), negative with Temp_SD (-0.44)
# PC2: strong positive with PREC_CV (0.62), negative with P_drie_quart (-0.56), and PREC_Sum (-0.38)
# PC3: strong positive with Tmean_mean (0.46) and Temp_SD (0.41). negative with PREC_CV (-0.55) and PREC_Sum (-0.51)


# Histograms PCA axes
hist(Frugivores$PC1)
hist(Frugivores$PC2)
hist(Frugivores$PC3)


# plot single relationships to check for non-linear response
plot(log(Frugivores$FruitSize) ~ Frugivores$BodySize)
plot(log(Frugivores$FruitSize) ~ Frugivores$PC1)
plot(log(Frugivores$FruitSize) ~ Frugivores$PC2)
plot(log(Frugivores$FruitSize) ~ Frugivores$PC3)
plot(log(Frugivores$FruitSize) ~ Frugivores$NPP_mean)
plot(log(Frugivores$FruitSize) ~ Frugivores$ensLGM_Tmean)
plot(log(Frugivores$FruitSize) ~ Frugivores$ensLGM_Pmean)
plot(log(Frugivores$FruitSize) ~ Frugivores$Soil)

plot(Frugivores$BodySize ~ Frugivores$PC1)
plot(Frugivores$BodySize ~ Frugivores$PC2)
plot(Frugivores$BodySize ~ Frugivores$PC3)
plot(Frugivores$BodySize ~ Frugivores$NPP_mean)
plot(Frugivores$BodySize ~ Frugivores$ensLGM_Tmean)
plot(Frugivores$BodySize ~ Frugivores$ensLGM_Pmean)
plot(Frugivores$BodySize ~ Frugivores$Soil)



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
#WDK: I think we need to explain this in the manuscript text. 

FruitSize <- log(Frugivores$FruitSize)
BodySize <- Frugivores$BodySize/1000000
NPP_mean <- Frugivores$NPP_mean/10000
ensLGM_Tmean <- Frugivores$ensLGM_Tmean/100
ensLGM_Pmean <- Frugivores$ensLGM_Pmean/1000
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


# delete PC3 from the Body size model
Frugivores_Model_1 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + ensLGM_Pmean + PC1 + PC2 + PC3
BodySize ~ NPP_mean + PC1 + PC2 + ensLGM_Tmean + ensLGM_Pmean
NPP_mean ~ Soil + PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean
'
sem.fit.1 <- sem(Frugivores_Model_1, data=Frugivores.dat)
summary(sem.fit.1, stand=T, rsq=T, fit.measures=T)


# delete ensLGM_Pmean from the body size model
Frugivores_Model_3 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + ensLGM_Pmean + PC1 + PC2 + PC3
BodySize ~ NPP_mean + PC1 + PC2 + ensLGM_Tmean
NPP_mean ~ Soil + PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean
'
sem.fit.3 <- sem(Frugivores_Model_3, data=Frugivores.dat)
summary(sem.fit.3, stand=T, rsq=T, fit.measures=T)


# delete PC2 from the Fruit size model
Frugivores_Model_4 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + ensLGM_Pmean + PC1 + PC3
BodySize ~ NPP_mean + PC1 + PC2 + ensLGM_Tmean
NPP_mean ~ Soil + PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean
'
sem.fit.4 <- sem(Frugivores_Model_4, data=Frugivores.dat)
summary(sem.fit.4, stand=T, rsq=T, fit.measures=T)


# delete ensLGM_Pmean from the Fruit size model
Frugivores_Model_5 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + PC1 + PC3
BodySize ~ NPP_mean + PC1 + PC2 + ensLGM_Tmean
NPP_mean ~ Soil + PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean
'
sem.fit.5 <- sem(Frugivores_Model_5, data=Frugivores.dat)
summary(sem.fit.5, stand=T, rsq=T, fit.measures=T)


# delete ensLGM_Tmean from the NPP model
Frugivores_Model_6 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + PC1 + PC3
BodySize ~ NPP_mean + PC1 + PC2 + ensLGM_Tmean
NPP_mean ~ Soil + PC1 + PC2 + PC3 + ensLGM_Pmean
'
sem.fit.6 <- sem(Frugivores_Model_6, data=Frugivores.dat)
summary(sem.fit.6, stand=T, rsq=T, fit.measures=T)


# delete PC1 from the Fruit size model
Frugivores_Model_7 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + PC3
BodySize ~ NPP_mean + PC1 + PC2 + ensLGM_Tmean
NPP_mean ~ Soil + PC1 + PC2 + PC3 + ensLGM_Pmean
'
sem.fit.7 <- sem(Frugivores_Model_7, data=Frugivores.dat)
summary(sem.fit.7, stand=T, rsq=T, fit.measures=T)


# delete PC1 from the body size model
Frugivores_Model_8 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + PC3
BodySize ~ NPP_mean + PC2 + ensLGM_Tmean
NPP_mean ~ Soil + PC1 + PC2 + PC3 + ensLGM_Pmean
'
sem.fit.8 <- sem(Frugivores_Model_8, data=Frugivores.dat)
summary(sem.fit.8, stand=T, rsq=T, fit.measures=T)


# delete body size from the fruit size model
Frugivores_Model_9 <- 'FruitSize ~ Soil + NPP_mean + ensLGM_Tmean + PC3
BodySize ~ NPP_mean + PC2 + ensLGM_Tmean
NPP_mean ~ Soil + PC1 + PC2 + PC3 + ensLGM_Pmean
'
sem.fit.9 <- sem(Frugivores_Model_9, data=Frugivores.dat)
summary(sem.fit.9, stand=T, rsq=T, fit.measures=T)


# delete ensLGM_Tmean from the body size model
Frugivores_Model_10 <- 'FruitSize ~ Soil + NPP_mean + ensLGM_Tmean + PC3
BodySize ~ NPP_mean + PC2
NPP_mean ~ Soil + PC1 + PC2 + PC3 + ensLGM_Pmean
'
sem.fit.10 <- sem(Frugivores_Model_10, data=Frugivores.dat)
summary(sem.fit.10, stand=T, rsq=T, fit.measures=T)
# FINAL model


#WDK: so, there is no body size effect on fruit size in the global model, right?



##### Neotropics ######

Frugivores_Neo <- Frugivores [ which(Frugivores$REALM_LONG =='Neotropics'), ]
str(Frugivores_Neo)

# PCA of current climate variables
all_climate <- Frugivores_Neo[,15:21]
str(all_climate)
all_climate
PCA_clim<-prcomp(all_climate, scale=TRUE)
summary(PCA_clim)    # first 3 axes explain ~96%
biplot(PCA_clim, cex=0.9, xlim=c(-0.23,0.23))

PCA_clim$x[,1:3]
Frugivores_Neo<-cbind(Frugivores_Neo, PCA_clim$x[,1:3])
is.na(Frugivores_Neo)
str(Frugivores_Neo)
summary(Frugivores_Neo)


hist(log(Frugivores_Neo$FruitSize))    # palm fruit size
hist(log(Frugivores_Neo$BodySize))                # mammal body size
# Histograms (predictor variables)
hist(Frugivores_Neo$PC1)               # mean annual precipitation
hist(Frugivores_Neo$PC2)                # precipitation seasonlity
hist(Frugivores_Neo$PC3)             # mean annual temperature
hist(Frugivores_Neo$NPP_mean)               # NPP (net prim. production)
hist(Frugivores_Neo$ensLGM_Tmean)           # mean LGM ANOM Temp
hist(log(Frugivores_Neo$ensLGM_Pmean))      # mean LGM ANOM Prec
hist(log(Frugivores_Neo$Soil))


colnames(Frugivores_Neo)
str(Frugivores_Neo)
summary(Frugivores_Neo)

# adapt data
FruitSize <- log(Frugivores_Neo$FruitSize)
BodySize <- log(Frugivores_Neo$BodySize)
NPP_mean <- Frugivores_Neo$NPP_mean/10000
ensLGM_Tmean <- log(Frugivores_Neo$ensLGM_Tmean)
ensLGM_Pmean <- log(Frugivores_Neo$ensLGM_Pmean)
Soil <- log(Frugivores_Neo$Soil)
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


# delete PC3 from the body size model
Frugivores_Neo_Model_2 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + ensLGM_Pmean + PC1 + PC2 + PC3
BodySize ~ NPP_mean + PC1 + PC2 + ensLGM_Tmean + ensLGM_Pmean
NPP_mean ~ Soil + PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean
'
sem.fit.2 <- sem(Frugivores_Neo_Model_2, data=Frugivores_Neo.dat)
summary(sem.fit.2, stand=T, rsq=T, fit.measures=T)


# delete ensLGM_Pmean from the NPP model
Frugivores_Neo_Model_3 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + ensLGM_Pmean + PC1 + PC2 + PC3
BodySize ~ NPP_mean + PC1 + PC2 + ensLGM_Tmean + ensLGM_Pmean
NPP_mean ~ Soil + PC1 + PC2 + PC3 + ensLGM_Tmean
'
sem.fit.3 <- sem(Frugivores_Neo_Model_3, data=Frugivores_Neo.dat)
summary(sem.fit.3, stand=T, rsq=T, fit.measures=T) 


# delete ensLGM_Tmean from the body size model
Frugivores_Neo_Model_4 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + ensLGM_Pmean + PC1 + PC2 + PC3
BodySize ~ NPP_mean + PC1 + PC2 + ensLGM_Pmean
NPP_mean ~ Soil + PC1 + PC2 + PC3 + ensLGM_Tmean
'
sem.fit.4 <- sem(Frugivores_Neo_Model_4, data=Frugivores_Neo.dat)
summary(sem.fit.4, stand=T, rsq=T, fit.measures=T)


# delete PC2 from the NPP model
Frugivores_Neo_Model_5 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + ensLGM_Pmean + PC1 + PC2 + PC3
BodySize ~ NPP_mean + PC1 + PC2 + ensLGM_Pmean
NPP_mean ~ Soil + PC1 + PC3 + ensLGM_Tmean
'
sem.fit.5 <- sem(Frugivores_Neo_Model_5, data=Frugivores_Neo.dat)
summary(sem.fit.5, stand=T, rsq=T, fit.measures=T)


# delete PC2 from the fruit size model
Frugivores_Neo_Model_6 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + ensLGM_Pmean + PC1 + PC3
BodySize ~ NPP_mean + PC1 + PC2 + ensLGM_Pmean
NPP_mean ~ Soil + PC1 + PC3 + ensLGM_Tmean
'
sem.fit.6 <- sem(Frugivores_Neo_Model_6, data=Frugivores_Neo.dat)
summary(sem.fit.6, stand=T, rsq=T, fit.measures=T)


# delete PC3 from the fruit size model
Frugivores_Neo_Model_7 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + ensLGM_Pmean + PC1
BodySize ~ NPP_mean + PC1 + PC2 + ensLGM_Pmean
NPP_mean ~ Soil + PC1 + PC3 + ensLGM_Tmean
'
sem.fit.7 <- sem(Frugivores_Neo_Model_7, data=Frugivores_Neo.dat)
summary(sem.fit.7, stand=T, rsq=T, fit.measures=T)


# delete ensLGM_Tmean from the fruit size model
Frugivores_Neo_Model_8 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Pmean + PC1
BodySize ~ NPP_mean + PC1 + PC2 + ensLGM_Pmean
NPP_mean ~ Soil + PC1 + PC3 + ensLGM_Tmean
'
sem.fit.8 <- sem(Frugivores_Neo_Model_8, data=Frugivores_Neo.dat)
summary(sem.fit.8, stand=T, rsq=T, fit.measures=T)


# delete PC1 from the fruit size model
Frugivores_Neo_Model_9 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Pmean
BodySize ~ NPP_mean + PC1 + PC2 + ensLGM_Pmean
NPP_mean ~ Soil + PC1 + PC3 + ensLGM_Tmean
'
sem.fit.9 <- sem(Frugivores_Neo_Model_9, data=Frugivores_Neo.dat)
summary(sem.fit.9, stand=T, rsq=T, fit.measures=T)


# delete PC2 from the body size model
Frugivores_Neo_Model_10 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Pmean
BodySize ~ NPP_mean + PC1 + ensLGM_Pmean
NPP_mean ~ Soil + PC1 + PC3 + ensLGM_Tmean
'
sem.fit.10 <- sem(Frugivores_Neo_Model_10, data=Frugivores_Neo.dat)
summary(sem.fit.10, stand=T, rsq=T, fit.measures=T)
# FINAL MODEL!

#WDK: result: effect of extinct frugivore body size on fruit size in the Neotropics (cool!)



##### Afrotropics ######

Frugivores_Afro <- Frugivores [ which(Frugivores$REALM_LONG =='Afrotropics'), ]
Frugivores_Afro

# PCA of current climate variables
all_climate <- Frugivores_Afro[,15:21]
names (all_climate)
str(all_climate)
all_climate
PCA_clim<-prcomp(all_climate, scale=TRUE)
summary(PCA_clim)    # first 3 axes explain ~96%
biplot(PCA_clim, cex=0.9, xlim=c(-0.23,0.23))

PCA_clim$x[,1:3]
Frugivores_Afro<-cbind(Frugivores_Afro, PCA_clim$x[,1:3])
is.na(Frugivores_Afro)
str(Frugivores_Afro)
summary(Frugivores_Afro)


hist(log(Frugivores_Afro$FruitSize))    # palm fruit size
hist(Frugivores_Afro$BodySize)                # mammal body size
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
BodySize <- Frugivores_Afro$BodySize/1000000
NPP_mean <- Frugivores_Afro$NPP_mean/10000
ensLGM_Tmean <- Frugivores_Afro$ensLGM_Tmean/100
ensLGM_Pmean <- Frugivores_Afro$ensLGM_Pmean/1000
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


# delete ensLGM_Tmean from the fruit size model
Frugivores_Afro_Model_2 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Pmean + PC1 + PC2 + PC3
BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean
NPP_mean ~ Soil + PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean
'
sem.fit.2 <- sem(Frugivores_Afro_Model_2, data=Frugivores_Afro.dat)
summary(sem.fit.2, stand=T, rsq=T, fit.measures=T)


# delete PC3 from the body size model
Frugivores_Afro_Model_3 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Pmean + PC1 + PC2 + PC3
BodySize ~ NPP_mean + PC1 + PC2 + ensLGM_Tmean + ensLGM_Pmean
NPP_mean ~ Soil + PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean
'
sem.fit.3 <- sem(Frugivores_Afro_Model_3, data=Frugivores_Afro.dat)
summary(sem.fit.3, stand=T, rsq=T, fit.measures=T)


# delete Soil from the NPP model
Frugivores_Afro_Model_4 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Pmean + PC1 + PC2 + PC3
BodySize ~ NPP_mean + PC1 + PC2 + ensLGM_Tmean + ensLGM_Pmean
NPP_mean ~ PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean
'
sem.fit.4 <- sem(Frugivores_Afro_Model_4, data=Frugivores_Afro.dat)
summary(sem.fit.4, stand=T, rsq=T, fit.measures=T)


# delete PC1 from the body size model
Frugivores_Afro_Model_5 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Pmean + PC1 + PC2 + PC3
BodySize ~ NPP_mean + PC2 + ensLGM_Tmean + ensLGM_Pmean
NPP_mean ~ PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean
'
sem.fit.5 <- sem(Frugivores_Afro_Model_5, data=Frugivores_Afro.dat)
summary(sem.fit.5, stand=T, rsq=T, fit.measures=T)


# delete body size model fruit size model
Frugivores_Afro_Model_6 <- 'FruitSize ~ Soil + NPP_mean + ensLGM_Pmean + PC1 + PC2 + PC3
BodySize ~ NPP_mean + PC2 + ensLGM_Tmean + ensLGM_Pmean
NPP_mean ~ PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean
'
sem.fit.6 <- sem(Frugivores_Afro_Model_6, data=Frugivores_Afro.dat)
summary(sem.fit.6, stand=T, rsq=T, fit.measures=T)


# delete PC2 from the fruit size model
Frugivores_Afro_Model_7 <- 'FruitSize ~ Soil + NPP_mean + ensLGM_Pmean + PC1 + PC3
BodySize ~ NPP_mean + PC2 + ensLGM_Tmean + ensLGM_Pmean
NPP_mean ~ PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean
'
sem.fit.7 <- sem(Frugivores_Afro_Model_7, data=Frugivores_Afro.dat)
summary(sem.fit.7, stand=T, rsq=T, fit.measures=T)


# delete Soil from the fruit size model
Frugivores_Afro_Model_8 <- 'FruitSize ~ NPP_mean + ensLGM_Pmean + PC1 + PC3
BodySize ~ NPP_mean + PC2 + ensLGM_Tmean + ensLGM_Pmean
NPP_mean ~ PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean
'
sem.fit.8 <- sem(Frugivores_Afro_Model_8, data=Frugivores_Afro.dat)
summary(sem.fit.8, stand=T, rsq=T, fit.measures=T)


# WDK: I added this
# delete PC1 from the NPP model
Frugivores_Afro_Model_9 <- 'FruitSize ~ NPP_mean + ensLGM_Pmean + PC1 + PC3
BodySize ~ NPP_mean + PC2 + ensLGM_Tmean + ensLGM_Pmean
NPP_mean ~ PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean
'
sem.fit.9 <- sem(Frugivores_Afro_Model_9, data=Frugivores_Afro.dat)
summary(sem.fit.9, stand=T, rsq=T, fit.measures=T)

#WDK: I deleted this
# delete ensLGM_Tmean from the NPP model
#Frugivores_Afro_Model_9 <- 'FruitSize ~ NPP_mean + ensLGM_Pmean + PC1 + PC3
#BodySize ~ NPP_mean + PC2 + ensLGM_Tmean + ensLGM_Pmean
#NPP_mean ~ PC1 + PC2 + PC3 + ensLGM_Pmean
#'
#sem.fit.9 <- sem(Frugivores_Afro_Model_9, data=Frugivores_Afro.dat)
#summary(sem.fit.9, stand=T, rsq=T, fit.measures=T)


# delete PC2 from the body size model
Frugivores_Afro_Model_10 <- 'FruitSize ~ NPP_mean + ensLGM_Pmean + PC1 + PC3
BodySize ~ NPP_mean + ensLGM_Tmean + ensLGM_Pmean
NPP_mean ~ PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean
'
sem.fit.10 <- sem(Frugivores_Afro_Model_10, data=Frugivores_Afro.dat)
summary(sem.fit.10, stand=T, rsq=T, fit.measures=T)
# FINAL Model

#WDK: result: no body size effect of extinct frugivore species on fruit size (as expected)



###### INDOMALAY ########

Frugivores_Indo <- Frugivores [ which(Frugivores$REALM_LONG =='IndoMalay'), ]
Frugivores_Indo


# PCA of current climate variables
all_climate <- Frugivores_Indo[,15:21]
names (all_climate)
str(all_climate)
all_climate
PCA_clim<-prcomp(all_climate, scale=TRUE)
summary(PCA_clim)    # first 3 axes explain ~92%
biplot(PCA_clim, cex=0.9, xlim=c(-0.23,0.23))

PCA_clim$x[,1:3]
Frugivores_Indo<-cbind(Frugivores_Indo, PCA_clim$x[,1:3])
is.na(Frugivores_Indo)
str(Frugivores_Indo)
summary(Frugivores_Indo)



hist(log(Frugivores_Indo$FruitSize))    # palm fruit size
hist(Frugivores_Indo$BodySize)                # mammal body size
# Histograms (predictor variables)
hist(Frugivores_Indo$PC1)               # mean annual precipitation
hist(Frugivores_Indo$PC2)                # precipitation seasonlity
hist(Frugivores_Indo$PC3)             # mean annual temperature
hist(Frugivores_Indo$NPP_mean)               # NPP (net prim. production)
hist(Frugivores_Indo$ensLGM_Tmean)           # mean LGM ANOM Temp
hist(Frugivores_Indo$ensLGM_Pmean)      # mean LGM ANOM Prec
hist(Frugivores_Indo$Soil)      

 
colnames(Frugivores_Indo)
str(Frugivores_Indo)
summary(Frugivores_Indo)

# adapt data        
FruitSize <- log(Frugivores_Indo$FruitSize)            
BodySize <- Frugivores_Indo$BodySize/1000000
NPP_mean <- Frugivores_Indo$NPP_mean/10000         
ensLGM_Tmean <- Frugivores_Indo$ensLGM_Tmean/100          
ensLGM_Pmean <- Frugivores_Indo$ensLGM_Pmean/1000     
Soil <- Frugivores_Indo$Soil/10
PC1 <- Frugivores_Indo$PC1
PC2 <- Frugivores_Indo$PC2
PC3 <- Frugivores_Indo$PC3

Frugivores_Indo.dat <- data.frame(FruitSize, BodySize, NPP_mean, ensLGM_Tmean, ensLGM_Pmean, Soil, PC1, PC2, PC3)
summary(Frugivores_Indo.dat)
varTable(Frugivores_Indo.dat)

#plot(Frugivores_Indo.dat$FruitSize ~ Frugivores_Indo.dat$NPP_mean)

# SEM approach
Frugivores_Indo_Model_1 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + ensLGM_Pmean + PC1 + PC2 + PC3
BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean  
NPP_mean ~ Soil + PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean                    
'
sem.fit.1 <- sem(Frugivores_Indo_Model_1, data=Frugivores_Indo.dat)
summary(sem.fit.1, stand=T, rsq=T, fit.measures=T)    

modindices(sem.fit.1)                               
resid(sem.fit.1, type="standardized")

#WDK: changed to
# delete ensLGM_Pmean from the fruit size model
Frugivores_Indo_Model_2 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + PC1 + PC2 + PC3
BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean  
NPP_mean ~ Soil + PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean                    
'
sem.fit.2 <- sem(Frugivores_Indo_Model_2, data=Frugivores_Indo.dat)
summary(sem.fit.2, stand=T, rsq=T, fit.measures=T) 

#WDK: changed to
# delete Soil from the fruit size model
Frugivores_Indo_Model_3 <- 'FruitSize ~ BodySize  + NPP_mean + ensLGM_Tmean + PC1 + PC2 + PC3
BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean  
NPP_mean ~ Soil + PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean                    
'
sem.fit.3 <- sem(Frugivores_Indo_Model_3, data=Frugivores_Indo.dat)
summary(sem.fit.3, stand=T, rsq=T, fit.measures=T) 

#WDK: changed to
# delete PC3 from the NPP model
Frugivores_Indo_Model_4 <- 'FruitSize ~ BodySize  + NPP_mean + ensLGM_Tmean + PC1 + PC2 + PC3
BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean  
NPP_mean ~ Soil + PC1 + PC2 + ensLGM_Pmean + ensLGM_Tmean                    
'

sem.fit.4 <- sem(Frugivores_Indo_Model_4, data=Frugivores_Indo.dat)
summary(sem.fit.4, stand=T, rsq=T, fit.measures=T) 

#WDK: changed to
# delete PC2 from the NPP model
Frugivores_Indo_Model_5 <- 'FruitSize ~ BodySize  + NPP_mean + ensLGM_Tmean + PC1 + PC2 + PC3
BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean  
NPP_mean ~ Soil + PC1 + ensLGM_Pmean + ensLGM_Tmean                    
'

sem.fit.5 <- sem(Frugivores_Indo_Model_5, data=Frugivores_Indo.dat)
summary(sem.fit.5, stand=T, rsq=T, fit.measures=T) 

#WDK: changed to
# delete ensLGM_Tmean from the NPP model
Frugivores_Indo_Model_6 <- 'FruitSize ~ BodySize  + NPP_mean + ensLGM_Tmean + PC1 + PC2 + PC3
BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean  
NPP_mean ~ Soil + PC1 + ensLGM_Pmean                    
'

sem.fit.6 <- sem(Frugivores_Indo_Model_6, data=Frugivores_Indo.dat)
summary(sem.fit.6, stand=T, rsq=T, fit.measures=T) 

#WDK: changed to
# delete NPP from the fruit size model
Frugivores_Indo_Model_7 <- 'FruitSize ~ BodySize + ensLGM_Tmean + PC1 + PC2 + PC3
BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean  
NPP_mean ~ Soil + PC1 + ensLGM_Pmean                    
'

sem.fit.7 <- sem(Frugivores_Indo_Model_7, data=Frugivores_Indo.dat)
summary(sem.fit.7, stand=T, rsq=T, fit.measures=T) 

#WDK: changed to
# delete PC3 from the fruit size model
Frugivores_Indo_Model_8 <- 'FruitSize ~ BodySize + ensLGM_Tmean + PC1 + PC2
BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean  
NPP_mean ~ Soil + PC1 + ensLGM_Pmean                    
'

sem.fit.8 <- sem(Frugivores_Indo_Model_8, data=Frugivores_Indo.dat)
summary(sem.fit.8, stand=T, rsq=T, fit.measures=T) 

#WDK: changed to
# delete body size from the fruit size model
Frugivores_Indo_Model_9 <- 'FruitSize ~ ensLGM_Tmean + PC1 + PC2
BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean  
NPP_mean ~ Soil + PC1 + ensLGM_Pmean                    
'

sem.fit.9 <- sem(Frugivores_Indo_Model_9, data=Frugivores_Indo.dat)
summary(sem.fit.9, stand=T, rsq=T, fit.measures=T) 

#WDK: changed to
# delete ensLGM_Pmean from the body size model
Frugivores_Indo_Model_10 <- 'FruitSize ~ ensLGM_Tmean + PC1 + PC2
BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tmean  
NPP_mean ~ Soil + PC1 + ensLGM_Pmean                    
'

sem.fit.10 <- sem(Frugivores_Indo_Model_10, data=Frugivores_Indo.dat)
summary(sem.fit.10, stand=T, rsq=T, fit.measures=T) 


#WDK: changed to
# delete PC2 from the fruit size model
Frugivores_Indo_Model_11 <- 'FruitSize ~ ensLGM_Tmean + PC1
BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tmean  
NPP_mean ~ Soil + PC1 + ensLGM_Pmean                    
'

sem.fit.11 <- sem(Frugivores_Indo_Model_11, data=Frugivores_Indo.dat)
summary(sem.fit.11, stand=T, rsq=T, fit.measures=T) 

#modindices(sem.fit.11)  

#WDK: changed to
# delete PC1 from the fruit size model
Frugivores_Indo_Model_12 <- 'FruitSize ~ ensLGM_Tmean
BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tmean  
NPP_mean ~ Soil + PC1 + ensLGM_Pmean                    
'

sem.fit.12 <- sem(Frugivores_Indo_Model_12, data=Frugivores_Indo.dat)
summary(sem.fit.12, stand=T, rsq=T, fit.measures=T) 

#WDK: changed to
# delete ensLGM_Tmean from the fruit size model
Frugivores_Indo_Model_13 <- 'FruitSize
BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tmean  
NPP_mean ~ Soil + PC1 + ensLGM_Pmean                    
'

sem.fit.13 <- sem(Frugivores_Indo_Model_13, data=Frugivores_Indo.dat)
summary(sem.fit.13, stand=T, rsq=T, fit.measures=T) 

#WDK: result is a bit weird. None of the predictor variables stays in the model.
# Is this result also coimng out when using present-day frugivore distributions?

#WDK: here your previous implementation, the end result is the same
# delete ensLGM_Pmean from the body size model
#Frugivores_Indo_Model_2 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + ensLGM_Pmean + PC1 + PC2 + PC3
#BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tmean 
#NPP_mean ~ Soil + PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean                    
#'
#sem.fit.2 <- sem(Frugivores_Indo_Model_2, data=Frugivores_Indo.dat)
#summary(sem.fit.2, stand=T, rsq=T, fit.measures=T) 

# delete PC3 from the NPP model
#Frugivores_Indo_Model_3 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + ensLGM_Pmean + PC1 + PC2 + PC3
#BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tmean 
#NPP_mean ~ Soil + PC1 + PC2 + ensLGM_Pmean + ensLGM_Tmean                    
#'
#sem.fit.3 <- sem(Frugivores_Indo_Model_3, data=Frugivores_Indo.dat)
#summary(sem.fit.3, stand=T, rsq=T, fit.measures=T) 


# delete PC2 from the NPP model
#Frugivores_Indo_Model_4 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + ensLGM_Pmean + PC1 + PC2 + PC3
#BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tmean 
#NPP_mean ~ Soil + PC1 + ensLGM_Pmean + ensLGM_Tmean                    
#'
#sem.fit.4 <- sem(Frugivores_Indo_Model_4, data=Frugivores_Indo.dat)
#summary(sem.fit.4, stand=T, rsq=T, fit.measures=T) 

# delete ensLGM_Pmean from the fruit size model
#Frugivores_Indo_Model_5 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + PC1 + PC2 + PC3
#BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tmean 
#NPP_mean ~ Soil + PC1 + ensLGM_Pmean + ensLGM_Tmean                    
#'
#sem.fit.5 <- sem(Frugivores_Indo_Model_5, data=Frugivores_Indo.dat)
#summary(sem.fit.5, stand=T, rsq=T, fit.measures=T) 

# delete Soil from the fruit size model
#Frugivores_Indo_Model_6 <- 'FruitSize ~ BodySize + NPP_mean + ensLGM_Tmean + PC1 + PC2 + PC3
#BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tmean 
#NPP_mean ~ Soil + PC1 + ensLGM_Pmean + ensLGM_Tmean                    
#'
#sem.fit.6 <- sem(Frugivores_Indo_Model_6, data=Frugivores_Indo.dat)
#summary(sem.fit.6, stand=T, rsq=T, fit.measures=T) 


# delete ensLGM_Tmean from the NPP model
#Frugivores_Indo_Model_7 <- 'FruitSize ~ BodySize + NPP_mean + ensLGM_Tmean + PC1 + PC2 + PC3
#BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tmean 
#NPP_mean ~ Soil + PC1 + ensLGM_Pmean                  
#'
#sem.fit.7 <- sem(Frugivores_Indo_Model_7, data=Frugivores_Indo.dat)
#summary(sem.fit.7, stand=T, rsq=T, fit.measures=T) 

# delete NPP from the fruit size model
#Frugivores_Indo_Model_8 <- 'FruitSize ~ BodySize + ensLGM_Tmean + PC1 + PC2 + PC3
#BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tmean 
#NPP_mean ~ Soil + PC1 + ensLGM_Pmean                  
#'
#sem.fit.8 <- sem(Frugivores_Indo_Model_8, data=Frugivores_Indo.dat)
#summary(sem.fit.8, stand=T, rsq=T, fit.measures=T)

# delete PC3 from the fruit size model
#Frugivores_Indo_Model_9 <- 'FruitSize ~ BodySize + ensLGM_Tmean + PC1 + PC2
#BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tmean 
#NPP_mean ~ Soil + PC1 + ensLGM_Pmean                  
#'
#sem.fit.9 <- sem(Frugivores_Indo_Model_9, data=Frugivores_Indo.dat)
#summary(sem.fit.9, stand=T, rsq=T, fit.measures=T)

# delete body size from the fruit size model
#Frugivores_Indo_Model_10 <- 'FruitSize ~ ensLGM_Tmean + PC1 + PC2
#BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tmean 
#NPP_mean ~ Soil + PC1 + ensLGM_Pmean                  
#'
#sem.fit.10 <- sem(Frugivores_Indo_Model_10, data=Frugivores_Indo.dat)
#summary(sem.fit.10, stand=T, rsq=T, fit.measures=T)

# delete PC2 from the fruit size model
#Frugivores_Indo_Model_11 <- 'FruitSize ~ ensLGM_Tmean + PC1
#BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tmean 
#NPP_mean ~ Soil + PC1 + ensLGM_Pmean                  
#'
#sem.fit.11 <- sem(Frugivores_Indo_Model_11, data=Frugivores_Indo.dat)
#summary(sem.fit.11, stand=T, rsq=T, fit.measures=T)

# delete PC1 from the fruit size model
#Frugivores_Indo_Model_12 <- 'FruitSize ~ ensLGM_Tmean
#BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tmean 
#NPP_mean ~ Soil + PC1 + ensLGM_Pmean                  
#'
#sem.fit.12 <- sem(Frugivores_Indo_Model_12, data=Frugivores_Indo.dat)
#summary(sem.fit.12, stand=T, rsq=T, fit.measures=T)
# delete ensLGM_Tmean from the fruit size model --> fruit size model empty now!







# ###### AUSTRALASIA ########      (Warning: small number of observations)

#WDK: my suggestion: could we combine Australia with INDOMALAY? Then we would have three
# major REALMS (New World, Old World West, Old World East)

Frugivores_Austral <- Frugivores [ which(Frugivores$REALM_LONG =='Australasia'), ]
Frugivores_Austral

# PCA of current climate variables
all_climate <- Frugivores_Austral[,15:21]
names (all_climate)
str(all_climate)
all_climate
PCA_clim<-prcomp(all_climate, scale=TRUE)
summary(PCA_clim)    # first 3 axes explain ~99%
biplot(PCA_clim, cex=0.9, xlim=c(-0.23,0.23))

PCA_clim$x[,1:3]
Frugivores_Austral<-cbind(Frugivores_Austral, PCA_clim$x[,1:3])
is.na(Frugivores_Austral)
str(Frugivores_Austral)
summary(Frugivores_Austral)


hist(log(Frugivores_Austral$FruitSize))    # palm fruit size
hist(log(Frugivores_Austral$BodySize))                # mammal body size
# Histograms (predictor variables)
hist(Frugivores_Austral$PC1)               # mean annual precipitation
hist(Frugivores_Austral$PC2)                # precipitation seasonlity
hist(log(Frugivores_Austral$PC3))             # mean annual temperature
hist(log(Frugivores_Austral$NPP_mean))               # NPP (net prim. production)
hist(log(Frugivores_Austral$ensLGM_Tmean))           # mean LGM ANOM Temp
hist(log(Frugivores_Austral$ensLGM_Pmean))      # mean LGM ANOM Prec
hist(log(Frugivores_Austral$Soil))      

 
colnames(Frugivores_Austral)
str(Frugivores_Austral)
summary(Frugivores_Austral)

# adapt data        
FruitSize <- log(Frugivores_Austral$FruitSize)            
BodySize <- log(Frugivores_Austral$BodySize)
NPP_mean <- log(Frugivores_Austral$NPP_mean)        
ensLGM_Tmean <- log(Frugivores_Austral$ensLGM_Tmean)     
ensLGM_Pmean <- log(Frugivores_Austral$ensLGM_Pmean)     
Soil <- log(Frugivores_Austral$Soil)
PC1 <- Frugivores_Austral$PC1
PC2 <- Frugivores_Austral$PC2
PC3 <- log(Frugivores_Austral$PC3)

Frugivores_Austral.dat <- data.frame(FruitSize, BodySize, NPP_mean, ensLGM_Tmean, ensLGM_Pmean, Soil, PC1, PC2, PC3)
summary(Frugivores_Austral.dat)
varTable(Frugivores_Austral.dat)


# SEM approach
Frugivores_Austral_Model_1 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + ensLGM_Pmean + PC1 + PC2 + PC3
BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean  
NPP_mean ~ Soil + PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean                    
'
sem.fit.1 <- sem(Frugivores_Austral_Model_1, data=Frugivores_Austral.dat)
summary(sem.fit.1, stand=T, rsq=T, fit.measures=T)    

modindices(sem.fit.1)               
resid(sem.fit.1, type="standardized")



# delete ensLGM_Tmean from the fruit size model
Frugivores_Austral_Model_2 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Pmean + PC1 + PC2 + PC3
BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean  
NPP_mean ~ Soil + PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean                    
'
sem.fit.2 <- sem(Frugivores_Austral_Model_2, data=Frugivores_Austral.dat)
summary(sem.fit.2, stand=T, rsq=T, fit.measures=T) 


# delete Soil from the fruit size model
Frugivores_Austral_Model_3 <- 'FruitSize ~ BodySize + NPP_mean + ensLGM_Pmean + PC1 + PC2 + PC3
BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean  
NPP_mean ~ Soil + PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean                    
'
sem.fit.3 <- sem(Frugivores_Austral_Model_3, data=Frugivores_Austral.dat)
summary(sem.fit.3, stand=T, rsq=T, fit.measures=T)                
# would be final model






###### SE-Asia ########


Frugivores_SEAsia <- rbind (Frugivores_Indo, Frugivores_Austral)
Frugivores_SEAsia

# PCA of current climate variables
all_climate <- Frugivores_SEAsia[,15:21]
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



hist(log(Frugivores_SEAsia$FruitSize))    # palm fruit size
hist(Frugivores_SEAsia$BodySize)                # mammal body size
# Histograms (predictor variables)
hist(Frugivores_SEAsia$PC1)               # mean annual precipitation
hist(Frugivores_SEAsia$PC2)                # precipitation seasonlity
hist(Frugivores_SEAsia$PC3)             # mean annual temperature
hist(Frugivores_SEAsia$NPP_mean)               # NPP (net prim. production)
hist(Frugivores_SEAsia$ensLGM_Tmean)           # mean LGM ANOM Temp
hist(Frugivores_SEAsia$ensLGM_Pmean)     # mean LGM ANOM Prec
hist(log(Frugivores_SEAsia$Soil))


colnames(Frugivores_SEAsia)
str(Frugivores_SEAsia)
summary(Frugivores_SEAsia)

# adapt data
FruitSize <- log(Frugivores_SEAsia$FruitSize)
BodySize <- Frugivores_SEAsia$BodySize/100000
NPP_mean <- Frugivores_SEAsia$NPP_mean/10000
ensLGM_Tmean <- Frugivores_SEAsia$ensLGM_Tmean/100
ensLGM_Pmean <- Frugivores_SEAsia$ensLGM_Pmean/1000
Soil <- log(Frugivores_SEAsia$Soil)
PC1 <- Frugivores_SEAsia$PC1
PC2 <- Frugivores_SEAsia$PC2
PC3 <- Frugivores_SEAsia$PC3

Frugivores_SEAsia.dat <- data.frame(FruitSize, BodySize, NPP_mean, ensLGM_Tmean, ensLGM_Pmean, Soil, PC1, PC2, PC3)
summary(Frugivores_SEAsia.dat)
varTable(Frugivores_SEAsia.dat)

plot(Frugivores_SEAsia.dat)

# SEM approach
Frugivores_SEA_Model_1 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + ensLGM_Pmean + PC1 + PC2 + PC3
BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean
NPP_mean ~ Soil + PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean
'
sem.fit.1 <- sem(Frugivores_SEA_Model_1, data=Frugivores_SEAsia.dat)
summary(sem.fit.1, stand=T, rsq=T, fit.measures=T)

modindices(sem.fit.1)
resid(sem.fit.1, type="standardized")


# delete PC3 from the fruit size model
Frugivores_SEA_Model_2 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + ensLGM_Pmean + PC1 + PC2
BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean
NPP_mean ~ Soil + PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean
'
sem.fit.2 <- sem(Frugivores_SEA_Model_2, data=Frugivores_SEAsia.dat)
summary(sem.fit.2, stand=T, rsq=T, fit.measures=T)


# delete NPP from the body size model
Frugivores_SEA_Model_3 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + ensLGM_Pmean + PC1 + PC2
BodySize ~ PC1 + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean
NPP_mean ~ Soil + PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean
'
sem.fit.3 <- sem(Frugivores_SEA_Model_3, data=Frugivores_SEAsia.dat)
summary(sem.fit.3, stand=T, rsq=T, fit.measures=T)


# delete Soil from the NPP model
Frugivores_SEA_Model_4 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + ensLGM_Pmean + PC1 + PC2
BodySize ~ PC1 + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean
NPP_mean ~ PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean
'
sem.fit.4 <- sem(Frugivores_SEA_Model_4, data=Frugivores_SEAsia.dat)
summary(sem.fit.4, stand=T, rsq=T, fit.measures=T)


#WDK changed to
# delete PC2 from the fruit size model
Frugivores_SEA_Model_5 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + ensLGM_Pmean + PC1
BodySize ~ PC1 + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean
NPP_mean ~ PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean
'
sem.fit.5 <- sem(Frugivores_SEA_Model_5, data=Frugivores_SEAsia.dat)
summary(sem.fit.5, stand=T, rsq=T, fit.measures=T)



#WDK changed to
# delete body size from the fruit size model
Frugivores_SEA_Model_6 <- 'FruitSize ~ Soil + NPP_mean + ensLGM_Tmean + ensLGM_Pmean + PC1
BodySize ~ PC1 + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean
NPP_mean ~ PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean
'
sem.fit.6 <- sem(Frugivores_SEA_Model_6, data=Frugivores_SEAsia.dat)
summary(sem.fit.6, stand=T, rsq=T, fit.measures=T)

#WDK changed to
# delete PC1 from the fruit size model
Frugivores_SEA_Model_7 <- 'FruitSize ~ Soil + NPP_mean + ensLGM_Tmean + ensLGM_Pmean
BodySize ~ PC1 + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean
NPP_mean ~ PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean
'
sem.fit.7 <- sem(Frugivores_SEA_Model_7, data=Frugivores_SEAsia.dat)
summary(sem.fit.7, stand=T, rsq=T, fit.measures=T)

#WDK changed to
# delete ensLGM_Tmean from the fruit size model
Frugivores_SEA_Model_8 <- 'FruitSize ~ Soil + NPP_mean + ensLGM_Pmean
BodySize ~ PC1 + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean
NPP_mean ~ PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean
'
sem.fit.8 <- sem(Frugivores_SEA_Model_8, data=Frugivores_SEAsia.dat)
summary(sem.fit.8, stand=T, rsq=T, fit.measures=T)

#WDK changed to
# delete NPP from the fruit size model
Frugivores_SEA_Model_9 <- 'FruitSize ~ Soil + ensLGM_Pmean
BodySize ~ PC1 + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean
NPP_mean ~ PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean
'
sem.fit.9 <- sem(Frugivores_SEA_Model_9, data=Frugivores_SEAsia.dat)
summary(sem.fit.9, stand=T, rsq=T, fit.measures=T)

#WDK changed to
# delete PC3 from the NPP model
Frugivores_SEA_Model_11 <- 'FruitSize ~ Soil + ensLGM_Pmean
BodySize ~ PC1 + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean
NPP_mean ~ PC1 + PC2 + ensLGM_Pmean + ensLGM_Tmean
'
sem.fit.11 <- sem(Frugivores_SEA_Model_11, data=Frugivores_SEAsia.dat)
summary(sem.fit.11, stand=T, rsq=T, fit.measures=T)

#WDK changed to
# delete ensLGM_Pmean from the NPP model
Frugivores_SEA_Model_12 <- 'FruitSize ~ Soil + ensLGM_Pmean
BodySize ~ PC1 + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean
NPP_mean ~ PC1 + PC2  + ensLGM_Tmean
'
sem.fit.12 <- sem(Frugivores_SEA_Model_12, data=Frugivores_SEAsia.dat)
summary(sem.fit.12, stand=T, rsq=T, fit.measures=T)

#WDK changed to
# delete ensLGM_Pmean from the fruit size model
Frugivores_SEA_Model_13 <- 'FruitSize ~ Soil
BodySize ~ PC1 + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean
NPP_mean ~ PC1 + PC2  + ensLGM_Tmean
'
sem.fit.13 <- sem(Frugivores_SEA_Model_13, data=Frugivores_SEAsia.dat)
summary(sem.fit.13, stand=T, rsq=T, fit.measures=T)

#WDK changed to
# delete Soil from the fruit size model
Frugivores_SEA_Model_14 <- 'FruitSize
BodySize ~ PC1 + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean
NPP_mean ~ PC1 + PC2  + ensLGM_Tmean
'
sem.fit.14 <- sem(Frugivores_SEA_Model_14, data=Frugivores_SEAsia.dat)
summary(sem.fit.14, stand=T, rsq=T, fit.measures=T)

#WDK: final result has no predictor for fruit size. A bit strange


#WDK: below your previous model implementation, same result


# delete Body size from the fruit size model
#Frugivores_SEA_Model_5 <- 'FruitSize ~ Soil + NPP_mean + ensLGM_Tmean + ensLGM_Pmean + PC1 + PC2
#BodySize ~ PC1 + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean
#NPP_mean ~ PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean
#'
#sem.fit.5 <- sem(Frugivores_SEA_Model_5, data=Frugivores_SEAsia.dat)
#summary(sem.fit.5, stand=T, rsq=T, fit.measures=T)


# delete ensLGM_Pmean from the fruit size model
Frugivores_SEA_Model_6 <- 'FruitSize ~ Soil + NPP_mean + ensLGM_Tmean + PC1 + PC2
BodySize ~ PC1 + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean
NPP_mean ~ PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean
'
sem.fit.6 <- sem(Frugivores_SEA_Model_6, data=Frugivores_SEAsia.dat)
summary(sem.fit.6, stand=T, rsq=T, fit.measures=T)


# delete NPP from the fruit size model
Frugivores_SEA_Model_7 <- 'FruitSize ~ Soil + ensLGM_Tmean + PC1 + PC2
BodySize ~ PC1 + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean
NPP_mean ~ PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean
'
sem.fit.7 <- sem(Frugivores_SEA_Model_7, data=Frugivores_SEAsia.dat)
summary(sem.fit.7, stand=T, rsq=T, fit.measures=T)


# delete PC2 from the fruit size model
Frugivores_SEA_Model_8 <- 'FruitSize ~ Soil + ensLGM_Tmean + PC1
BodySize ~ PC1 + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean
NPP_mean ~ PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean
'
sem.fit.8 <- sem(Frugivores_SEA_Model_8, data=Frugivores_SEAsia.dat)
summary(sem.fit.8, stand=T, rsq=T, fit.measures=T)


# delete Soil from the fruit size model
Frugivores_SEA_Model_9 <- 'FruitSize ~ ensLGM_Tmean + PC1
BodySize ~ PC1 + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean
NPP_mean ~ PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean
'
sem.fit.9 <- sem(Frugivores_SEA_Model_9, data=Frugivores_SEAsia.dat)
summary(sem.fit.9, stand=T, rsq=T, fit.measures=T)


# delete PC1 from the fruit size model
Frugivores_SEA_Model_10 <- 'FruitSize ~ ensLGM_Tmean
BodySize ~ PC1 + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean
NPP_mean ~ PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean
'
sem.fit.10 <- sem(Frugivores_SEA_Model_10, data=Frugivores_SEAsia.dat)
summary(sem.fit.10, stand=T, rsq=T, fit.measures=T)


# delete ensLGM_Tmean from the fruit size model --> then fruit size model is empty!