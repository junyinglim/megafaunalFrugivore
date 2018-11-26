## 
# Created by: Bastian Golden, Modified by: Jun Ying Lim

# Bastian Goldel
# Project: Global Palm and Mammals Functional traits distributions and relation  (category fruits = all frugivoreData)
# added: extinct frugivoreData & replaced current ranges of present frugivoreData by natural ranges
######## CATEGORY FRUIT = ALL (also extinct and natural) frugivoreData #########

## DIRECTORIES ===============
main.dir <- "~/Dropbox/Projects/2019/palms/projects/megafaunalFrugivore/"
data.dir <- file.path(main.dir,"data")
src.dir <- file.path(main.dir, "src")
fig.dir <- file.path(main.dir, "figs")

## PACKAGES ===============
library(lavaan); library(semPlot)
source(file.path(src.dir, "semTools.R"))

## IMPORT DATA ===============
# Includes all mammals, current and present-natural
frugivoreData <- read.table(file.path(data.dir, "ALL_MAMMALS_COMPLETELY_230915.txt"), header = T, sep = "\t")

str(frugivoreData)
names(frugivoreData)
frugivoreDataSubset <- na.omit(frugivoreData) # ignores NAs
setdiff(frugivoreData$LEVEL_NAME,frugivoreDataSubset$LEVEL_NAME) # 2 levels are removed
subset(frugivoreData, LEVEL_NAME %in% c("Marquesas", "Tuamotu"))


## MODEL A: All biogeographic regions  ===============
# Principal component analysis
currClimVar <- c("PREC_Sum", "PREC_CV", "P_drie_quart", "Tmean_mean", "T_cold_quart", "Temp_SD")
currClim <- frugivoreDataSubset[currClimVar]
pcaClim <-prcomp(currClim, scale=TRUE)
summary(pcaClim)    # first 3 axes explain ~93.8%
#biplot(pcaClim, cex=0.9, xlim=c(-0.23,0.23))

frugivoreDataFinal<- cbind(frugivoreDataSubset, pcaClim$x[,1:3])

# Histograms
hist(log(frugivoreDataFinal$FruitSize))         # palm fruit size
hist(log(frugivoreDataFinal$BodySize))          # mammal body size
hist(frugivoreDataFinal$NPP_mean)               # NPP (net prim. production)
hist(frugivoreDataFinal$ensLGM_Tmean)           # mean LGM ANOM Temp
hist(log(frugivoreDataFinal$ensLGM_Pmean))      # mean LGM ANOM Prec
hist(frugivoreDataFinal$Soil)                   # % of sandy soil (depth 45cm)
hist(frugivoreDataFinal$PC1)
hist(frugivoreDataFinal$PC2)
hist(frugivoreDataFinal$PC3)

# plot single relationships to check for non-linear response
plot(log(frugivoreDataFinal$FruitSize) ~ log(frugivoreDataFinal$BodySize))
plot(log(frugivoreDataFinal$FruitSize) ~ frugivoreDataFinal$PC1)
plot(log(frugivoreDataFinal$FruitSize) ~ frugivoreDataFinal$PC2)
plot(log(frugivoreDataFinal$FruitSize) ~ frugivoreDataFinal$PC3)
plot(log(frugivoreDataFinal$FruitSize) ~ frugivoreDataFinal$NPP_mean)
plot(log(frugivoreDataFinal$FruitSize) ~ frugivoreDataFinal$ensLGM_Tmean)
plot(log(frugivoreDataFinal$FruitSize) ~ log(frugivoreDataFinal$ensLGM_Pmean))
plot(log(frugivoreDataFinal$FruitSize) ~ frugivoreDataFinal$Soil)

plot(log(frugivoreDataFinal$BodySize) ~ frugivoreDataFinal$PC1)
plot(log(frugivoreDataFinal$BodySize) ~ frugivoreDataFinal$PC2)
plot(log(frugivoreDataFinal$BodySize) ~ frugivoreDataFinal$PC3)
plot(log(frugivoreDataFinal$BodySize) ~ frugivoreDataFinal$NPP_mean)
plot(log(frugivoreDataFinal$BodySize) ~ frugivoreDataFinal$ensLGM_Tmean)
plot(log(frugivoreDataFinal$BodySize) ~ log(frugivoreDataFinal$ensLGM_Pmean))
plot(log(frugivoreDataFinal$BodySize) ~ frugivoreDataFinal$Soil)

# SEM models
# Scale variables (recode vars to roughly same scale (recommended!)
#frugivoreDataSEM <- frugivoreDataFinal[c("FruitSize", "BodySize", "NPP_mean", "ensLGM_Tmean", "ensLGM_Pmean","PC1", "PC2", "PC3")]

FruitSize <- log(frugivoreDataFinal$FruitSize)
BodySize <- log(frugivoreDataFinal$BodySize)
NPP_mean <- frugivoreDataFinal$NPP_mean/10000
ensLGM_Tmean <- frugivoreDataFinal$ensLGM_Tmean/100
ensLGM_Pmean <- log(frugivoreDataFinal$ensLGM_Pmean)
Soil <- frugivoreDataFinal$Soil/10
PC1 <- frugivoreDataFinal$PC1
PC2 <- frugivoreDataFinal$PC2
PC3 <- frugivoreDataFinal$PC3

frugivoreDataSEM <- data.frame(FruitSize, BodySize, NPP_mean, ensLGM_Tmean, ensLGM_Pmean, Soil, PC1, PC2, PC3)

frugivoreDataSEM$group <- as.vector(frugivoreDataFinal$THREEREALM)

varTable(frugivoreDataSEM)

# Criteria for fit measures:
# Model Chi-Square with its df and p-value: prefer p-value greater than 0.05
# Root Mean Square Error of Approximation (RMSEA): prefer lower 90%CI to be < 0.05
# Comparative Fit Index (CFI): prefer value greater than 0.90
# Mod. indexes should be << 3.84; if higher then may suggest a missing path


# JYL notes:
# Not sure what to look out for in resid() or modindices
#modindices(sem.fit.1)
#resid(sem.fit.1, type="standardized")

semMod1 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + ensLGM_Pmean + PC1 + PC2 + PC3\n BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean\n NPP_mean ~ Soil + PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean'

sem.fit.full <- sem(semMod1, data = frugivoreDataSEM)
sem.fit.1 <- sem(semMod1, data=frugivoreDataSEM, group = "group")
summary(sem.fit.full, stand=T, rsq=T, fit.measures=T, modindices = T)

semPaths(sem.fit.full, whatLabels = "est", what = "est", intercepts = FALSE, residuals = FALSE, exoCov = FALSE, rotation =3)
semPaths(sem.fit.1, whatLabels = "est", what = "est", intercepts = FALSE, residuals = FALSE, exoCov = FALSE, rotation =3, panelGroups = TRUE, title = T)
?semPaths
semPaths(sem.fit.1, whatLabels = "est", what = "est", intercepts = FALSE, residuals = FALSE, exoCov = FALSE, rotation =3, mar = c(3,1,5,1))



semModGroupEqual <- sem(semMod1, data=frugivoreDataSEM, group = "group")

semModNW <- sem(semMod1, data=subset(frugivoreDataSEM, group == "NewWorld"))
semModOWEast <- sem(semMod1, data=subset(frugivoreDataSEM, group == "OWEast"))
semModOWWest <- sem(semMod1, data=subset(frugivoreDataSEM, group == "OWWest"))

semPaths(semModNW, whatLabels = "est", what = "est", intercepts = FALSE, residuals = FALSE, exoCov = FALSE, rotation =3)
semPaths(semModOWEast, whatLabels = "est", what = "est", intercepts = FALSE, residuals = FALSE, exoCov = FALSE, rotation =3)
semPaths(semModOWWest, whatLabels = "est", what = "est", intercepts = FALSE, residuals = FALSE, exoCov = FALSE, rotation =3)

summary(sem.fit.full, stand=T, rsq=T, fit.measures=T, modindices = T)






semMod2 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + ensLGM_Pmean + PC1 + PC2 + PC3\n BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Pmean\n NPP_mean ~ Soil + PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean' # Remove BodySIze ~ ensLGM_Tmean
sem.fit.2 <- sem(semMod2, data=frugivoreDataSEM); summary(sem.fit.2, stand=T, rsq=T, fit.measures=T, modindices = T)
lavTestLRT(sem.fit.1, sem.fit.2)

semMod3 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + ensLGM_Pmean + PC1 + PC3\n BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Pmean\n NPP_mean ~ Soil + PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean' # Remove FruitSize ~ PC2
sem.fit.3 <- sem(semMod3, data=frugivoreDataSEM); summary(sem.fit.3, stand=T, rsq=T, fit.measures=T, modindices = TRUE)
lavTestLRT(sem.fit.2, sem.fit.3)

semMod4 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + ensLGM_Pmean + PC1 + PC3\n BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Pmean\n NPP_mean ~ PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean' # Remove NPP_mean ~ Soil
sem.fit.4 <- sem(semMod4, data=frugivoreDataSEM); summary(sem.fit.4, stand=T, rsq=T, fit.measures=T, modindices = TRUE)
lavTestLRT(sem.fit.3, sem.fit.4)

semMod5 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + ensLGM_Pmean + PC1\n BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Pmean\n NPP_mean ~ PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean' # Remove FruitSize ~ PC3
sem.fit.5 <- sem(semMod5, data=frugivoreDataSEM); summary(sem.fit.5, stand=T, rsq=T, fit.measures=T, modindices=T)
lavTestLRT(sem.fit.4, sem.fit.5)

semMod6 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + ensLGM_Pmean + PC1\n BodySize ~  PC1 + PC2 + PC3 + ensLGM_Pmean\n NPP_mean ~ PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean' # Remove BidySize ~ NPP_mean
sem.fit.6 <- sem(semMod6, data=frugivoreDataSEM); summary(sem.fit.6, stand=T, rsq=T, fit.measures=T, modindices=T)
lavTestLRT(sem.fit.5, sem.fit.6)

semMod7 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + ensLGM_Pmean + PC1\n BodySize ~  PC2 + PC3 + ensLGM_Pmean\n NPP_mean ~ PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean' # Remove BodySize ~ PC1
sem.fit.7 <- sem(semMod7, data=frugivoreDataSEM); summary(sem.fit.7, stand=T, rsq=T, fit.measures=T, modindices = T)
lavTestLRT(sem.fit.1, sem.fit.7)

semMod8 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + ensLGM_Pmean + PC1\n BodySize ~  PC2 + PC3\n NPP_mean ~ PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean' # Remove BodySize ~ ensLGM_Pmean
sem.fit.8 <- sem(semMod8, data=frugivoreDataSEM); summary(sem.fit.8, stand=T, rsq=T, fit.measures=T, modindices=T)
lavTestLRT(sem.fit.1, sem.fit.8)

semMod9 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + ensLGM_Pmean\n BodySize ~  PC2 + PC3\n NPP_mean ~ PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean' # Remove Fruitsize ~ PC1
sem.fit.9 <- sem(semMod9, data=frugivoreDataSEM); summary(sem.fit.9, stand=T, rsq=T, fit.measures=T)

semMod10 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean\n BodySize ~  PC2 + PC3\n NPP_mean ~ PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean' # Remove Fruitsize ~ ensLGM_Pmean
sem.fit.10 <- sem(semMod10, data=frugivoreDataSEM); summary(sem.fit.10, stand=T, rsq=T, fit.measures=T)

semMod11 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean\n BodySize ~  PC2 + PC3\n NPP_mean ~ PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean' # Remove NPP_mean ~ PC1
sem.fit.11 <- sem(semMod11, data=frugivoreDataSEM); summary(sem.fit.11, stand=T, rsq=T, fit.measures=T)
modindices(sem.fit.11)

pdf(file.path(fig.dir, "globalSEM.pdf"))
semPaths(sem.fit.11, what = "est", residuals = FALSE, intercepts = FALSE, nCharNodes = 0, rotation = 3)
dev.off()

## MODEL B: New World only  ===============
# Principal component analysis (don't do a separate PCA)
frugivoreData_NW <- subset(frugivoreDataSubset, THREEREALM == "NewWorld") # PCA of new world only

currClimVar <- c("PREC_Sum", "PREC_CV", "P_drie_quart", "Tmean_mean", "T_cold_quart", "Temp_SD")
currClimNW <- frugivoreData_NW[currClimVar]
pcaClim <-prcomp(currClimNW, scale=TRUE)
summary(pcaClim)    # first 3 axes explain ~90.3%

frugivoreDataFinal_NW<- cbind(frugivoreData_NW, pcaClim$x[,1:3])

# Rescale variables 
FruitSize <- log(frugivoreDataFinal_NW$FruitSize)
BodySize <- log(frugivoreDataFinal_NW$BodySize)
NPP_mean <- frugivoreDataFinal_NW$NPP_mean/10000
ensLGM_Tmean <- frugivoreDataFinal_NW$ensLGM_Tmean/100
ensLGM_Pmean <- log(frugivoreDataFinal_NW$ensLGM_Pmean)
Soil <- frugivoreDataFinal_NW$Soil/10
PC1 <- frugivoreDataFinal_NW$PC1
PC2 <- frugivoreDataFinal_NW$PC2
PC3 <- frugivoreDataFinal_NW$PC3

frugivoreDataSEM_NW <- data.frame(FruitSize, BodySize, NPP_mean, ensLGM_Tmean, ensLGM_Pmean, Soil, PC1, PC2, PC3)

# Structural equation models
semNWMod_1 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + ensLGM_Pmean + PC1 + PC2 + PC3\n BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean\n NPP_mean ~ Soil + PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean'
semNW.fit.1 <- sem(semNWMod_1, data=frugivoreDataSEM_NW); summary(semNW.fit.1, stand=T, rsq=T, fit.measures=T)

semNWMod_2 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + ensLGM_Pmean + PC2 + PC3\n BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean\n NPP_mean ~ Soil + PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean' # FruitSize ~ PC1
semNW.fit.2 <- sem(semNWMod_2, data=frugivoreDataSEM_NW); summary(semNW.fit.2, stand=T, rsq=T, fit.measures=T)

semNWMod_3 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + ensLGM_Pmean + PC2 + PC3\n BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean\n NPP_mean ~ Soil + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean' # NPP_mean ~ PC1
semNW.fit.3 <- sem(semNWMod_3, data=frugivoreDataSEM_NW); summary(semNW.fit.3, stand=T, rsq=T, fit.measures=T)

semNWMod_4 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + ensLGM_Pmean + PC2 + PC3\n BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean\n NPP_mean ~ Soil + PC2 + PC3 + ensLGM_Tmean' # NPP_mean ~ ensLGM_Pmean
semNW.fit.4 <- sem(semNWMod_4, data=frugivoreDataSEM_NW); summary(semNW.fit.4, stand=T, rsq=T, fit.measures=T)

semNWMod_5 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + ensLGM_Pmean + PC2 + PC3\n BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean\n NPP_mean ~ Soil + PC2 + ensLGM_Tmean' # NPP_mean ~ PC3
semNW.fit.5 <- sem(semNWMod_5, data=frugivoreDataSEM_NW); summary(semNW.fit.5, stand=T, rsq=T, fit.measures=T)

semNWMod_6 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + ensLGM_Pmean + PC2 + PC3\n BodySize ~ NPP_mean + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean\n NPP_mean ~ Soil + PC2 + ensLGM_Tmean' # BodySize ~ PC1
semNW.fit.6 <- sem(semNWMod_6, data=frugivoreDataSEM_NW); summary(semNW.fit.6, stand=T, rsq=T, fit.measures=T)

semNWMod_7 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Pmean + PC2 + PC3\n BodySize ~ NPP_mean + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean\n NPP_mean ~ Soil + PC2 + ensLGM_Tmean' # FruitSize ~ ensLGM_Tmean
semNW.fit.7 <- sem(semNWMod_7, data=frugivoreDataSEM_NW); summary(semNW.fit.7, stand=T, rsq=T, fit.measures=T)

semNWMod_8 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Pmean + PC2\n BodySize ~ NPP_mean + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean\n NPP_mean ~ Soil + PC2 + ensLGM_Tmean' # FruitSize ~ PC3
semNW.fit.8 <- sem(semNWMod_8, data=frugivoreDataSEM_NW); summary(semNW.fit.8, stand=T, rsq=T, fit.measures=T)

semNWMod_9 <- 'FruitSize ~ BodySize + Soil + NPP_mean + PC2\n BodySize ~ NPP_mean + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean\n NPP_mean ~ Soil + PC2 + ensLGM_Tmean' # FruitSize ~ ensLGM_Pmean
semNW.fit.9 <- sem(semNWMod_9, data=frugivoreDataSEM_NW); summary(semNW.fit.9, stand=T, rsq=T, fit.measures=T)

semNWMod_10 <- 'FruitSize ~ BodySize + NPP_mean + PC2\n BodySize ~ NPP_mean + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean\n NPP_mean ~ Soil + PC2 + ensLGM_Tmean' # FruitSize ~ Soil
semNW.fit.10 <- sem(semNWMod_10, data=frugivoreDataSEM_NW); summary(semNW.fit.10, stand=T, rsq=T, fit.measures=T)

semNWMod_11 <- 'FruitSize ~ BodySize + NPP_mean\n BodySize ~ NPP_mean + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean\n NPP_mean ~ Soil + PC2 + ensLGM_Tmean' # FruitSize ~ PC2
semNW.fit.11 <- sem(semNWMod_11, data=frugivoreDataSEM_NW); summary(semNW.fit.11, stand=T, rsq=T, fit.measures=T)

semNWMod_12 <- 'FruitSize ~ BodySize + NPP_mean\n BodySize ~ PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean\n NPP_mean ~ Soil + PC2 + ensLGM_Tmean' # BodySize ~ NPP_mean
semNW.fit.12 <- sem(semNWMod_12, data=frugivoreDataSEM_NW); summary(semNW.fit.12, stand=T, rsq=T, fit.measures=T)

semNWMod_13 <- 'FruitSize ~ BodySize + NPP_mean\n BodySize ~ PC3 + ensLGM_Tmean + ensLGM_Pmean\n NPP_mean ~ Soil + PC2 + ensLGM_Tmean' # BodySize ~ PC2
semNW.fit.13 <- sem(semNWMod_13, data=frugivoreDataSEM_NW); summary(semNW.fit.13, stand=T, rsq=T, fit.measures=T)

semNWMod_14 <- 'FruitSize ~ BodySize + NPP_mean\n BodySize ~ PC3 + ensLGM_Tmean\n NPP_mean ~ Soil + PC2 + ensLGM_Tmean' # BodySize ~ ensLGM_Pmean
semNW.fit.14 <- sem(semNWMod_14, data=frugivoreDataSEM_NW); summary(semNW.fit.14, stand=T, rsq=T, fit.measures=T)

semNWMod_15 <- 'FruitSize ~ BodySize + NPP_mean\n BodySize ~ ensLGM_Tmean\n NPP_mean ~ Soil + PC2 + ensLGM_Tmean' # BodySize ~ PC3
semNW.fit.15 <- sem(semNWMod_15, data=frugivoreDataSEM_NW); summary(semNW.fit.15, stand=T, rsq=T, fit.measures=T)


pdf(file.path(fig.dir, "newworldSEM.pdf"))
semPaths(semNW.fit.15, what = "est", residuals = FALSE, intercepts = FALSE, nCharNodes = 0, rotation = 3)
dev.off()


