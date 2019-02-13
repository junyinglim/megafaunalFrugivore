
## SEM analysis ========================
library(lavaan); library(semTools); library(semPlot)

# PCA of climate globally
currClimVar <- c("PREC_Sum", "PREC_CV", "P_drie_quart", "Tmean_mean", "T_cold_quart", "Temp_SD")
tdwg_env_subset <- tdwg_env[complete.cases(tdwg_env[currClimVar]), ]
pcaClim <- prcomp(tdwg_env_subset[currClimVar], scale. = TRUE, center = TRUE)
summary(pcaClim)    # first 3 axes explain ~93.8%
tdwg_env_subset <- cbind(tdwg_env_subset, pcaClim$x[,1:3])

tdwg_res_sem <- merge(tdwg_res, tdwg_env_subset, by = "LEVEL_3_CO", all.x = TRUE)
# Merge environmental variables 
write.csv(tdwg_res_sem, file.path(res.dir, "semAnalysis.csv"), row.names = FALSE)

# 
mod1 <- lm(log(meanFruitLengthFilled) ~ scale(log(presNat_meanBodySize)), data = subset(tdwg_res_sem, THREEREALM == "NewWorld"))
summary(mod1)

mod2 <- lm(log(meanFruitLengthFilled) ~ scale(log(curr_meanBodySize)), data = subset(tdwg_res_sem, THREEREALM == "NewWorld"))
summary(mod2)


# Normalize some variables (copied exactly what Goldel did here)
THREEREALM <- tdwg_res_sem$THREEREALM 
FruitSize <- log(tdwg_res_sem$meanFruitLength)
BodySize <- log(tdwg_res_sem$presNat_meanBodySize)
NPP_mean <- tdwg_res_sem$NPP_mean/10000
ensLGM_Tmean <- tdwg_res_sem$ensLGM_Tmean/100
ensLGM_Pmean <- log(tdwg_res_sem$ensLGM_Pmean)
Soil <- tdwg_res_sem$soilcount/10
PC1 <- tdwg_res_sem$PC1
PC2 <- tdwg_res_sem$PC2
PC3 <- tdwg_res_sem$PC3

frugivoreDataSEM <- data.frame(FruitSize, BodySize, NPP_mean, ensLGM_Tmean, ensLGM_Pmean, Soil, PC1, PC2, PC3, THREEREALM)

semMod1 <- 'FruitSize ~ BodySize + Soil + NPP_mean + ensLGM_Tmean + ensLGM_Pmean + PC1 + PC2 + PC3\n BodySize ~ NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tmean + ensLGM_Pmean\n NPP_mean ~ Soil + PC1 + PC2 + PC3 + ensLGM_Pmean + ensLGM_Tmean'
semfit.1 <- sem(semMod1, data = frugivoreDataSEM)
summary(semfit.1, stand=T, rsq=T, fit.measures=T, modindices = T)

fitUpdateSEM(semMod1, data = frugivoreDataSEM)
semPaths(semfit.1, what = "stand", residuals = FALSE, intercepts = FALSE, nCharNodes = 0, rotation = 3)

semfit.2 <- sem(semMod1, data = frugivoreDataSEM)
semfit.3 <- sem(semMod1, data = subset(frugivoreDataSEM, THREEREALM == "OWEast"))
semfit.4 <- sem(semMod1, data = subset(frugivoreDataSEM, THREEREALM == "OWWest"))
semfit.5 <- sem(semMod1, data = subset(frugivoreDataSEM, THREEREALM == "NewWorld"))
par(mfrow = c(4,1))
semPaths(semfit.5, what = "stand", residuals = FALSE, intercepts = FALSE, nCharNodes = 0, rotation = 3)

