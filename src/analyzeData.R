## Analyze Data
# Fit statistical models
# Author: Jun Ying Lim
# Reference: Lim, J.Y., Svenning, J.-C., Göldel, B., Faurby, S. & Kissling, W.D. Past and future extinctions shape the body size - fruit size relationship between palms and mammalian frugivores.

## Directories ========================
rm(list = ls())
main.dir <- "/Users/junyinglim/Dropbox/Projects/2019/palms/projects/megafaunalFrugivore"
data.dir <- file.path(main.dir, "data")
res.dir <- file.path(main.dir, "results")
src.dir <- file.path(main.dir, "src")
fig.dir <- file.path(main.dir, "figs")
rawdata.dir <- "~/Dropbox/Projects/2019/palms/data"
frug.dir <- file.path(rawdata.dir, "frugivores")
options(stringsAsFactors =FALSE)

## Packages ========================
library(plyr); library(dplyr); library(reshape2)
library(ggplot2); library(sp); library(rgdal); library(viridis); library(gridExtra); library(ggrepel); library(adespatial); library(RANN); library(spatialreg)
library(car); library(MuMIn); library(stringr)
source(file.path(src.dir, "ggmodavg.R"))

## Import palm dataset ========================
tdwg_final <- read.csv(file.path(data.dir,"tdwg_final.csv"))

## Prepare data ========================
tdwg_final$logMedFS <- log(tdwg_final$medianFruitLengthFilled)
tdwg_final$logMax95FS <- log(tdwg_final$max95FruitLengthFilled)

tdwg_final$curr_logMedBS <- log(tdwg_final$curr_medianBodySize)
tdwg_final$curr_logMax95BS <- log(tdwg_final$curr_max95BodySize)

tdwg_final$pnat_logMedBS <- log(tdwg_final$presNat_medianBodySize)
tdwg_final$pnat_logMax95BS <- log(tdwg_final$presNat_max95BodySize)

tdwg_final$pnat_logMedBS_lib <- log(tdwg_final$presNat_medianBodySize_lib)
tdwg_final$pnat_logMax95BS_lib <- log(tdwg_final$presNat_max95BodySize_lib)

tdwg_final$pnat_logMedBS_cons <- log(tdwg_final$presNat_medianBodySize_cons)
tdwg_final$pnat_logMax95BS_cons <- log(tdwg_final$presNat_max95BodySize_cons)

table(tdwg_final$REALM_LONG)
tdwg_final_nw <- subset(tdwg_final, REALM_LONG == "Neotropics")
tdwg_final_oww <- subset(tdwg_final, REALM_LONG == "Afrotropics")
tdwg_final_owe <- subset(tdwg_final, REALM_LONG %in% c("Australasia","IndoMalay"))
tdwg_final_owe_noAus <- subset(tdwg_final_owe, !(REGION_NAM == "Australia" | LEVEL_NAME == "New Guinea"))

scaleVars <- function(x, col, suffix){
  # Scales target columns and creates new columns with specified suffix
  x[paste0(col, suffix)] <- scale(x[col], scale = T, center = T)
  return(x)
}

scale.col <- c("logMedFS", "logMax95FS",
               "curr_logMedBS", "curr_logMax95BS",
               "pnat_logMedBS", "pnat_logMax95BS",
               "pnat_logMedBS_lib", "pnat_logMax95BS_lib",
               "pnat_logMedBS_cons","pnat_logMax95BS_cons",
               "globalPC1", "globalPC2", "globalPC3",
               "regionalPC1", "regionalPC2", "regionalPC3",
               "regionalPC1_noAus", "regionalPC2_noAus", "regionalPC3_noAus",
               "bio12_mean", "bio15_mean", "bio17_mean", "bio1_mean", "bio11_mean", "bio4_mean",
               "lgm_ens_Pano", "lgm_ens_Tano")
tdwg_final <- scaleVars(x = tdwg_final, col = scale.col, suffix = "_scl")

# Individual regions need to be scaled separately 
tdwg_final_nw <- scaleVars(tdwg_final_nw, col = scale.col, suffix = "_scl")
tdwg_final_oww <- scaleVars(tdwg_final_oww, col = scale.col, suffix = "_scl")
tdwg_final_owe <- scaleVars(tdwg_final_owe, col = scale.col, suffix = "_scl")
tdwg_final_owe_noAus <- scaleVars(tdwg_final_owe_noAus, col = scale.col, suffix = "_scl")

## Summary statistics ========================
mean(tdwg_final$curr_max95BodySize) # 52.2 kg
mean(tdwg_final$presNat_max95BodySize) # 601.6 kg
mean(tdwg_final$presNat_max95BodySize_lib) # 630.5 kg
mean(tdwg_final$presNat_max95BodySize_cons) # 352.8 kg

mean(tdwg_final$futr_lib_max95BodySize) # 44.1 kg
mean(tdwg_final$futr_cons_max95BodySize) # 51.8 kg

## Fitting OLS models for median fruit size ========================
# Default present-natural definition
glob_curr_medBS_mod <- lm(logMedFS_scl ~ curr_logMedBS_scl + globalPC1_scl + globalPC2_scl + globalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data =tdwg_final, na.action = "na.fail")
glob_pnat_medBS_mod <- update(glob_curr_medBS_mod, ~.-curr_logMedBS_scl + pnat_logMedBS_scl)
vif(glob_curr_medBS_mod)
vif(glob_pnat_medBS_mod)

nw_curr_medBS_mod <- lm(logMedFS_scl ~ curr_logMedBS_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_nw, na.action = "na.fail")
nw_pnat_medBS_mod <- update(nw_curr_medBS_mod, ~.-curr_logMedBS_scl  + pnat_logMedBS_scl)
vif(nw_curr_medBS_mod)
vif(nw_pnat_medBS_mod)

oww_curr_medBS_mod <- lm(logMedFS_scl ~ curr_logMedBS_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_oww, na.action = "na.fail")
oww_pnat_medBS_mod <- update(oww_curr_medBS_mod, ~.-curr_logMedBS_scl + pnat_logMedBS_scl)
vif(oww_curr_medBS_mod)
vif(oww_pnat_medBS_mod)

owe_curr_medBS_mod <- lm(logMedFS_scl ~ curr_logMedBS_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_owe, na.action = "na.fail")
owe_pnat_medBS_mod <- update(owe_curr_medBS_mod, ~.-curr_logMedBS_scl + pnat_logMedBS_scl )
vif(owe_curr_medBS_mod)
vif(owe_pnat_medBS_mod)

owe_noAus_curr_medBS_mod <- lm(logMedFS_scl ~ curr_logMedBS_scl + regionalPC1_noAus_scl + regionalPC2_noAus_scl + regionalPC3_noAus_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_owe_noAus, na.action = "na.fail")
owe_noAus_pnat_medBS_mod <- update(owe_noAus_curr_medBS_mod, ~.-curr_logMedBS_scl + pnat_logMedBS_scl)
vif(owe_noAus_curr_medBS_mod)
vif(owe_noAus_pnat_medBS_mod)

# Conservative present-natural definition
glob_pnat_medBS_cons_mod <- update(glob_curr_medBS_mod, ~.-curr_logMedBS_scl + pnat_logMedBS_cons_scl)
nw_pnat_medBS_cons_mod <- update(nw_curr_medBS_mod, ~.-curr_logMedBS_scl + pnat_logMedBS_cons_scl)
oww_pnat_medBS_cons_mod <- update(oww_curr_medBS_mod, ~.-curr_logMedBS_scl + pnat_logMedBS_cons_scl)
owe_pnat_medBS_cons_mod <- update(owe_curr_medBS_mod, ~.-curr_logMedBS_scl + pnat_logMedBS_cons_scl)
owe_noAus_pnat_medBS_cons_mod <- update(owe_noAus_curr_medBS_mod, ~.-curr_logMedBS_scl + pnat_logMedBS_cons_scl)
vif(glob_pnat_medBS_cons_mod)
vif(nw_pnat_medBS_cons_mod)
vif(oww_pnat_medBS_cons_mod)
vif(owe_pnat_medBS_cons_mod)
vif(owe_noAus_pnat_medBS_cons_mod)

# Liberal present-natural definition
glob_pnat_medBS_lib_mod <- update(glob_curr_medBS_mod, ~.-curr_logMedBS_scl + pnat_logMedBS_lib_scl)
nw_pnat_medBS_lib_mod <- update(nw_curr_medBS_mod, ~.-curr_logMedBS_scl + pnat_logMedBS_lib_scl)
oww_pnat_medBS_lib_mod <- update(oww_curr_medBS_mod, ~.-curr_logMedBS_scl + pnat_logMedBS_lib_scl)
owe_pnat_medBS_lib_mod <- update(owe_curr_medBS_mod, ~.-curr_logMedBS_scl + pnat_logMedBS_lib_scl)
owe_noAus_pnat_medBS_lib_mod <- update(owe_noAus_curr_medBS_mod, ~.-curr_logMedBS_scl + pnat_logMedBS_lib_scl)
vif(glob_pnat_medBS_lib_mod)
vif(nw_pnat_medBS_lib_mod)
vif(oww_pnat_medBS_lib_mod)
vif(owe_pnat_medBS_lib_mod)
vif(owe_noAus_pnat_medBS_lib_mod)

# Implement model-averaging
medBS_ols_modlist <- list(glob_curr_medBS_mod, glob_pnat_medBS_mod, glob_pnat_medBS_cons_mod, glob_pnat_medBS_lib_mod,
                          nw_curr_medBS_mod, nw_pnat_medBS_mod, nw_pnat_medBS_cons_mod, nw_pnat_medBS_lib_mod,
                          oww_curr_medBS_mod, oww_pnat_medBS_mod, oww_pnat_medBS_cons_mod, oww_pnat_medBS_lib_mod,
                          owe_curr_medBS_mod, owe_pnat_medBS_mod, owe_pnat_medBS_cons_mod, owe_pnat_medBS_lib_mod,
                          owe_noAus_curr_medBS_mod, owe_noAus_pnat_medBS_mod, owe_noAus_pnat_medBS_cons_mod, owe_noAus_pnat_medBS_lib_mod)

medBS_ols_modavglist <- lapply(medBS_ols_modlist,
                               FUN = computeModelAvg)
medBS_ols_modavgres <- do.call("rbind", medBS_ols_modavglist)
medBS_ols_modavgres$Geographic.Scale <- rep(c("Global", "Neotropics", "Afrotropics", "Indotropics", "Indotropics (No Australia + New Guinea)"), each = 28)
medBS_ols_modavgres$Scenario <- rep(c("Current", "Present-Natural", "Present-Natural (Conservative)", "Present-Natural (Liberal)"), each = 7)
medBS_ols_modavgres$Method <- "OLS"

write.csv(roundNumbers(medBS_ols_modavgres), file.path(res.dir, "medBS_ols_modavg.csv"),
          row.names = FALSE)

# Implement-model averaging where coefficients are standardized by partial standard dev. before averaging
medBS_ols_cade_modavglist <- lapply(medBS_ols_modlist,
                                    FUN = computeModelAvg, beta = "partial.sd")

medBS_ols_cade_modavgres <- do.call("rbind", medBS_ols_cade_modavglist)
medBS_ols_cade_modavgres$Geographic.Scale <- rep(c("Global", "Neotropics", "Afrotropics", "Indotropics", "Indotropics (No Australia + New Guinea)"), each = 28)
medBS_ols_cade_modavgres$Scenario <- rep(c("Current", "Present-Natural", "Present-Natural (Conservative)", "Present-Natural (Liberal)"), each = 7)
medBS_ols_cade_modavgres$Method <- "OLS"

write.csv(roundNumbers(medBS_ols_cade_modavgres), file.path(res.dir, "medBS_ols_cade_modavg.csv"),
          row.names = FALSE)

## Fitting OLS models for maximum fruit size ========================
glob_curr_maxBS_mod <- lm(logMax95FS_scl ~ curr_logMax95BS_scl + globalPC1_scl + globalPC2_scl + globalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data =tdwg_final, na.action = "na.fail")
glob_pnat_maxBS_mod <- update(glob_curr_maxBS_mod, ~.-curr_logMax95BS_scl + pnat_logMax95BS_scl)
vif(glob_curr_maxBS_mod)
vif(glob_pnat_maxBS_mod)

nw_curr_maxBS_mod <- lm(logMax95FS_scl ~ curr_logMax95BS_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_nw, na.action = "na.fail")
nw_pnat_maxBS_mod <- update(nw_curr_maxBS_mod, ~.- curr_logMax95BS_scl + pnat_logMax95BS_scl)
vif(nw_curr_maxBS_mod)
vif(nw_pnat_maxBS_mod)

oww_curr_maxBS_mod <- lm(logMax95FS_scl ~ curr_logMax95BS_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_oww, na.action = "na.fail")
oww_pnat_maxBS_mod <- update(oww_curr_maxBS_mod, ~.-curr_logMax95BS_scl + pnat_logMax95BS_scl)
vif(oww_curr_maxBS_mod)
vif(oww_pnat_maxBS_mod)

owe_curr_maxBS_mod <- lm(logMax95FS_scl ~ curr_logMax95BS_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_owe, na.action = "na.fail")
owe_pnat_maxBS_mod <- update(owe_curr_maxBS_mod, ~.-curr_logMax95BS_scl + pnat_logMax95BS_scl)
vif(owe_curr_maxBS_mod)
vif(owe_pnat_maxBS_mod)

owe_noAus_curr_maxBS_mod <- lm(logMax95FS_scl ~ curr_logMax95BS_scl + regionalPC1_noAus_scl + regionalPC2_noAus_scl + regionalPC3_noAus_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_owe_noAus, na.action = "na.fail")
owe_noAus_pnat_maxBS_mod <- update(owe_noAus_curr_maxBS_mod, ~.-curr_logMax95BS_scl  + pnat_logMax95BS_scl)
vif(owe_noAus_curr_maxBS_mod)
vif(owe_noAus_pnat_maxBS_mod)

# Conservative present-natural definition
glob_pnat_maxBS_cons_mod <- update(glob_curr_maxBS_mod, ~.-curr_logMax95BS_scl + pnat_logMax95BS_cons_scl)
nw_pnat_maxBS_cons_mod <- update(nw_curr_maxBS_mod, ~.-curr_logMax95BS_scl + pnat_logMax95BS_cons_scl)
oww_pnat_maxBS_cons_mod <- update(oww_curr_maxBS_mod, ~.-curr_logMax95BS_scl + pnat_logMax95BS_cons_scl)
owe_pnat_maxBS_cons_mod <- update(owe_curr_maxBS_mod, ~.-curr_logMax95BS_scl + pnat_logMax95BS_cons_scl)
owe_noAus_pnat_maxBS_cons_mod <- update(owe_noAus_curr_maxBS_mod, ~.-curr_logMax95BS_scl + pnat_logMax95BS_cons_scl)
vif(glob_pnat_maxBS_cons_mod)
vif(nw_pnat_maxBS_cons_mod)
vif(oww_pnat_maxBS_cons_mod)
vif(owe_pnat_maxBS_cons_mod)
vif(owe_noAus_pnat_maxBS_cons_mod)

# Liberal present-natural definition
glob_pnat_maxBS_lib_mod <- update(glob_curr_maxBS_mod, ~.-curr_logMax95BS_scl + pnat_logMax95BS_lib_scl)
nw_pnat_maxBS_lib_mod <- update(nw_curr_maxBS_mod, ~.-curr_logMax95BS_scl + pnat_logMax95BS_lib_scl)
oww_pnat_maxBS_lib_mod <- update(oww_curr_maxBS_mod, ~.-curr_logMax95BS_scl + pnat_logMax95BS_lib_scl)
owe_pnat_maxBS_lib_mod <- update(owe_curr_maxBS_mod, ~.-curr_logMax95BS_scl + pnat_logMax95BS_lib_scl)
owe_noAus_pnat_maxBS_lib_mod <- update(owe_noAus_curr_maxBS_mod, ~.-curr_logMax95BS_scl + pnat_logMax95BS_lib_scl)
vif(glob_pnat_maxBS_lib_mod)
vif(nw_pnat_maxBS_lib_mod)
vif(oww_pnat_maxBS_lib_mod)
vif(owe_pnat_maxBS_lib_mod)
vif(owe_noAus_pnat_maxBS_lib_mod)

# Implement model-averaging
maxBS_ols_modlist <- list(glob_curr_maxBS_mod, glob_pnat_maxBS_mod, glob_pnat_maxBS_cons_mod, glob_pnat_maxBS_lib_mod,
                          nw_curr_maxBS_mod, nw_pnat_maxBS_mod, nw_pnat_maxBS_cons_mod, nw_pnat_maxBS_lib_mod,
                          oww_curr_maxBS_mod, oww_pnat_maxBS_mod, oww_pnat_maxBS_cons_mod, oww_pnat_maxBS_lib_mod,
                          owe_curr_maxBS_mod, owe_pnat_maxBS_mod, owe_pnat_maxBS_cons_mod, owe_pnat_maxBS_lib_mod,
                          owe_noAus_curr_maxBS_mod, owe_noAus_pnat_maxBS_mod, owe_noAus_pnat_maxBS_cons_mod, owe_noAus_pnat_maxBS_lib_mod)
maxBS_ols_modavglist <- lapply(maxBS_ols_modlist,
                               FUN = computeModelAvg)
maxBS_ols_modavgres <- do.call("rbind", maxBS_ols_modavglist)
maxBS_ols_modavgres$Geographic.Scale <- rep(c("Global", "Neotropics", "Afrotropics", "Indotropics", "Indotropics (No Australia + New Guinea)"), each = 28)
maxBS_ols_modavgres$Scenario <- rep(c("Current", "Present-Natural", "Present-Natural (Conservative)", "Present-Natural (Liberal)"), each = 7)
maxBS_ols_modavgres$Method <- "OLS"
write.csv(roundNumbers(maxBS_ols_modavgres), file.path(res.dir, "maxBS_ols_modavg.csv"), row.names = FALSE)

# Implement-model averaging where coefficients are standardized by partial standard dev. before averaging
maxBS_ols_cade_modavglist <- lapply(maxBS_ols_modlist,
                                    FUN = computeModelAvg, beta = "partial.sd" )
maxBS_ols_cade_modavgres <- do.call("rbind", maxBS_ols_cade_modavglist)
maxBS_ols_cade_modavgres$Geographic.Scale <- rep(c("Global", "Neotropics", "Afrotropics", "Indotropics", "Indotropics (No Australia + New Guinea)"), each = 28)
maxBS_ols_cade_modavgres$Scenario <- rep(c("Current", "Present-Natural", "Present-Natural (Conservative)", "Present-Natural (Liberal)"), each = 7)
maxBS_ols_cade_modavgres$Method <- "OLS"
write.csv(roundNumbers(maxBS_ols_cade_modavgres), file.path(res.dir, "maxBS_ols_cade_modavg.csv"), row.names = FALSE)

# Generate partial residuals 
generatePartialResiduals <- function(mod, modcoeff){
  target_col <- c("curr_logMax95BS_scl", "pnat_logMax95BS_scl", "pnat_logMax95BS_cons_scl", "pnat_logMax95BS_lib_scl")
  temp_mod <- mod
  modelavgcoeff <- setNames(modcoeff$fullAvgCoef, modcoeff$coefficient)
  # Replace model coefficients with model averaged coefficients
  temp_mod$coefficients <- modelavgcoeff[match(names(modelavgcoeff), names(temp_mod$coefficients))]
  presid_temp <- resid(temp_mod, type = "partial")
  list("points" =
         data.frame(presid = presid_temp[,colnames(presid_temp) %in% target_col],
                    x = mod$model[,names(mod$model) %in% target_col]), 
       "intercept" = modcoeff$fullAvgCoef[modcoeff$coefficient=="(Intercept)"],
       "slope" = modcoeff$fullAvgCoef[modcoeff$coefficient %in% target_col])
}

maxBS_ModPartialResid <- mapply(FUN = generatePartialResiduals, mod = maxBS_ols_modlist, modcoeff = maxBS_ols_modavglist, SIMPLIFY = F)
saveRDS(maxBS_ModPartialResid, file.path(res.dir, "maxBS_partialresid.rds"))

# Making sure the partial residuals were extracted correctly
crPlot(model = glob_curr_maxBS_mod, variable = "curr_logMax95BS_scl")
plot(presid ~ x, data = maxBS_ModPartialResid[[1]]$points)
abline(maxBS_ModPartialResid[[1]]$intercept, maxBS_ModPartialResid[[1]]$slope)

crPlot(model = glob_pnat_maxBS_mod, variable = "pnat_logMax95BS_scl")
plot(presid ~ x, data = maxBS_ModPartialResid[[2]]$points)
abline(maxBS_ModPartialResid[[2]]$intercept, maxBS_ModPartialResid[[2]]$slope)

## Defining spatial neighbourhoods ========================
# Extract centroids for each botanical countries
glob_coords <- as.matrix(tdwg_final[c("LONG","LAT")])
nw_coords <- as.matrix(tdwg_final_nw[c("LONG","LAT")])
owe_coords <- as.matrix(tdwg_final_owe[c("LONG","LAT")])
oww_coords <- as.matrix(tdwg_final_oww[c("LONG","LAT")])

# Define neighbourhoods
# nb_dnear <- dnearneigh(glob_coords, longlat = TRUE, d1 = 0, d2 = 1550) # Min dist for all units to have at least one neighbour
# nb_knear <- knn2nb(knearneigh(glob_coords, k = 1), sym = TRUE) # k-nearest neighbours
nb_dlny <- tri2nb(glob_coords) 
nb_dlny_nw <- tri2nb(nw_coords)
nb_dlny_oww <- tri2nb(oww_coords)
nb_dlny_owe <- tri2nb(owe_coords)

nb_soi_glob <- graph2nb(soi.graph(nb_dlny, glob_coords)) # Sphere of influence (symmetrical)
nb_soi_nw <- graph2nb(soi.graph(nb_dlny_nw, nw_coords))
nb_soi_oww <- graph2nb(soi.graph(nb_dlny_oww, oww_coords))
nb_soi_owe <- graph2nb(soi.graph(nb_dlny_owe, owe_coords))

listw_soi_glob <- nb2listw(nb_soi_glob, style = "W")
listw_soi_nw <- nb2listw(nb_soi_nw, style = "W")
listw_soi_oww <- nb2listw(nb_soi_oww, style = "W")
listw_soi_owe <- nb2listw(nb_soi_owe, style = "W")

nb_def <- list(nb_soi_glob, glob_coords, nb_soi_oww, oww_coords, nb_soi_owe, owe_coords, nb_soi_nw, nw_coords)
saveRDS(nb_def, file.path(res.dir, "nb_def.rds"))

## Fitting SAR models for median fruit size ========================
# Apparently update doesn't work on the new errorsarlm
glob_curr_medBS_sar_mod <- spatialreg::errorsarlm(logMedFS_scl ~ curr_logMedBS_scl + globalPC1_scl + globalPC2_scl + globalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, listw = listw_soi_glob, na.action = "na.fail", data = tdwg_final)
glob_pnat_medBS_sar_mod <- spatialreg::errorsarlm(logMedFS_scl ~ pnat_logMedBS_scl + globalPC1_scl + globalPC2_scl + globalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, listw = listw_soi_glob, na.action = "na.fail", data = tdwg_final)
glob_pnat_medBS_cons_sar_mod <- spatialreg::errorsarlm(logMedFS_scl ~ pnat_logMedBS_cons_scl + globalPC1_scl + globalPC2_scl + globalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, listw = listw_soi_glob, na.action = "na.fail", data = tdwg_final)
glob_pnat_medBS_lib_sar_mod <- spatialreg::errorsarlm(logMedFS_scl ~ pnat_logMedBS_lib_scl + globalPC1_scl + globalPC2_scl + globalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, listw = listw_soi_glob, na.action = "na.fail", data = tdwg_final)

nw_curr_medBS_sar_mod <- spatialreg::errorsarlm(logMedFS_scl ~ curr_logMedBS_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_nw, listw = listw_soi_nw, na.action = "na.fail")
nw_pnat_medBS_sar_mod <- spatialreg::errorsarlm(logMedFS_scl ~ pnat_logMedBS_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_nw, listw = listw_soi_nw, na.action = "na.fail")
nw_pnat_medBS_cons_sar_mod <- spatialreg::errorsarlm(logMedFS_scl ~ pnat_logMedBS_cons_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_nw, listw = listw_soi_nw, na.action = "na.fail")
nw_pnat_medBS_lib_sar_mod <- spatialreg::errorsarlm(logMedFS_scl ~ pnat_logMedBS_lib_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_nw, listw = listw_soi_nw, na.action = "na.fail")

oww_curr_medBS_sar_mod <- spatialreg::errorsarlm(logMedFS_scl ~ curr_logMedBS_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_oww, listw = listw_soi_oww, na.action = "na.fail")
oww_pnat_medBS_sar_mod <- spatialreg::errorsarlm(logMedFS_scl ~ pnat_logMedBS_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_oww, listw = listw_soi_oww, na.action = "na.fail")
oww_pnat_medBS_cons_sar_mod <- spatialreg::errorsarlm(logMedFS_scl ~ pnat_logMedBS_cons_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_oww, listw = listw_soi_oww, na.action = "na.fail")
oww_pnat_medBS_lib_sar_mod <- spatialreg::errorsarlm(logMedFS_scl ~ pnat_logMedBS_lib_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_oww, listw = listw_soi_oww, na.action = "na.fail")

owe_curr_medBS_sar_mod <- spatialreg::errorsarlm(logMedFS_scl ~ curr_logMedBS_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_owe, listw = listw_soi_owe, na.action = "na.fail")
owe_pnat_medBS_sar_mod <- spatialreg::errorsarlm(logMedFS_scl ~ pnat_logMedBS_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_owe, listw = listw_soi_owe, na.action = "na.fail")
owe_pnat_medBS_cons_sar_mod <- spatialreg::errorsarlm(logMedFS_scl ~ pnat_logMedBS_cons_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_owe, listw = listw_soi_owe, na.action = "na.fail")
owe_pnat_medBS_lib_sar_mod <- spatialreg::errorsarlm(logMedFS_scl ~ pnat_logMedBS_lib_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_owe, listw = listw_soi_owe, na.action = "na.fail")

# Implement model-averaging
medBS_sar_modlist <- list(glob_curr_medBS_sar_mod, glob_pnat_medBS_sar_mod, glob_pnat_medBS_cons_sar_mod, glob_pnat_medBS_lib_sar_mod,
                          nw_curr_medBS_sar_mod, nw_pnat_medBS_sar_mod, nw_pnat_medBS_cons_sar_mod, nw_pnat_medBS_lib_sar_mod,
                          oww_curr_medBS_sar_mod,oww_pnat_medBS_sar_mod, oww_pnat_medBS_cons_sar_mod, oww_pnat_medBS_lib_sar_mod,
                          owe_curr_medBS_sar_mod,owe_pnat_medBS_sar_mod, owe_pnat_medBS_cons_sar_mod, owe_pnat_medBS_lib_sar_mod)
medBS_sar_modelavglist <- lapply(medBS_sar_modlist,
                                 FUN = computeModelAvg)
medBS_sar_modelavgdf <- Reduce( medBS_sar_modelavglist, f = "rbind")
medBS_sar_modelavgdf$Geographic.Scale <- rep(c("Global", "Neotropics", "Afrotropics", "Indotropics"), each = 32)
medBS_sar_modelavgdf$Scenario <- rep(c("Current", "Present-Natural", "Present-Natural (Conservative)", "Present-Natural (Liberal)") , each = 8 )
write.csv(medBS_sar_modelavgdf,
          file.path(res.dir, "medBS_sar_modavg.csv"), row.names = FALSE)

## Fitting SAR models for maximum fruit size ========================
glob_curr_maxBS_sar_mod <- spatialreg::errorsarlm(logMax95FS_scl ~ curr_logMax95BS_scl + globalPC1_scl + globalPC2_scl + globalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final, listw = listw_soi_glob, na.action = "na.fail")
glob_pnat_maxBS_sar_mod <- spatialreg::errorsarlm(logMax95FS_scl ~ pnat_logMax95BS_scl + globalPC1_scl + globalPC2_scl + globalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final, listw = listw_soi_glob, na.action = "na.fail")
glob_pnat_maxBS_cons_sar_mod <- spatialreg::errorsarlm(logMax95FS_scl ~ pnat_logMax95BS_cons_scl + globalPC1_scl + globalPC2_scl + globalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final, listw = listw_soi_glob, na.action = "na.fail")
glob_pnat_maxBS_lib_sar_mod <- spatialreg::errorsarlm(logMax95FS_scl ~ pnat_logMax95BS_lib_scl + globalPC1_scl + globalPC2_scl + globalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final, listw = listw_soi_glob, na.action = "na.fail")

nw_curr_maxBS_sar_mod <- spatialreg::errorsarlm(logMax95FS_scl ~ curr_logMax95BS_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_nw, listw = listw_soi_nw, na.action = "na.fail")
nw_pnat_maxBS_sar_mod <- spatialreg::errorsarlm(logMax95FS_scl ~ pnat_logMax95BS_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_nw, listw = listw_soi_nw, na.action = "na.fail")
nw_pnat_maxBS_cons_sar_mod <- spatialreg::errorsarlm(logMax95FS_scl ~ pnat_logMax95BS_cons_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_nw, listw = listw_soi_nw, na.action = "na.fail")
nw_pnat_maxBS_lib_sar_mod <- spatialreg::errorsarlm(logMax95FS_scl ~ pnat_logMax95BS_lib_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_nw, listw = listw_soi_nw, na.action = "na.fail")

oww_curr_maxBS_sar_mod <- spatialreg::errorsarlm(logMax95FS_scl ~ curr_logMax95BS_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_oww, listw = listw_soi_oww, na.action = "na.fail")
oww_pnat_maxBS_sar_mod <- spatialreg::errorsarlm(logMax95FS_scl ~ pnat_logMax95BS_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_oww, listw = listw_soi_oww, na.action = "na.fail")
oww_pnat_maxBS_cons_sar_mod <- spatialreg::errorsarlm(logMax95FS_scl ~ pnat_logMax95BS_cons_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_oww, listw = listw_soi_oww, na.action = "na.fail")
oww_pnat_maxBS_lib_sar_mod <- spatialreg::errorsarlm(logMax95FS_scl ~ pnat_logMax95BS_lib_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_oww, listw = listw_soi_oww, na.action = "na.fail")

owe_curr_maxBS_sar_mod <- spatialreg::errorsarlm(logMax95FS_scl ~ curr_logMax95BS_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_owe, listw = listw_soi_owe, na.action = "na.fail")
owe_pnat_maxBS_sar_mod <- spatialreg::errorsarlm(logMax95FS_scl ~ pnat_logMax95BS_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_owe, listw = listw_soi_owe, na.action = "na.fail")
owe_pnat_maxBS_cons_sar_mod <- spatialreg::errorsarlm(logMax95FS_scl ~ pnat_logMax95BS_cons_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_owe, listw = listw_soi_owe, na.action = "na.fail")
owe_pnat_maxBS_lib_sar_mod <- spatialreg::errorsarlm(logMax95FS_scl ~ pnat_logMax95BS_lib_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_owe, listw = listw_soi_owe, na.action = "na.fail")

# Implement model-averaging
maxBS_sar_modlist <- list(glob_curr_maxBS_sar_mod, glob_pnat_maxBS_sar_mod, glob_pnat_maxBS_cons_sar_mod, glob_pnat_maxBS_lib_sar_mod,
                          nw_curr_maxBS_sar_mod, nw_pnat_maxBS_sar_mod, nw_pnat_maxBS_cons_sar_mod, nw_pnat_maxBS_lib_sar_mod,
                          oww_curr_maxBS_sar_mod, oww_pnat_maxBS_sar_mod, oww_pnat_maxBS_cons_sar_mod, oww_pnat_maxBS_lib_sar_mod,
                          owe_curr_maxBS_sar_mod, owe_pnat_maxBS_sar_mod, owe_pnat_maxBS_cons_sar_mod, owe_pnat_maxBS_lib_sar_mod)
maxBS_sar_modelavglist <- lapply(maxBS_sar_modlist,
                                    FUN = computeModelAvg)
maxBS_sar_modelavgdf <- Reduce( maxBS_sar_modelavglist, f = "rbind")
maxBS_sar_modelavgdf$Geographic.Scale <- rep(c("Global", "Neotropics", "Afrotropics", "Indotropics"), each = 32)
maxBS_sar_modelavgdf$Scenario <- rep(c("Current", "Present-Natural", "Present-Natural (Conservative)", "Present-Natural (Liberal)") , each = 8 )
write.csv(maxBS_sar_modelavgdf,
          file.path(res.dir, "maxBS_sar_modavg.csv"), row.names = FALSE)

## Moran's I tests of residuals =============================
moran.test(x = residuals(glob_curr_medBS_sar_mod), listw = listw_soi_glob)
moran.test(x = residuals(glob_pnat_medBS_sar_mod), listw = listw_soi_glob)
moran.test(x = residuals(glob_pnat_medBS_cons_sar_mod), listw = listw_soi_glob)
moran.test(x = residuals(glob_pnat_medBS_lib_sar_mod), listw = listw_soi_glob)

moran.test(x = residuals(nw_curr_medBS_sar_mod), listw = listw_soi_nw)
moran.test(x = residuals(nw_pnat_medBS_sar_mod), listw = listw_soi_nw)
moran.test(x = residuals(nw_pnat_medBS_cons_sar_mod), listw = listw_soi_nw)
moran.test(x = residuals(nw_pnat_medBS_lib_sar_mod), listw = listw_soi_nw)

moran.test(x = residuals(oww_curr_medBS_sar_mod), listw = listw_soi_oww)
moran.test(x = residuals(oww_pnat_medBS_sar_mod), listw = listw_soi_oww)
moran.test(x = residuals(oww_pnat_medBS_cons_sar_mod), listw = listw_soi_oww)
moran.test(x = residuals(oww_pnat_medBS_lib_sar_mod), listw = listw_soi_oww)

moran.test(x = residuals(owe_curr_medBS_sar_mod), listw = listw_soi_owe)
moran.test(x = residuals(owe_pnat_medBS_sar_mod), listw = listw_soi_owe)
moran.test(x = residuals(owe_pnat_medBS_cons_sar_mod), listw = listw_soi_owe)
moran.test(x = residuals(owe_pnat_medBS_lib_sar_mod), listw = listw_soi_owe)

moran.test(x = residuals(glob_curr_maxBS_sar_mod), listw = listw_soi_glob)
moran.test(x = residuals(glob_pnat_maxBS_sar_mod), listw = listw_soi_glob)
moran.test(x = residuals(glob_pnat_maxBS_cons_sar_mod), listw = listw_soi_glob)
moran.test(x = residuals(glob_pnat_maxBS_lib_sar_mod), listw = listw_soi_glob)

moran.test(x = residuals(nw_curr_maxBS_sar_mod), listw = listw_soi_nw)
moran.test(x = residuals(nw_pnat_maxBS_sar_mod), listw = listw_soi_nw)
moran.test(x = residuals(nw_pnat_maxBS_cons_sar_mod), listw = listw_soi_nw)
moran.test(x = residuals(nw_pnat_maxBS_lib_sar_mod), listw = listw_soi_nw)

moran.test(x = residuals(oww_curr_maxBS_sar_mod), listw = listw_soi_oww)
moran.test(x = residuals(oww_pnat_maxBS_sar_mod), listw = listw_soi_oww)
moran.test(x = residuals(oww_pnat_maxBS_cons_sar_mod), listw = listw_soi_oww)
moran.test(x = residuals(oww_pnat_maxBS_lib_sar_mod), listw = listw_soi_oww)

moran.test(x = residuals(owe_curr_maxBS_sar_mod), listw = listw_soi_owe)
moran.test(x = residuals(owe_pnat_maxBS_sar_mod), listw = listw_soi_owe)
moran.test(x = residuals(owe_pnat_maxBS_cons_sar_mod), listw = listw_soi_owe)
moran.test(x = residuals(owe_pnat_maxBS_lib_sar_mod), listw = listw_soi_owe)

## Estimating defaunation impact ====================
mean(tdwg_final$curr_max95BodySize) # 52.2kg
mean(tdwg_final$presNat_max95BodySize) # 601.5 kg
mean(tdwg_final$presNat_max95BodySize_lib) # 630.5 kg
mean(tdwg_final$presNat_max95BodySize_cons) # 352.8 kg

mean(tdwg_final$curr_max95BodySize) # 52.2 kg
mean(tdwg_final$futr_lib_max95BodySize) # 44.1 kg (now 43.1)
mean(tdwg_final$futr_cons_max95BodySize) # 51.7 kg (50.8 kg)

(mean(tdwg_final$curr_max95BodySize) - mean(tdwg_final$futr_lib_max95BodySize)) /1000
(mean(tdwg_final$curr_max95BodySize) - mean(tdwg_final$futr_cons_max95BodySize)) /1000

mean(tdwg_final$curr_maxBodySize) # 1,916 kg
mean(tdwg_final$futr_lib_maxBodySize) # 1,548 kg
mean(tdwg_final$futr_cons_maxBodySize) # 1,857 kg

# Percent decline in body size
(mean(tdwg_final$curr_max95BodySize - tdwg_final$futr_lib_max95BodySize)/mean(tdwg_final$curr_max95BodySize)) * 100 # 15.5 (now 17.4)
(mean(tdwg_final$curr_max95BodySize - tdwg_final$futr_cons_max95BodySize)/mean(tdwg_final$curr_max95BodySize)) * 100 # 0.8 (now 2.68)

(mean(tdwg_final$curr_maxBodySize - tdwg_final$futr_lib_maxBodySize)/mean(tdwg_final$curr_maxBodySize)) * 100
(mean(tdwg_final$curr_maxBodySize - tdwg_final$futr_cons_maxBodySize)/mean(tdwg_final$curr_maxBodySize)) * 100

# Fit a global model with all predictor variables
glob_curr_max95BS_delta_mod <- lm(logMax95FS ~ curr_logMax95BS + globalPC1_scl + 
                            globalPC2_scl + globalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, 
                          data = tdwg_final, na.action = "na.fail")
summary(glob_curr_max95BS_delta_mod)

# Create a new prediction dataframes with lower and upper bound estimate of extinction rate
predict_max95BS_lib <- with(data.frame(globalPC1_scl, globalPC2_scl, globalPC3_scl, 
                                       lgm_ens_Tano_scl, lgm_ens_Pano_scl,
                                       curr_logMax95BS = log(futr_lib_max95BodySize)),
                            data = tdwg_final)
predict_max95BS_cons <- with(data.frame(globalPC1_scl, globalPC2_scl, globalPC3_scl, 
                             lgm_ens_Tano_scl, lgm_ens_Pano_scl,
                             curr_logMax95BS = log(futr_cons_max95BodySize) ),
                     data = tdwg_final)

glob_futr_delta_max95BS_lib_mod <- predict(glob_curr_max95BS_delta_mod, newdata = predict_max95BS_lib) # fitted FS taking into account future extinctions
glob_futr_delta_max95BS_cons_mod <- predict(glob_curr_max95BS_delta_mod, newdata = predict_max95BS_cons) # fitted FS taking into account future extinctions

# Calculate absolute change in 95th percentile fruit size
deltaFS_max95BS_lib <- exp(fitted(glob_curr_max95BS_delta_mod)) - exp(glob_futr_delta_max95BS_lib_mod)
deltaFS_max95BS_cons <- exp(fitted(glob_curr_max95BS_delta_mod)) - exp(glob_futr_delta_max95BS_cons_mod)

tdwg_final$changeInMax95FruitSize <- deltaFS_max95BS_lib
tdwg_final$changeInMax95FruitSize_cons <- deltaFS_max95BS_cons

# Calculate percentage change in 95th percentile fruit size
deltaFSprop_max95BS_lib <- round((1 - exp(glob_futr_delta_max95BS_lib_mod)/exp(fitted(glob_curr_max95BS_delta_mod))) * 100, 2)
deltaFSprop_max95BS_cons <- round((1 - exp(glob_futr_delta_max95BS_cons_mod)/exp(fitted(glob_curr_max95BS_delta_mod))) * 100, 2)

tdwg_final$percentChangeInMax95FruitSize <- deltaFSprop_max95BS_lib
tdwg_final$percentChangeInMax95FruitSize_cons <- deltaFSprop_max95BS_cons

## Alternative projection using true maximum
tdwg_final$logMFS <- log(tdwg_final$maxFruitLengthFilled)
tdwg_final$logMBS <- log(tdwg_final$curr_maxBodySize)
glob_curr_maxBS_delta_mod <- lm(logMFS ~ logMBS + globalPC1_scl + 
                                  globalPC2_scl + globalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, 
                                data = tdwg_final, na.action = "na.fail")

# Create a new prediction dataframes with lower and upper bound estimate of extinction rate
predict_maxBS_lib <- with(data.frame(globalPC1_scl, globalPC2_scl, globalPC3_scl, 
                                     lgm_ens_Tano_scl, lgm_ens_Pano_scl,
                                     logMBS = log(futr_lib_maxBodySize)),
                          data = tdwg_final)
predict_maxBS_cons <- with(data.frame(globalPC1_scl, globalPC2_scl, globalPC3_scl, 
                                lgm_ens_Tano_scl, lgm_ens_Pano_scl,
                                logMBS = log(futr_cons_maxBodySize) ),
                     data = tdwg_final)

glob_futr_delta_maxBS_lib_mod <- predict(glob_curr_maxBS_delta_mod, newdata = predict_maxBS_lib) # fitted FS taking into account future extinctions
glob_futr_delta_maxBS_cons_mod <- predict(glob_curr_maxBS_delta_mod, newdata = predict_maxBS_cons) # fitted FS taking into account future extinctions

# Calculate absolute change in true maximum fruit size
deltaFS_maxBS_lib <- exp(fitted(glob_curr_maxBS_delta_mod)) - exp(glob_futr_delta_maxBS_lib_mod)
deltaFS_maxBS_cons <- exp(fitted(glob_curr_maxBS_delta_mod)) - exp(glob_futr_delta_maxBS_cons_mod)

tdwg_final$changeInMaxFruitSize <- deltaFS_maxBS_lib
tdwg_final$changeInMaxFruitSize_cons <- deltaFS_maxBS_cons

# Calculate percentage change in true maximum fruit size
deltaFSprop_maxBS_lib <- round((1 - exp(glob_futr_delta_maxBS_lib_mod)/exp(fitted(glob_curr_maxBS_delta_mod))) * 100, 2)
deltaFSprop_maxBS_cons <- round((1 - exp(glob_futr_delta_maxBS_cons_mod)/exp(fitted(glob_curr_maxBS_delta_mod))) * 100, 2)

tdwg_final$percentChangeInMaxFruitSize <- deltaFSprop_maxBS_lib
tdwg_final$percentChangeInMaxFruitSize_cons <- deltaFSprop_maxBS_cons

# Export fruit size change results
write.csv(tdwg_final, file.path(res.dir, "tdwgFruitSizeChange.csv"), row.names = F)

# Fruit size change summary statistics
tdwg_final <- read.csv(file.path(res.dir, "tdwgFruitSizeChange.csv"))
mean(tdwg_final$changeInMax95FruitSize) # 0.22 cm (now 0.25 cm)
sd(tdwg_final$changeInMax95FruitSize) # 0.42 (now 0.437)
mean(tdwg_final$changeInMax95FruitSize_cons) # 0.01 cm (now 0.03)
sd(tdwg_final$changeInMax95FruitSize_cons) # 0.08 (now 0.05)

mean(tdwg_final$percentChangeInMax95FruitSize) # 3.81 % (now 4.1%)
sd(tdwg_final$percentChangeInMax95FruitSize) # 8.25
mean(tdwg_final$percentChangeInMax95FruitSize_cons) # 0.17 % (now 0.45%)
sd(tdwg_final$percentChangeInMax95FruitSize_cons) # 1.51

mean(tdwg_final$changeInMaxFruitSize) # 0.32 cm (now 0.33 cm)
sd(tdwg_final$changeInMaxFruitSize) # 0.47
mean(tdwg_final$changeInMaxFruitSize_cons) # 0.03 cm (0.04 cm)
sd(tdwg_final$changeInMaxFruitSize_cons) # 0.03

mean(tdwg_final$percentChangeInMaxFruitSize) # 3.71 %
sd(tdwg_final$percentChangeInMaxFruitSize) # 6.63
mean(tdwg_final$percentChangeInMaxFruitSize_cons) # 0.36 %
sd(tdwg_final$percentChangeInMaxFruitSize_cons) # 0.35


mean(tdwg_final$max95FruitLengthFilled)
median(tdwg_final$max95FruitLengthFilled)
dim(tdwg_final)
