## Analyze Data
# Generate summary statistics at the TDWG-level scale

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

tdwg_final$lgm_ens_Tano <- (tdwg_final$lgm_ens_Tmean/10 - 273.15) - tdwg_final$bio1_mean/10
tdwg_final$lgm_ens_Pano <- (tdwg_final$lgm_ens_Pmean/10 ) - tdwg_final$bio12_mean

## Subset data ========================
tdwg_final_glob <- subset(tdwg_final, REALM_LONG %in% c("Afrotropics", "Neotropics", "IndoMalay", "Australasia"))
tdwg_final_glob$logMedFS <- log(tdwg_final_glob$medianFruitLengthFilled)
tdwg_final_glob$logMax95FS <- log(tdwg_final_glob$max95FruitLengthFilled)

tdwg_final_glob$curr_logMedBS <- log(tdwg_final_glob$curr_medianBodySize)
tdwg_final_glob$curr_logMax95BS <- log(tdwg_final_glob$curr_max95BodySize)

tdwg_final_glob$pnat_logMedBS <- log(tdwg_final_glob$presNat_medianBodySize)
tdwg_final_glob$pnat_logMax95BS <- log(tdwg_final_glob$presNat_max95BodySize)

tdwg_final_glob$pnat_logMedBS_cons <- log(tdwg_final_glob$presNat_medianBodySize_cons)
tdwg_final_glob$pnat_logMax95BS_cons <- log(tdwg_final_glob$presNat_max95BodySize_cons)

tdwg_final_glob$pnat_logMedBS_supercons <- log(tdwg_final_glob$presNat_medianBodySize_supercons)
tdwg_final_glob$pnat_logMax95BS_supercons <- log(tdwg_final_glob$presNat_max95BodySize_supercons)

tdwg_final_nw <- subset(tdwg_final_glob, REALM_LONG == "Neotropics")
tdwg_final_oww <- subset(tdwg_final_glob, REALM_LONG == "Afrotropics")
tdwg_final_owe <- subset(tdwg_final_glob, REALM_LONG %in% c("Australasia","IndoMalay"))

scaleVars <- function(x, col, suffix){
  # Scales target columns and creates new columns with specified suffix
  x[paste0(col, suffix)] <- scale(x[col], scale = T, center = T)
  return(x)
}

scale.col <- c("logMedFS", "logMax95FS",
               "curr_logMedBS", "curr_logMax95BS",
               "pnat_logMedBS", "pnat_logMax95BS",
               "pnat_logMedBS_cons", "pnat_logMax95BS_cons",
               "pnat_logMedBS_supercons","pnat_logMax95BS_supercons",
               "globalPC1", "globalPC2", "globalPC3",
               "regionalPC1", "regionalPC2", "regionalPC3",
               "lgm_ens_Pano", "lgm_ens_Tano")

tdwg_final_glob <- scaleVars(x = tdwg_final_glob, col = scale.col, suffix = "_scl")

# Individual regions need to be scaled separately 
tdwg_final_nw <- scaleVars(tdwg_final_nw, col = scale.col, suffix = "_scl")
tdwg_final_oww <- scaleVars(tdwg_final_oww, col = scale.col, suffix = "_scl")
tdwg_final_owe <- scaleVars(tdwg_final_owe, col = scale.col, suffix = "_scl")

## OLS - Median body size ========================
glob_curr_medBS_mod <- lm(logMedFS_scl ~ curr_logMedBS_scl + globalPC1_scl + globalPC2_scl + globalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data =tdwg_final_glob, na.action = "na.fail")
glob_pnat_medBS_mod <- update(glob_curr_medBS_mod, ~.-curr_logMedBS_scl + pnat_logMedBS_scl)
vif(glob_curr_medBS_mod); vif(glob_pnat_medBS_mod)

nw_curr_medBS_mod <- lm(logMedFS_scl ~ curr_logMedBS_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_nw, na.action = "na.fail")
nw_pnat_medBS_mod <- update(nw_curr_medBS_mod, ~.-curr_logMedBS_scl  + pnat_logMedBS_scl)
vif(nw_curr_medBS_mod); vif(nw_pnat_medBS_mod)

oww_curr_medBS_mod <- lm(logMedFS_scl ~ curr_logMedBS_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_oww, na.action = "na.fail")
oww_pnat_medBS_mod <- update(oww_curr_medBS_mod, ~.-curr_logMedBS_scl + pnat_logMedBS_scl)
vif(oww_curr_medBS_mod); vif(oww_pnat_medBS_mod)

owe_curr_medBS_mod <- lm(logMedFS_scl ~ curr_logMedBS_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_owe, na.action = "na.fail")
owe_pnat_medBS_mod <- update(owe_curr_medBS_mod, ~.-curr_logMedBS_scl + pnat_logMedBS_scl )
vif(owe_curr_medBS_mod); vif(owe_pnat_medBS_mod)

glob_pnat_medBS_cons_mod <- update(glob_curr_medBS_mod, ~.-curr_logMedBS_scl + pnat_logMedBS_cons_scl)
nw_pnat_medBS_cons_mod <- update(nw_curr_medBS_mod, ~.-curr_logMedBS_scl + pnat_logMedBS_cons_scl)
oww_pnat_medBS_cons_mod <- update(oww_curr_medBS_mod, ~.-curr_logMedBS_scl + pnat_logMedBS_cons_scl)
owe_pnat_medBS_cons_mod <- update(owe_curr_medBS_mod, ~.-curr_logMedBS_scl + pnat_logMedBS_cons_scl)
vif(glob_pnat_medBS_cons_mod); vif(nw_pnat_medBS_cons_mod); vif(oww_pnat_medBS_cons_mod); vif(owe_pnat_medBS_cons_mod)

glob_pnat_medBS_supercons_mod <- update(glob_curr_medBS_mod, ~.-curr_logMedBS_scl + pnat_logMedBS_supercons_scl)
nw_pnat_medBS_supercons_mod <- update(nw_curr_medBS_mod, ~.-curr_logMedBS_scl + pnat_logMedBS_supercons_scl)
oww_pnat_medBS_supercons_mod <- update(oww_curr_medBS_mod, ~.-curr_logMedBS_scl + pnat_logMedBS_supercons_scl)
owe_pnat_medBS_supercons_mod <- update(owe_curr_medBS_mod, ~.-curr_logMedBS_scl + pnat_logMedBS_supercons_scl)
vif(glob_pnat_medBS_supercons_mod); vif(nw_pnat_medBS_supercons_mod); vif(oww_pnat_medBS_supercons_mod); vif(owe_pnat_medBS_supercons_mod)

medBS_ols_modlist <- list(glob_curr_medBS_mod, glob_pnat_medBS_mod, glob_pnat_medBS_cons_mod, glob_pnat_medBS_supercons_mod,
                          nw_curr_medBS_mod, nw_pnat_medBS_mod, nw_pnat_medBS_cons_mod, nw_pnat_medBS_supercons_mod,
                          oww_curr_medBS_mod, oww_pnat_medBS_mod, oww_pnat_medBS_cons_mod, oww_pnat_medBS_supercons_mod,
                          owe_curr_medBS_mod, owe_pnat_medBS_mod, owe_pnat_medBS_cons_mod, owe_pnat_medBS_supercons_mod)

medBS_ols_modavglist <- lapply(medBS_ols_modlist,
                               FUN = computeModelAvg)
medBS_ols_modavgres <- do.call("rbind", medBS_ols_modavglist)
medBS_ols_modavgres$Geographic.Scale <- rep(c("Global", "Neotropics", "Afrotropics", "Indotropics"), each = 28)
medBS_ols_modavgres$Scenario <- rep(c("Current", "Present-Natural", "Present-Natural (Conservative)", "Present-Natural (Very Conservative)"), each = 7)
medBS_ols_modavgres$Method <- "OLS"

write.csv(roundNumbers(medBS_ols_modavgres), file.path(res.dir, "medBS_ols_modavg.csv"),
          row.names = FALSE)

## OLS - Median body size (Cade 2015) ============
medBS_ols_cade_modavglist <- lapply(medBS_ols_modlist,
                                    FUN = computeModelAvg, beta = "partial.sd")

medBS_ols_cade_modavgres <- do.call("rbind", medBS_ols_cade_modavglist)
medBS_ols_cade_modavgres$Geographic.Scale <- rep(c("Global", "Neotropics", "Afrotropics", "Indotropics"), each = 28)
medBS_ols_cade_modavgres$Scenario <- rep(c("Current", "Present-Natural", "Present-Natural (Conservative)", "Present-Natural (Very Conservative)"), each = 7)
medBS_ols_cade_modavgres$Method <- "OLS"

write.csv(roundNumbers(medBS_ols_cade_modavgres), file.path(res.dir, "medBS_ols_cade_modavg.csv"),
          row.names = FALSE)

## OLS - Maximum body size ============
glob_curr_maxBS_mod <- lm(logMax95FS_scl ~ curr_logMax95BS_scl + globalPC1_scl + globalPC2_scl + globalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data =tdwg_final_glob, na.action = "na.fail")
glob_pnat_maxBS_mod <- update(glob_curr_maxBS_mod, ~.-curr_logMax95BS_scl + pnat_logMax95BS_scl)
vif(glob_curr_maxBS_mod); vif(glob_pnat_maxBS_mod)

nw_curr_maxBS_mod <- lm(logMax95FS_scl ~ curr_logMax95BS_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_nw, na.action = "na.fail")
nw_pnat_maxBS_mod <- update(nw_curr_maxBS_mod, ~.- curr_logMax95BS_scl + pnat_logMax95BS_scl)
vif(nw_curr_maxBS_mod); vif(nw_pnat_maxBS_mod)

oww_curr_maxBS_mod <- lm(logMax95FS_scl ~ curr_logMax95BS_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_oww, na.action = "na.fail")
oww_pnat_maxBS_mod <- update(oww_curr_maxBS_mod, ~.-curr_logMax95BS_scl + pnat_logMax95BS_scl)
vif(oww_curr_maxBS_mod); vif(oww_pnat_maxBS_mod)

owe_curr_maxBS_mod <- lm(logMax95FS_scl ~ curr_logMax95BS_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_owe, na.action = "na.fail")
owe_pnat_maxBS_mod <- update(owe_curr_maxBS_mod, ~.-curr_logMax95BS_scl + pnat_logMax95BS_scl)
vif(owe_curr_maxBS_mod); vif(owe_pnat_maxBS_mod)

glob_pnat_maxBS_cons_mod <- update(glob_curr_maxBS_mod, ~.-curr_logMax95BS_scl + pnat_logMax95BS_cons_scl)
nw_pnat_maxBS_cons_mod <- update(nw_curr_maxBS_mod, ~.-curr_logMax95BS_scl + pnat_logMax95BS_cons_scl)
oww_pnat_maxBS_cons_mod <- update(oww_curr_maxBS_mod, ~.-curr_logMax95BS_scl + pnat_logMax95BS_cons_scl)
owe_pnat_maxBS_cons_mod <- update(owe_curr_maxBS_mod, ~.-curr_logMax95BS_scl + pnat_logMax95BS_cons_scl)
vif(glob_pnat_maxBS_cons_mod); vif(nw_pnat_maxBS_cons_mod); vif(oww_pnat_maxBS_cons_mod); vif(owe_pnat_maxBS_cons_mod)

glob_pnat_maxBS_supercons_mod <- update(glob_curr_maxBS_mod, ~.-curr_logMax95BS_scl + pnat_logMax95BS_supercons_scl)
nw_pnat_maxBS_supercons_mod <- update(nw_curr_maxBS_mod, ~.-curr_logMax95BS_scl + pnat_logMax95BS_supercons_scl)
oww_pnat_maxBS_supercons_mod <- update(oww_curr_maxBS_mod, ~.-curr_logMax95BS_scl + pnat_logMax95BS_supercons_scl)
owe_pnat_maxBS_supercons_mod <- update(owe_curr_maxBS_mod, ~.-curr_logMax95BS_scl + pnat_logMax95BS_supercons_scl)
vif(glob_pnat_maxBS_supercons_mod); vif(nw_pnat_maxBS_supercons_mod); vif(oww_pnat_maxBS_supercons_mod); vif(owe_pnat_maxBS_supercons_mod)

maxBS_ols_modlist <- list(glob_curr_maxBS_mod, glob_pnat_maxBS_mod, glob_pnat_maxBS_cons_mod, glob_pnat_maxBS_supercons_mod,
                          nw_curr_maxBS_mod, nw_pnat_maxBS_mod, nw_pnat_maxBS_cons_mod, nw_pnat_maxBS_supercons_mod,
                          oww_curr_maxBS_mod, oww_pnat_maxBS_mod, oww_pnat_maxBS_cons_mod, oww_pnat_maxBS_supercons_mod,
                          owe_curr_maxBS_mod, owe_pnat_maxBS_mod, owe_pnat_maxBS_cons_mod, owe_pnat_maxBS_supercons_mod)
maxBS_ols_modavglist <- lapply(maxBS_ols_modlist,
                               FUN = computeModelAvg)
maxBS_ols_modavgres <- do.call("rbind", maxBS_ols_modavglist)
maxBS_ols_modavgres$Geographic.Scale <- rep(c("Global", "Neotropics", "Afrotropics", "Indotropics"), each = 28)
maxBS_ols_modavgres$Scenario <- rep(c("Current", "Present-Natural", "Present-Natural (Conservative)", "Present-Natural (Very Conservative"), each = 7)
maxBS_ols_modavgres$Method <- "OLS"

write.csv(roundNumbers(maxBS_ols_modavgres), file.path(res.dir, "maxBS_ols_modavg.csv"), row.names = FALSE)

## OLS - Maximum body size (Cade) ============
maxBS_ols_cade_modavglist <- lapply(maxBS_ols_modlist,
                                    FUN = computeModelAvg, beta = "partial.sd" )
maxBS_ols_cade_modavgres <- do.call("rbind", maxBS_ols_cade_modavglist)
maxBS_ols_cade_modavgres$Geographic.Scale <- rep(c("Global", "Neotropics", "Afrotropics", "Indotropics"), each = 28)
maxBS_ols_cade_modavgres$Scenario <- rep(c("Current", "Present-Natural", "Present-Natural (Conservative)", "Present-Natural (Very Conservative"), each = 7)
maxBS_ols_cade_modavgres$Method <- "OLS"
write.csv(roundNumbers(maxBS_ols_cade_modavgres), file.path(res.dir, "maxBS_ols_cade_modavg.csv"), row.names = FALSE)

## Generate global partial residuals ============
generatePartialResiduals <- function(mod, modcoeff){
  target_col <- c("curr_logMax95BS_scl", "pnat_logMax95BS_scl", "pnat_logMax95BS_cons_scl", "pnat_logMax95BS_supercons_scl")
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

## SPATIAL AUTOREGRESSIVE MODELLING ============================================================
# Extract centroids for each botanical countries
glob_coords <- as.matrix(tdwg_final_glob[c("LONG","LAT")])
nw_coords <- as.matrix(tdwg_final_nw[c("LONG","LAT")])
owe_coords <- as.matrix(tdwg_final_owe[c("LONG","LAT")])
oww_coords <- as.matrix(tdwg_final_oww[c("LONG","LAT")])

# Define neighbourhoods
#nb_dnear <- dnearneigh(glob_coords, longlat = TRUE, d1 = 0, d2 = 1550) # Min dist for all units to have at least one neighbour
#nb_knear <- knn2nb(knearneigh(glob_coords, k = 1), sym = TRUE) # k-nearest neighbours
nb_dlny <- tri2nb(glob_coords) 
nb_soi_glob <- graph2nb(soi.graph(nb_dlny, glob_coords)) # Sphere of influence (symmetrical)

# Define realm based neighbourhood
nb_dlny_nw <- tri2nb(nw_coords)
nb_dlny_oww <- tri2nb(oww_coords)
nb_dlny_owe <- tri2nb(owe_coords)

nb_soi_nw <- graph2nb(soi.graph(nb_dlny_nw, nw_coords))
nb_soi_oww <- graph2nb(soi.graph(nb_dlny_oww, oww_coords))
nb_soi_owe <- graph2nb(soi.graph(nb_dlny_owe, owe_coords))

listw_soi_glob <- nb2listw(nb_soi_glob, style = "W")
listw_soi_nw <- nb2listw(nb_soi_nw, style = "W")
listw_soi_oww <- nb2listw(nb_soi_oww, style = "W")
listw_soi_owe <- nb2listw(nb_soi_owe, style = "W")

nb_def <- list(nb_soi_glob, glob_coords, nb_soi_oww, oww_coords, nb_soi_owe, owe_coords, nb_soi_nw, nw_coords)
saveRDS(nb_def, file.path(res.dir, "nb_def.rds"))

# Fit models for median 
# Apparently update doesn't work on the new errorsarlm
glob_curr_medBS_sar_mod <- spatialreg::errorsarlm(logMedFS_scl ~ curr_logMedBS_scl + globalPC1_scl + globalPC2_scl + globalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, listw = listw_soi_glob, na.action = "na.fail", data = tdwg_final_glob)
glob_pnat_medBS_sar_mod <- spatialreg::errorsarlm(logMedFS_scl ~ pnat_logMedBS_scl + globalPC1_scl + globalPC2_scl + globalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, listw = listw_soi_glob, na.action = "na.fail", data = tdwg_final_glob)
glob_pnat_medBS_cons_sar_mod <- spatialreg::errorsarlm(logMedFS_scl ~ pnat_logMedBS_cons_scl + globalPC1_scl + globalPC2_scl + globalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, listw = listw_soi_glob, na.action = "na.fail", data = tdwg_final_glob)
glob_pnat_medBS_supercons_sar_mod <- spatialreg::errorsarlm(logMedFS_scl ~ pnat_logMedBS_supercons_scl + globalPC1_scl + globalPC2_scl + globalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, listw = listw_soi_glob, na.action = "na.fail", data = tdwg_final_glob)

nw_curr_medBS_sar_mod <- spatialreg::errorsarlm(logMedFS_scl ~ curr_logMedBS_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_nw, listw = listw_soi_nw, na.action = "na.fail")
nw_pnat_medBS_sar_mod <- spatialreg::errorsarlm(logMedFS_scl ~ pnat_logMedBS_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_nw, listw = listw_soi_nw, na.action = "na.fail")
nw_pnat_medBS_cons_sar_mod <- spatialreg::errorsarlm(logMedFS_scl ~ pnat_logMedBS_cons_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_nw, listw = listw_soi_nw, na.action = "na.fail")
nw_pnat_medBS_supercons_sar_mod <- spatialreg::errorsarlm(logMedFS_scl ~ pnat_logMedBS_supercons_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_nw, listw = listw_soi_nw, na.action = "na.fail")

oww_curr_medBS_sar_mod <- spatialreg::errorsarlm(logMedFS_scl ~ curr_logMedBS_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_oww, listw = listw_soi_oww, na.action = "na.fail")
oww_pnat_medBS_sar_mod <- spatialreg::errorsarlm(logMedFS_scl ~ pnat_logMedBS_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_oww, listw = listw_soi_oww, na.action = "na.fail")
oww_pnat_medBS_cons_sar_mod <- spatialreg::errorsarlm(logMedFS_scl ~ pnat_logMedBS_cons_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_oww, listw = listw_soi_oww, na.action = "na.fail")
oww_pnat_medBS_supercons_sar_mod <- spatialreg::errorsarlm(logMedFS_scl ~ pnat_logMedBS_supercons_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_oww, listw = listw_soi_oww, na.action = "na.fail")

owe_curr_medBS_sar_mod <- spatialreg::errorsarlm(logMedFS_scl ~ curr_logMedBS_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_owe, listw = listw_soi_owe, na.action = "na.fail")
owe_pnat_medBS_sar_mod <- spatialreg::errorsarlm(logMedFS_scl ~ pnat_logMedBS_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_owe, listw = listw_soi_owe, na.action = "na.fail")
owe_pnat_medBS_cons_sar_mod <- spatialreg::errorsarlm(logMedFS_scl ~ pnat_logMedBS_cons_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_owe, listw = listw_soi_owe, na.action = "na.fail")
owe_pnat_medBS_supercons_sar_mod <- spatialreg::errorsarlm(logMedFS_scl ~ pnat_logMedBS_supercons_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_owe, listw = listw_soi_owe, na.action = "na.fail")

# SAR - max body size
glob_curr_maxBS_sar_mod <- spatialreg::errorsarlm(logMax95FS_scl ~ curr_logMax95BS_scl + globalPC1_scl + globalPC2_scl + globalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_glob, listw = listw_soi_glob, na.action = "na.fail")
glob_pnat_maxBS_sar_mod <- spatialreg::errorsarlm(logMax95FS_scl ~ pnat_logMax95BS_scl + globalPC1_scl + globalPC2_scl + globalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_glob, listw = listw_soi_glob, na.action = "na.fail")
glob_pnat_maxBS_cons_sar_mod <- spatialreg::errorsarlm(logMax95FS_scl ~ pnat_logMax95BS_cons_scl + globalPC1_scl + globalPC2_scl + globalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_glob, listw = listw_soi_glob, na.action = "na.fail")
glob_pnat_maxBS_supercons_sar_mod <- spatialreg::errorsarlm(logMax95FS_scl ~ pnat_logMax95BS_supercons_scl + globalPC1_scl + globalPC2_scl + globalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_glob, listw = listw_soi_glob, na.action = "na.fail")

nw_curr_maxBS_sar_mod <- spatialreg::errorsarlm(logMax95FS_scl ~ curr_logMax95BS_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_nw, listw = listw_soi_nw, na.action = "na.fail")
nw_pnat_maxBS_sar_mod <- spatialreg::errorsarlm(logMax95FS_scl ~ pnat_logMax95BS_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_nw, listw = listw_soi_nw, na.action = "na.fail")
nw_pnat_maxBS_cons_sar_mod <- spatialreg::errorsarlm(logMax95FS_scl ~ pnat_logMax95BS_cons_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_nw, listw = listw_soi_nw, na.action = "na.fail")
nw_pnat_maxBS_supercons_sar_mod <- spatialreg::errorsarlm(logMax95FS_scl ~ pnat_logMax95BS_supercons_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_nw, listw = listw_soi_nw, na.action = "na.fail")

oww_curr_maxBS_sar_mod <- spatialreg::errorsarlm(logMax95FS_scl ~ curr_logMax95BS_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_oww, listw = listw_soi_oww, na.action = "na.fail")
oww_pnat_maxBS_sar_mod <- spatialreg::errorsarlm(logMax95FS_scl ~ pnat_logMax95BS_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_oww, listw = listw_soi_oww, na.action = "na.fail")
oww_pnat_maxBS_cons_sar_mod <- spatialreg::errorsarlm(logMax95FS_scl ~ pnat_logMax95BS_cons_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_oww, listw = listw_soi_oww, na.action = "na.fail")
oww_pnat_maxBS_supercons_sar_mod <- spatialreg::errorsarlm(logMax95FS_scl ~ pnat_logMax95BS_supercons_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_oww, listw = listw_soi_oww, na.action = "na.fail")

owe_curr_maxBS_sar_mod <- spatialreg::errorsarlm(logMax95FS_scl ~ curr_logMax95BS_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_owe, listw = listw_soi_owe, na.action = "na.fail")
owe_pnat_maxBS_sar_mod <- spatialreg::errorsarlm(logMax95FS_scl ~ pnat_logMax95BS_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_owe, listw = listw_soi_owe, na.action = "na.fail")
owe_pnat_maxBS_cons_sar_mod <- spatialreg::errorsarlm(logMax95FS_scl ~ pnat_logMax95BS_cons_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_owe, listw = listw_soi_owe, na.action = "na.fail")
owe_pnat_maxBS_supercons_sar_mod <- spatialreg::errorsarlm(logMax95FS_scl ~ pnat_logMax95BS_supercons_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_owe, listw = listw_soi_owe, na.action = "na.fail")

# Perform model averaging
medBS_sar_modlist <- list(glob_curr_medBS_sar_mod, glob_pnat_medBS_sar_mod, glob_pnat_medBS_cons_sar_mod, glob_pnat_medBS_supercons_sar_mod,
                          nw_curr_medBS_sar_mod, nw_pnat_medBS_sar_mod, nw_pnat_medBS_cons_sar_mod, nw_pnat_medBS_supercons_sar_mod,
                          oww_curr_medBS_sar_mod,oww_pnat_medBS_sar_mod, oww_pnat_medBS_cons_sar_mod, oww_pnat_medBS_supercons_sar_mod,
                          owe_curr_medBS_sar_mod,owe_pnat_medBS_sar_mod, owe_pnat_medBS_cons_sar_mod, owe_pnat_medBS_supercons_sar_mod)
medBS_sar_modelavglist <- lapply(medBS_sar_modlist,
                                 FUN = computeModelAvg)
medBS_sar_modelavgdf <- Reduce( medBS_sar_modelavglist, f = "rbind")
medBS_sar_modelavgdf$Geographic.Scale <- rep(c("Global", "Neotropics", "Afrotropics", "Indotropics"), each = 32)
medBS_sar_modelavgdf$Scenario <- rep(c("Current", "Present-Natural", "Present-Natural (Conservative)", "Present-Natural (Very Conservative)") , each = 8 )
write.csv(medBS_sar_modelavgdf,
          file.path(res.dir, "medBS_sar_modavg.csv"), row.names = FALSE)


maxBS_sar_modlist <- list(glob_curr_maxBS_sar_mod, glob_pnat_maxBS_sar_mod, glob_pnat_maxBS_cons_sar_mod, glob_pnat_maxBS_supercons_sar_mod,
                          nw_curr_maxBS_sar_mod, nw_pnat_maxBS_sar_mod, nw_pnat_maxBS_cons_sar_mod, nw_pnat_maxBS_supercons_sar_mod,
                          oww_curr_maxBS_sar_mod, oww_pnat_maxBS_sar_mod, oww_pnat_maxBS_cons_sar_mod, oww_pnat_maxBS_supercons_sar_mod,
                          owe_curr_maxBS_sar_mod, owe_pnat_maxBS_sar_mod, owe_pnat_maxBS_cons_sar_mod, owe_pnat_maxBS_supercons_sar_mod)
maxBS_sar_modelavglist <- lapply(maxBS_sar_modlist,
                                    FUN = computeModelAvg)
maxBS_sar_modelavgdf <- Reduce( maxBS_sar_modelavglist, f = "rbind")
maxBS_sar_modelavgdf$Geographic.Scale <- rep(c("Global", "Neotropics", "Afrotropics", "Indotropics"), each = 32)
maxBS_sar_modelavgdf$Scenario <- rep(c("Current", "Present-Natural", "Present-Natural (Conservative)", "Present-Natural (Very Conservative)") , each = 8 )
write.csv(maxBS_sar_modelavgdf,
          file.path(res.dir, "maxBS_sar_modavg.csv"), row.names = FALSE)

# Making sure that model fits do not have Moran I's higher than 0
moran.test(x = residuals(glob_curr_medBS_sar_mod), listw = listw_soi_glob)
moran.test(x = residuals(glob_pnat_medBS_sar_mod), listw = listw_soi_glob)
moran.test(x = residuals(glob_pnat_medBS_cons_sar_mod), listw = listw_soi_glob)
moran.test(x = residuals(glob_pnat_medBS_supercons_sar_mod), listw = listw_soi_glob)

moran.test(x = residuals(nw_curr_medBS_sar_mod), listw = listw_soi_nw)
moran.test(x = residuals(nw_pnat_medBS_sar_mod), listw = listw_soi_nw)
moran.test(x = residuals(nw_pnat_medBS_cons_sar_mod), listw = listw_soi_nw)
moran.test(x = residuals(nw_pnat_medBS_supercons_sar_mod), listw = listw_soi_nw)

moran.test(x = residuals(oww_curr_medBS_sar_mod), listw = listw_soi_oww)
moran.test(x = residuals(oww_pnat_medBS_sar_mod), listw = listw_soi_oww)
moran.test(x = residuals(oww_pnat_medBS_cons_sar_mod), listw = listw_soi_oww)
moran.test(x = residuals(oww_pnat_medBS_supercons_sar_mod), listw = listw_soi_oww)

moran.test(x = residuals(owe_curr_medBS_sar_mod), listw = listw_soi_owe)
moran.test(x = residuals(owe_pnat_medBS_sar_mod), listw = listw_soi_owe)
moran.test(x = residuals(owe_pnat_medBS_cons_sar_mod), listw = listw_soi_owe)
moran.test(x = residuals(owe_pnat_medBS_supercons_sar_mod), listw = listw_soi_owe)

moran.test(x = residuals(glob_curr_maxBS_sar_mod), listw = listw_soi_glob)
moran.test(x = residuals(glob_pnat_maxBS_sar_mod), listw = listw_soi_glob)
moran.test(x = residuals(glob_pnat_maxBS_cons_sar_mod), listw = listw_soi_glob)
moran.test(x = residuals(glob_pnat_maxBS_supercons_sar_mod), listw = listw_soi_glob)

moran.test(x = residuals(nw_curr_maxBS_sar_mod), listw = listw_soi_nw)
moran.test(x = residuals(nw_pnat_maxBS_sar_mod), listw = listw_soi_nw)
moran.test(x = residuals(nw_pnat_maxBS_cons_sar_mod), listw = listw_soi_nw)
moran.test(x = residuals(nw_pnat_maxBS_supercons_sar_mod), listw = listw_soi_nw)

moran.test(x = residuals(oww_curr_maxBS_sar_mod), listw = listw_soi_oww)
moran.test(x = residuals(oww_pnat_maxBS_sar_mod), listw = listw_soi_oww)
moran.test(x = residuals(oww_pnat_maxBS_cons_sar_mod), listw = listw_soi_oww)
moran.test(x = residuals(oww_pnat_maxBS_supercons_sar_mod), listw = listw_soi_oww)

moran.test(x = residuals(owe_curr_maxBS_sar_mod), listw = listw_soi_owe)
moran.test(x = residuals(owe_pnat_maxBS_sar_mod), listw = listw_soi_owe)
moran.test(x = residuals(owe_pnat_maxBS_cons_sar_mod), listw = listw_soi_owe)
moran.test(x = residuals(owe_pnat_maxBS_supercons_sar_mod), listw = listw_soi_owe)

## ESTIMATING DEFAUNATION IMPACT ====================

# Fit a global model with all predictor variables
glob_curr_maxBS_delta_mod <- lm(logMax95FS ~ curr_logMax95BS + globalPC1_scl + 
                            globalPC2_scl + globalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, 
                          data = tdwg_final_glob, na.action = "na.fail")
summary(glob_curr_maxBS_delta_mod)

# Create a new prediction dataframes with lower and upper bound estimate of extinction rate
predict_cons <- with(data.frame(globalPC1_scl, globalPC2_scl, globalPC3_scl, 
                             lgm_ens_Tano_scl, lgm_ens_Pano_scl,
                             curr_logMax95BS = log(futr_cons_maxBodySize) ),
                  data = tdwg_final_glob)

predict_lib <- with(data.frame(globalPC1_scl, globalPC2_scl, globalPC3_scl, 
                                lgm_ens_Tano_scl, lgm_ens_Pano_scl,
                                curr_logMax95BS = log(futr_lib_maxBodySize) ),
                     data = tdwg_final_glob)

glob_futr_delta_cons_mod <- predict(glob_curr_maxBS_delta_mod, newdata = predict_cons) # fitted FS taking into account future extinctions
glob_futr_delta_lib_mod <- predict(glob_curr_maxBS_delta_mod, newdata = predict_lib) # fitted FS taking into account future extinctions

deltaFS_lib <- exp(fitted(glob_curr_maxBS_delta_mod)) - exp(glob_futr_delta_lib_mod)
deltaFS_cons <- exp(fitted(glob_curr_maxBS_delta_mod)) - exp(glob_futr_delta_cons_mod)

mean(deltaFS_lib)
sd(deltaFS_lib)
mean(deltaFS_cons)
sd(deltaFS_cons)

deltaFSprop_lib <- round((1 - exp(glob_futr_delta_lib_mod)/exp(fitted(glob_curr_maxBS_delta_mod))) * 100, 2)
deltaFSprop_cons <- round((1 - exp(glob_futr_delta_cons_mod)/exp(fitted(glob_curr_maxBS_delta_mod))) * 100, 2)

mean(deltaFSprop_lib)
sd(deltaFSprop_lib)
mean(deltaFSprop_cons)
sd(deltaFSprop_cons)

tdwg_final_glob$changeInMaxFruitSize <- deltaFS_lib
tdwg_final_glob$changeInMaxFruitSize_cons <- deltaFS_cons
write.csv(tdwg_final_glob, file.path(res.dir, "tdwgFruitSizeChange.csv"), row.names = F)


# Difference between the fitted log median fruit size (with current body size) and the fitted log median fruit size (with future body size)
mean(tdwg_final_glob$curr_max95BodySize - tdwg_final_glob$futr_lib_maxBodySize)
mean(tdwg_final_glob$curr_max95BodySize - tdwg_final_glob$futr_cons_maxBodySize)

mean(tdwg_final_glob$curr_max95BodySize - tdwg_final_glob$futr_lib_maxBodySize)/mean(tdwg_final_glob$curr_max95BodySize)
mean(tdwg_final_glob$curr_max95BodySize - tdwg_final_glob$futr_cons_maxBodySize)/mean(tdwg_final_glob$curr_max95BodySize)

mean(tdwg_final_glob$presNat_max95BodySize_cons)
mean(tdwg_final_glob$presNat_max95BodySize) / mean(tdwg_final_glob$presNat_max95BodySize_supercons)