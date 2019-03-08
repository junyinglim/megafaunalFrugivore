## Analyze Data
# Generate summary statistics at the TDWG-level scale

## Directories ========================
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
library(ggplot2); library(sp); library(rgdal); library(viridis); library(gridExtra); library(ggrepel); library(adespatial)
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

tdwg_final_nw <- subset(tdwg_final_glob, REALM_LONG == "Neotropics")
tdwg_final_oww <- subset(tdwg_final_glob, REALM_LONG == "Afrotropics")
tdwg_final_owe <- subset(tdwg_final_glob, REALM_LONG %in% c("Australasia","IndoMalay"))

scaleVars <- function(x, col, suffix){
  # Scales target columns and creates new columns with specified suffix
  x[paste0(col, suffix)] <- scale(x[col], scale = T, center = T)
  return(x)
}
scale.col <- c("logMedFS", "logMax95FS", "dispFruitLengthFilled",
               "curr_logMedBS", "curr_logMax95BS", "curr_dispBodySize",
               "pnat_logMedBS", "pnat_logMax95BS", "presNat_dispBodySize",
               "globalPC1", "globalPC2", "globalPC3",
               "regionalPC1", "regionalPC2", "regionalPC3",
               "lgm_ens_Pano", "lgm_ens_Tano", "soilcount")

tdwg_final_glob <- scaleVars(x = tdwg_final_glob, col = scale.col, suffix = "_scl")

# Individual regions need to be scaled separately 
tdwg_final_nw <- scaleVars(tdwg_final_nw, col = scale.col, suffix = "_scl")
tdwg_final_oww <- scaleVars(tdwg_final_oww, col = scale.col, suffix = "_scl")
tdwg_final_owe <- scaleVars(tdwg_final_owe, col = scale.col, suffix = "_scl")

## Model averaging (Median Body Size) ========================
glob_curr_medBS_mod <- lm(logMedFS_scl ~ curr_logMedBS_scl + globalPC1_scl + globalPC2_scl + globalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data =tdwg_final_glob, na.action = "na.fail")
glob_pnat_medBS_mod <- update(glob_curr_medBS_mod, ~.-curr_logMedBS_scl + pnat_logMedBS_scl)
vif(glob_curr_medBS_mod); vif(glob_pnat_medBS_mod)

glob_curr_medBS_moddr <- dredge2(glob_curr_medBS_mod, data = tdwg_final_glob, model.type = "OLS")
glob_curr_medBS_modavg_summ <- model.avg2(glob_curr_medBS_moddr)
glob_pnat_medBS_moddr <- dredge2(glob_pnat_medBS_mod)
glob_pnat_medBS_modavg_summ <- model.avg2(glob_pnat_medBS_moddr)

# Just comparing coefficients using the MuMIn package
glob_curr_medBS_moddr2 <- dredge(glob_curr_medBS_mod)
glob_curr_medBS_modavg2 <- model.avg(glob_curr_medBS_moddr2)

nw_curr_medBS_mod <- lm(logMedFS_scl ~ curr_logMedBS_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_nw, na.action = "na.fail")
nw_pnat_medBS_mod <- update(nw_curr_medBS_mod, ~.-curr_logMedBS_scl  + pnat_logMedBS_scl)
vif(nw_curr_medBS_mod); vif(nw_pnat_medBS_mod)

nw_curr_medBS_moddr <- dredge2(nw_curr_medBS_mod)
nw_curr_medBS_modavg_summ <- model.avg2(nw_curr_medBS_moddr)
nw_pnat_medBS_moddr <- dredge2(nw_pnat_medBS_mod)
nw_pnat_medBS_modavg_summ <- model.avg2(nw_pnat_medBS_moddr)

# nw_curr_medBS_moddr <- dredge(nw_curr_medBS_mod)
# nw_curr_medBS_modavg <- model.avg(nw_curr_medBS_moddr)
# nw_curr_medBS_modavg_summ <- summarizeRelImportance(nw_curr_medBS_modavg)
# nw_pnat_medBS_moddr <- dredge(nw_pnat_medBS_mod)
# nw_pnat_medBS_modavg <- model.avg(nw_pnat_medBS_moddr)
# nw_pnat_medBS_modavg_summ <- summarizeRelImportance(nw_pnat_medBS_modavg)

oww_curr_medBS_mod <- lm(logMedFS_scl ~ curr_logMedBS_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_oww, na.action = "na.fail")
oww_pnat_medBS_mod <- update(oww_curr_medBS_mod, ~.-curr_logMedBS_scl + pnat_logMedBS_scl)
vif(oww_curr_medBS_mod); vif(oww_pnat_medBS_mod)

oww_curr_medBS_moddr <- dredge2(oww_curr_medBS_mod)
oww_curr_medBS_modavg_summ <- model.avg2(oww_curr_medBS_moddr)
oww_pnat_medBS_moddr <- dredge2(oww_pnat_medBS_mod)
oww_pnat_medBS_modavg_summ <- model.avg2(oww_pnat_medBS_moddr)

owe_curr_medBS_mod <- lm(logMedFS_scl ~ curr_logMedBS_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_owe, na.action = "na.fail")
owe_pnat_medBS_mod <- update(owe_curr_medBS_mod, ~.-curr_logMedBS_scl + pnat_logMedBS_scl )
vif(owe_curr_medBS_mod); vif(owe_pnat_medBS_mod)

owe_curr_medBS_moddr <- dredge2(owe_curr_medBS_mod)
owe_curr_medBS_modavg_summ <- model.avg2(owe_curr_medBS_moddr)
owe_pnat_medBS_moddr <- dredge2(owe_pnat_medBS_mod)
owe_pnat_medBS_modavg_summ <- model.avg2(owe_pnat_medBS_moddr)

# Export data
glob_curr_medBS_modavg_summ$GeographicScale <- "Global"
glob_pnat_medBS_modavg_summ$GeographicScale <- "Global"
glob_curr_medBS_modavg_summ$Scenario <- "Current"
glob_pnat_medBS_modavg_summ$Scenario <- "Present-natural"

nw_curr_medBS_modavg_summ$GeographicScale <- "Neotropics"
nw_pnat_medBS_modavg_summ$GeographicScale <- "Neotropics"
nw_curr_medBS_modavg_summ$Scenario <- "Current"
nw_pnat_medBS_modavg_summ$Scenario <- "Present-natural"

oww_curr_medBS_modavg_summ$GeographicScale <- "Afrotropics"
oww_pnat_medBS_modavg_summ$GeographicScale <- "Afrotropics"
oww_curr_medBS_modavg_summ$Scenario <- "Current"
oww_pnat_medBS_modavg_summ$Scenario <- "Present-natural"

owe_curr_medBS_modavg_summ$GeographicScale <- "Indotropics"
owe_pnat_medBS_modavg_summ$GeographicScale <- "Indotropics"
owe_curr_medBS_modavg_summ$Scenario <- "Current"
owe_pnat_medBS_modavg_summ$Scenario <- "Present-natural"

medBS_modavg_res <- Reduce("rbind", list(glob_curr_medBS_modavg_summ,
                                         glob_pnat_medBS_modavg_summ,
                                         nw_curr_medBS_modavg_summ,
                                         nw_pnat_medBS_modavg_summ,
                                         oww_curr_medBS_modavg_summ,
                                         oww_pnat_medBS_modavg_summ,
                                         owe_curr_medBS_modavg_summ,
                                         owe_pnat_medBS_modavg_summ))
write.csv(medBS_modavg_res, file.path(res.dir, "medBS_modavg_res.csv"), row.names = FALSE)

## MAXIMUM BODY SIZE ============
glob_curr_maxBS_mod <- lm(logMax95FS_scl ~ curr_logMax95BS_scl + globalPC1_scl + globalPC2_scl + globalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data =tdwg_final_glob, na.action = "na.fail")
glob_pnat_maxBS_mod <- update(glob_curr_maxBS_mod, ~.- curr_logMax95BS_scl + pnat_logMax95BS_scl)
vif(glob_curr_maxBS_mod); vif(glob_pnat_maxBS_mod)

glob_curr_maxBS_moddr <- dredge2(glob_curr_maxBS_mod)
glob_curr_maxBS_modavg_summ <- model.avg2(glob_curr_maxBS_moddr)
glob_pnat_maxBS_moddr <- dredge2(glob_pnat_maxBS_mod)
glob_pnat_maxBS_modavg_summ <- model.avg2(glob_pnat_maxBS_moddr)

nw_curr_maxBS_mod <- lm(logMax95FS_scl ~ curr_logMax95BS_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_nw, na.action = "na.fail")
nw_pnat_maxBS_mod <- update(nw_curr_maxBS_mod, ~.- curr_logMax95BS_scl + pnat_logMax95BS_scl)
vif(nw_curr_maxBS_mod); vif(nw_pnat_maxBS_mod)

nw_curr_maxBS_moddr <- dredge2(nw_curr_maxBS_mod)
nw_curr_maxBS_modavg_summ <- model.avg2(nw_curr_maxBS_moddr)
nw_pnat_maxBS_moddr <- dredge2(nw_pnat_maxBS_mod)
nw_pnat_maxBS_modavg_summ <- model.avg2(nw_pnat_maxBS_moddr)

oww_curr_maxBS_mod <- lm(logMax95FS_scl ~ curr_logMax95BS_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_oww, na.action = "na.fail")
oww_pnat_maxBS_mod <- update(oww_curr_maxBS_mod, ~.-curr_logMax95BS_scl + pnat_logMax95BS_scl)
vif(oww_curr_maxBS_mod); vif(oww_pnat_maxBS_mod)

oww_curr_maxBS_moddr <- dredge2(oww_curr_maxBS_mod)
oww_curr_maxBS_modavg_summ <- model.avg2(oww_curr_maxBS_moddr)
oww_pnat_maxBS_moddr <- dredge2(oww_pnat_maxBS_mod)
oww_pnat_maxBS_modavg_summ <- model.avg2(oww_pnat_maxBS_moddr)

owe_curr_maxBS_mod <- lm(logMax95FS ~ curr_logMax95BS_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_owe, na.action = "na.fail")
owe_pnat_maxBS_mod <- update(owe_curr_maxBS_mod, ~.-curr_logMax95BS_scl + pnat_logMax95BS_scl)
vif(owe_curr_maxBS_mod); vif(owe_pnat_maxBS_mod)

owe_curr_maxBS_moddr <- dredge2(owe_curr_maxBS_mod)
owe_curr_maxBS_modavg_summ <- model.avg2(owe_curr_maxBS_moddr)
owe_pnat_maxBS_moddr <- dredge2(owe_pnat_maxBS_mod)
owe_pnat_maxBS_modavg_summ <- model.avg2(owe_pnat_maxBS_moddr)

# Export data
glob_curr_maxBS_modavg_summ$GeographicScale <- "Global"
glob_pnat_maxBS_modavg_summ$GeographicScale <- "Global"
glob_curr_maxBS_modavg_summ$Scenario <- "Current"
glob_pnat_maxBS_modavg_summ$Scenario <- "Present-natural"

nw_curr_maxBS_modavg_summ$GeographicScale <- "Neotropics"
nw_pnat_maxBS_modavg_summ$GeographicScale <- "Neotropics"
nw_curr_maxBS_modavg_summ$Scenario <- "Current"
nw_pnat_maxBS_modavg_summ$Scenario <- "Present-natural"

oww_curr_maxBS_modavg_summ$GeographicScale <- "Afrotropics"
oww_pnat_maxBS_modavg_summ$GeographicScale <- "Afrotropics"
oww_curr_maxBS_modavg_summ$Scenario <- "Current"
oww_pnat_maxBS_modavg_summ$Scenario <- "Present-natural"

owe_curr_maxBS_modavg_summ$GeographicScale <- "Indotropics"
owe_pnat_maxBS_modavg_summ$GeographicScale <- "Indotropics"
owe_curr_maxBS_modavg_summ$Scenario <- "Current"
owe_pnat_maxBS_modavg_summ$Scenario <- "Present-natural"

maxBS_modavg_res <- Reduce("rbind", list(glob_curr_maxBS_modavg_summ,
                                         glob_pnat_maxBS_modavg_summ,
                                         nw_curr_maxBS_modavg_summ,
                                         nw_pnat_maxBS_modavg_summ,
                                         oww_curr_maxBS_modavg_summ,
                                         oww_pnat_maxBS_modavg_summ,
                                         owe_curr_maxBS_modavg_summ,
                                         owe_pnat_maxBS_modavg_summ))
write.csv(maxBS_modavg_res, file.path(res.dir, "maxBS_modavg_res.csv"), row.names = FALSE)

## DISPERSION IN BODY SIZE ============
tdwg_final_glob_disp <- subset(tdwg_final_glob, palm_nSp > 1 & curr_nSp > 1 & presNat_nSp > 1)
tdwg_final_nw_disp <- subset(tdwg_final_glob_disp, REALM_LONG == "Neotropics")
tdwg_final_oww_disp <- subset(tdwg_final_glob_disp, REALM_LONG == "Afrotropics")
tdwg_final_owe_disp <- subset(tdwg_final_glob_disp, REALM_LONG %in% c("Australasia","IndoMalay"))

glob_curr_dispBS_mod <- lm(dispFruitLengthFilled_scl ~ curr_dispBodySize_scl + globalPC1_scl + globalPC2_scl + globalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_glob_disp, na.action = "na.fail")
glob_pnat_dispBS_mod <- update(glob_curr_dispBS_mod, ~.- curr_dispBodySize_scl + presNat_dispBodySize_scl)
vif(glob_curr_dispBS_mod); vif(glob_pnat_dispBS_mod)

glob_curr_dispBS_moddr <- dredge2(glob_curr_dispBS_mod)
glob_curr_dispBS_modavg_summ <- model.avg2(glob_curr_dispBS_moddr)

glob_pnat_dispBS_moddr <- dredge2(glob_pnat_dispBS_mod)
glob_pnat_dispBS_modavg_summ <- model.avg2(glob_pnat_dispBS_moddr)

nw_curr_dispBS_mod <- lm(dispFruitLengthFilled_scl ~ curr_dispBodySize_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data =tdwg_final_nw_disp, na.action = "na.fail")
nw_pnat_dispBS_mod <- update(nw_curr_dispBS_mod,~.- curr_dispBodySize_scl + presNat_dispBodySize_scl)
vif(nw_curr_dispBS_mod); vif(nw_pnat_dispBS_mod)

nw_curr_dispBS_moddr <- dredge2(nw_curr_dispBS_mod)
nw_curr_dispBS_modavg_summ <- model.avg2(nw_curr_dispBS_moddr)
nw_pnat_dispBS_moddr <- dredge2(nw_pnat_dispBS_mod)
nw_pnat_dispBS_modavg_summ <- model.avg2(nw_pnat_dispBS_moddr)

oww_curr_dispBS_mod <- lm(dispFruitLengthFilled ~ curr_dispBodySize_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_oww_disp, na.action = "na.fail")
oww_pnat_dispBS_mod <- update(oww_curr_dispBS_mod,~.- curr_dispBodySize_scl + presNat_dispBodySize_scl)
vif(oww_curr_dispBS_mod); vif(oww_pnat_dispBS_mod)

oww_curr_dispBS_moddr <- dredge2(oww_curr_dispBS_mod)
oww_curr_dispBS_modavg_summ <- model.avg2(oww_curr_dispBS_moddr)
oww_pnat_dispBS_moddr <- dredge2(oww_pnat_dispBS_mod)
oww_pnat_dispBS_modavg_summ <- model.avg2(oww_pnat_dispBS_moddr)

owe_curr_dispBS_mod <- lm(dispFruitLengthFilled ~ curr_dispBodySize_scl + regionalPC1_scl + regionalPC2_scl + regionalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data =tdwg_final_owe_disp, na.action = "na.fail")
owe_pnat_dispBS_mod <- update(owe_curr_dispBS_mod,~.- curr_dispBodySize_scl + presNat_dispBodySize_scl)
vif(owe_curr_dispBS_mod); vif(owe_pnat_dispBS_mod)

owe_curr_dispBS_moddr <- dredge2(owe_curr_dispBS_mod)
owe_curr_dispBS_modavg_summ <- model.avg2(owe_curr_dispBS_moddr)
owe_pnat_dispBS_moddr <- dredge2(owe_pnat_dispBS_mod)
owe_pnat_dispBS_modavg_summ <- model.avg2(owe_pnat_dispBS_moddr)

# Export data
glob_curr_dispBS_modavg_summ$GeographicScale <- "Global"
glob_pnat_dispBS_modavg_summ$GeographicScale <- "Global"
glob_curr_dispBS_modavg_summ$Scenario <- "Current"
glob_pnat_dispBS_modavg_summ$Scenario <- "Present-natural"

nw_curr_dispBS_modavg_summ$GeographicScale <- "Neotropics"
nw_pnat_dispBS_modavg_summ$GeographicScale <- "Neotropics"
nw_curr_dispBS_modavg_summ$Scenario <- "Current"
nw_pnat_dispBS_modavg_summ$Scenario <- "Present-natural"

oww_curr_dispBS_modavg_summ$GeographicScale <- "Afrotropics"
oww_pnat_dispBS_modavg_summ$GeographicScale <- "Afrotropics"
oww_curr_dispBS_modavg_summ$Scenario <- "Current"
oww_pnat_dispBS_modavg_summ$Scenario <- "Present-natural"

owe_curr_dispBS_modavg_summ$GeographicScale <- "Indotropics"
owe_pnat_dispBS_modavg_summ$GeographicScale <- "Indotropics"
owe_curr_dispBS_modavg_summ$Scenario <- "Current"
owe_pnat_dispBS_modavg_summ$Scenario <- "Present-natural"

dispBS_modavg_res <- Reduce("rbind", list(glob_curr_dispBS_modavg_summ,
                                         glob_pnat_dispBS_modavg_summ,
                                         nw_curr_dispBS_modavg_summ,
                                         nw_pnat_dispBS_modavg_summ,
                                         oww_curr_dispBS_modavg_summ,
                                         oww_pnat_dispBS_modavg_summ,
                                         owe_curr_dispBS_modavg_summ,
                                         owe_pnat_dispBS_modavg_summ))
write.csv(dispBS_modavg_res, file.path(res.dir, "dispBS_modavg_res.csv"), row.names = FALSE)

## Spatial eigenvectors modelling ========
library(spdep)
tdwg_shp_raw <- readOGR("~/Dropbox/Projects/2019/palms/data/TDWG/level3/level3.shp")
glob_list <- tdwg_final_glob$LEVEL_3_CO
tdwg_shp_glob <- subset(tdwg_shp_raw, LEVEL_3_CO %in% glob_list)

# Generate neighbourhoods
# `poly2nb` not used as some polygons are not adjacent to anything else (i.e., islands)
k1_mat <- knearneigh(as.matrix(tdwg_final_glob[c("LONG","LAT")]), k = 1, longlat = TRUE)
k1_nb <- knn2nb(k1_mat)
n.comp.nb(k1_nb) # number of disjoint connected subgraphs

k2_mat <- knearneigh(as.matrix(tdwg_final_glob[c("LONG","LAT")]), k = 2, longlat = TRUE)
k2_nb <- knn2nb(k2_mat)
n.comp.nb(k2_nb) # number of disjoint connected subgraphs

k3_mat <- knearneigh(as.matrix(tdwg_final_glob[c("LONG","LAT")]), k = 3, longlat = TRUE)
k3_nb <- knn2nb(k3_mat)
n.comp.nb(k3_nb)

k4_mat <- knearneigh(as.matrix(tdwg_final_glob[c("LONG","LAT")]), k = 4, longlat = TRUE)
k4_nb <- knn2nb(k4_mat)
n.comp.nb(k4_nb)

k5_mat <- knearneigh(as.matrix(tdwg_final_glob[c("LONG","LAT")]), k = 5, longlat = TRUE)
k5_nb <- knn2nb(k5_mat)
n.comp.nb(k5_nb) # 3 subgroups are the different realms

k6_mat <- knearneigh(as.matrix(tdwg_final_glob[c("LONG","LAT")]), k = 6, longlat = TRUE)
k6_nb <- knn2nb(k6_mat)
n.comp.nb(k6_nb) # Number of subgraphs

pdf(file.path(fig.dir, "knn.pdf"), width = 12, height = 6)
par(mfrow = c(2,3))
plot(tdwg_shp_glob, col = "grey", border = "white", main = "k-nearest neighbours = 1")
plot(k1_nb, tdwg_final_glob[c("LONG","LAT")], add = T, col = "red")
plot(tdwg_shp_glob, col = "grey", border = "white", main = "k-nearest neighbours = 2")
plot(k2_nb, tdwg_final_glob[c("LONG","LAT")], add = T, col = "blue")
plot(tdwg_shp_glob, col = "grey", border = "white", main = "k-nearest neighbours = 3")
plot(k3_nb, tdwg_final_glob[c("LONG","LAT")], add = T, col = "gold")
plot(tdwg_shp_glob, col = "grey", border = "white", main = "k-nearest neighbours = 4")
plot(k4_nb, tdwg_final_glob[c("LONG","LAT")], add = T, col = "purple")
plot(tdwg_shp_glob, col = "grey", border = "white", main = "k-nearest neighbours = 5")
plot(k5_nb, tdwg_final_glob[c("LONG","LAT")], add = T, col = "green")
plot(tdwg_shp_glob, col = "grey", border = "white", main = "k-nearest neighbours = 6")
plot(k6_nb, tdwg_final_glob[c("LONG","LAT")], add = T, col = "pink")
dev.off()

# Generate spatial weighting matrices based on great circle distances
k6_nb_dist <- nbdists(k6_nb, coords = as.matrix(tdwg_final_glob[c("LONG","LAT")]), longlat = TRUE) # calculate distances between neighbours
fdist <- lapply(k6_nb_dist, function(x) 1 - x/max(dist(as.matrix(tdwg_final_glob[c("LONG","LAT")])))) # convert into spatial weights
listwgab <- nb2listw(neighbours = k6_nb, glist = fdist, style = "B") # essentially global weighting
mem.gab <- mem(listwgab) # generates Moran's eigenvector maps where correspionding eigenvalues are linearly related to Moran's index of spatial autocorrelation
barplot(attr(mem.gab, "values"), main = "Eigenvalues of spatial weighting matrix")

# Perform Moran's I on each eigenvector
round(attr(mem.gab, "values"), 3)
moranItest <- moran.randtest(x = mem.gab, 
                             list = listwgab, nrepet = 999)
signi <- which(moranItest$pvalue < 0.05) # identify eigenvectors that are positively significant (i.e., positive spatial autocorrelation)
plot(mem.gab[,signi[35:36]], SpORcoords = as.matrix(tdwg_final_glob[c("LONG","LAT")]), nb = k6_nb )

# F
tdwg_final_glob_sp <- cbind(tdwg_final_glob, mem.gab)
glob_med_f <- formula(paste0("logMedFS_scl ~ curr_logMedBS_scl + globalPC1_scl + globalPC2_scl + globalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl", "+", paste0(names(mem.gab)[signi], collapse = "+")))

glob_curr_medBS_spmod <- lm(logMedFS_scl ~ curr_logMedBS_scl + globalPC1_scl + globalPC2_scl + globalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl + MEM1, data =tdwg_final_glob_sp, na.action = "na.fail")
#glob_pnat_medBS_spmod <- update(glob_curr_medBS_spmod, ~.-curr_logMedBS_scl + pnat_logMedBS_scl)
summary(glob_curr_medBS_spmod)

# 
sp.cor <- sp.correlogram(k6_nb, tdwg_final_glob$logMedFS_scl, order=5,
                         method="I", randomisation=T)
sp.cor <- sp.correlogram(k6_nb, residuals(glob_curr_medBS_mod), order=5,
                         method="I", randomisation=T)
sp.cor <- sp.correlogram(k6_nb, residuals(glob_curr_medBS_spmod), order=5,
                          method="I", randomisation=T)
par(mfrow = c(1,1))
plot(sp.cor, ylim = c(-1,1))
abline(h = 0.1, col = "red")
# https://cran.r-project.org/web/packages/adespatial/vignettes/tutorial.html#ref-Dray2012


## Spatial autoregressive modelling ========
glob_coords <- as.matrix(tdwg_final_glob[c("LONG","LAT")])

# Distance based neighbourhood, minimum required for all units to have at least one neighbour, no distance-based spatial weighting
dnear1 <- dnearneigh(glob_coords, longlat = TRUE, d1 = 0, d2 = 1550) # no units with no links
n.comp.nb(dnear1)
dnear1_nsw <- nb2listw(dnear1, style = "W")

# Distance based neighbourhood distance weighting
fdist <- lapply(dnear1, function(x) 1 - x/max(dist(glob_coords)))
dnear1_sw <- nb2listw(neighbours = dnear1, glist = fdist, style = "W")

# k-nearest neighbours
knear1 <- knearneigh(glob_coords, k = 1)
knear1 <- knn2nb(knear1, sym = TRUE)
n.comp.nb(knear1)
knear1_nsw <- nb2listw(knear1, style = "W")

# k-nearest neighbourhood distance weighting
fdist <- lapply(knear1, function(x) { 1 - x/max(dist(glob_coords))})
knear1_sw <- nb2listw(knear1, style = "W", glist = fdist)

# Fit error SARs
dnear1_nsw_mod <- errorsarlm(logMedFS_scl ~ curr_logMedBS_scl + globalPC1_scl + globalPC2_scl + globalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_glob, listw = dnear1_nsw)
summary(dnear1_nsw_mod)

dnear1_sw_mod <- errorsarlm(logMedFS_scl ~ curr_logMedBS_scl + globalPC1_scl + globalPC2_scl + globalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_glob, listw = dnear1_sw)
summary(dnear1_sw_mod)

knear1_nsw_mod <- errorsarlm(logMedFS_scl ~ curr_logMedBS_scl + globalPC1_scl + globalPC2_scl + globalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_glob, listw = knear1_nsw)
summary(knear1_nsw_mod)

knear1_sw_mod <- errorsarlm(logMedFS_scl ~ curr_logMedBS_scl + globalPC1_scl + globalPC2_scl + globalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl, data = tdwg_final_glob, listw = knear1_sw)
summary(knear1_sw_mod)

# Calculate Moran's I to test for spatial autocorrelation in residuals
moran.mc(x = residuals(knear1_sw_mod), listw= knear1_sw, nsim = 999, alternative = "less")
moran.mc(x = residuals(knear1_sw_mod), listw= dnear1_sw, nsim = 999, alternative = 
           "greater") # test with more inclusive definition
moran.mc(x = tdwg_final_glob$logMedFS_scl, listw= knear1_nsw, nsim = 999)

