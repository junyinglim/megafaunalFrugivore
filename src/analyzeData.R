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
library(ggplot2); library(sp); library(rgdal); library(viridis); library(gridExtra); library(ggrepel)
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

glob_curr_medBS_moddr <- dredge2(glob_curr_medBS_mod)
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

# 
# ggplot(aes(y = dispFruitLengthFilled, x=  palm_nSp), data = tdwg_final_glob) + geom_point()+geom_smooth(method = "lm")
# names(tdwg_final_glob)
# 
# sum(tdwg_final_glob$curr_nSp == 0)

## SPECIES DIVERSITY ============
# library(pscl)
# head(tdwg_final_glob)
# sum(tdwg_final_glob$megapalm_nsp == 0, na.rm = TRUE)
# 
# glob_sp_mod <- zeroinfl(megapalm_nsp ~ presNat_megaHerb_nSp, data = tdwg_final_glob, dist = "poisson")
# 
# 
# glob_sp_mod2 <- glm(megapalm_nsp ~ presNat_megaHerb_nSp+ scale(globalPC1) + scale(globalPC2) + scale(globalPC3) + scale(lgm_ens_Tano) + scale(lgm_ens_Pano)+ THREEREALM, data = tdwg_final_glob, family = "poisson")
# glob_sp_mod3 <- glm(megapalm_nsp ~ curr_megaHerb_nSp+ scale(globalPC1) + scale(globalPC2) + scale(globalPC3) + scale(lgm_ens_Tano) + scale(lgm_ens_Pano) + THREEREALM, data = tdwg_final_glob, family = "poisson")
# AIC(glob_sp_mod2)
# 
# nw_sp_mod1 <- glm(megapalm_nsp ~ presNat_megaHerb_nSp + scale(regionalPC1) + scale(regionalPC2) + scale(regionalPC3) + scale(lgm_ens_Tano) + scale(lgm_ens_Pano), data = tdwg_final_nw, family = "poisson")
# nw_sp_mod2 <- glm(megapalm_nsp ~ curr_megaHerb_nSp+ scale(regionalPC1) + scale(regionalPC2) + scale(regionalPC3) + scale(lgm_ens_Tano) + scale(lgm_ens_Pano), data = tdwg_final_nw, family = "poisson")
# AIC(nw_sp_mod1)
# AIC(nw_sp_mod2)
# 
# 
# nw_prop_mod <- glm(propMegaPalm ~ propMegaMam_presnat+ scale(regionalPC1) + scale(regionalPC2) + scale(regionalPC3) + scale(lgm_ens_Tano) + scale(lgm_ens_Pano), data = tdwg_final_nw, family = "binomial")
# nw_prop_mod2 <- glm(propMegaPalm ~ propMegaMam_curr+ scale(regionalPC1) + scale(regionalPC2) + scale(regionalPC3) + scale(lgm_ens_Tano) + scale(lgm_ens_Pano), data = tdwg_final_nw, family = "binomial")
# summary(step(nw_prop_mod2))


# ## CATEGORIES
# hist(tdwg_final$megaSpLoss)
# tdwg_final$propMegaLoss <- tdwg_final$megaSpLoss / tdwg_final$presNat_mega_nSp
# summary(lm(medianFruitLengthFilled ~ log(curr_medianBodySize)*megaSpLoss, data = tdwg_final))
# 
# plot(medianFruitLengthFilled ~ log(curr_medianBodySize)*megaSpLoss, data = tdwg_final)
# 
# summary(lm(propMegaMam_presnat ~ propMegaPalm, data = subset(tdwg_final, propMegaPalm > 0)))
# summary(lm(propMegaMam_curr ~ propMegaPalm, data = subset(tdwg_final, propMegaPalm > 0)))
# 
# plot(propMegaMam_curr ~ propMegaPalm, data = subset(tdwg_final, propMegaPalm > 0))
# abline(lm(propMegaMam_curr ~ propMegaPalm, data = subset(tdwg_final, propMegaPalm > 0)))
# 
# plot(propMegaMam_presnat ~ propMegaPalm, data = subset(tdwg_final, propMegaPalm > 0))
# abline(lm(propMegaMam_presnat ~ propMegaPalm, data = subset(tdwg_final, propMegaPalm > 0))) # areas with more megapalms will also have higher palm fruit sizes
# #subset(tdwg_final, propMegaPalm > 0)
# ggplot(aes(y = propMegaMam_presnat, x = propMegaPalm, colour = THREEREALM), data = tdwg_final) + geom_point() + geom_smooth(method = "lm")

# For areas with any megafaunal palms, areas with a greater proportion of megafaunal palms also correlated with a greater proportion of megafaunal mammals