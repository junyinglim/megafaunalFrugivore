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
tdwg_final_nw <- subset(tdwg_final, REALM_LONG == "Neotropics")
tdwg_final_oww <- subset(tdwg_final, REALM_LONG == "Afrotropics")
tdwg_final_owe <- subset(tdwg_final, REALM_LONG %in% c("Australasia","IndoMalay"))

## Model averaging (Median Body Size) ========================
glob_curr_medBS_mod <- lm(log(medianFruitLengthFilled) ~ scale(log(curr_medianBodySize)) + scale(globalPC1) + scale(globalPC2) + scale(globalPC3) + scale(soilcount) + scale(lgm_ens_Tano) + scale(lgm_ens_Pano), data =tdwg_final_glob, na.action = "na.fail")
glob_pnat_medBS_mod <- update(glob_curr_medBS_mod, ~.-scale(log(curr_medianBodySize))+
                          scale(log(presNat_medianBodySize)))

glob_curr_medBS_moddr <- dredge(glob_curr_medBS_mod)
glob_curr_medBS_modavg <- model.avg(glob_curr_medBS_moddr)
glob_curr_medBS_modavg_summ <- summarizeRelImportance(glob_curr_medBS_modavg)

glob_pnat_medBS_moddr <- dredge(glob_pnat_medBS_mod)
glob_pnat_medBS_modavg <- model.avg(glob_pnat_medBS_moddr)
glob_pnat_medBS_modavg_summ <- summarizeRelImportance(glob_pnat_medBS_modavg)

# Create placeholder linear model to insert model averaged coefficients
# z2<- glob_curr_medBS_mod
# z2$coefficients <-  coefficients(glob_curr_medBS_modavg)
# library(piecewiseSEM)
# glob_curr_medBS_mod_presid <- data.frame(presid = resid(z2, type = "partial")[,1], var = scale(log(tdwg_final_glob$curr_medianBodySize)))
# ggplot(aes(y = presid, x = var), data = glob_curr_medBS_mod_presid) + geom_point() + geom_smooth(method = "lm")
# car::crPlot(model = z2, variable = "scale(log(curr_medianBodySize))")
# car::crPlot(model = glob_curr_medBS_mod, variable = "scale(log(curr_medianBodySize))")
# z3 <- partialResid(medianFruitLengthFilled ~ curr_medianBodySize, modelList = z2, data = tdwg_final_glob)

nw_curr_medBS_mod <- lm(log(medianFruitLengthFilled) ~ scale(log(curr_medianBodySize)) + scale(regionalPC1) + scale(regionalPC2) + scale(regionalPC3) + scale(soilcount) + scale(lgm_ens_Tano) + scale(lgm_ens_Pano), data = tdwg_final_nw, na.action = "na.fail")
nw_pnat_medBS_mod <- update(nw_curr_medBS_mod, ~.-scale(log(curr_medianBodySize))  + scale(log(presNat_medianBodySize)))
vif(nw_curr_medBS_mod)
vif(nw_pnat_medBS_mod)

nw_curr_medBS_moddr <- dredge(nw_curr_medBS_mod)
nw_curr_medBS_modavg <- model.avg(nw_curr_medBS_moddr)
nw_curr_medBS_modavg_summ <- summarizeRelImportance(nw_curr_medBS_modavg)

nw_pnat_medBS_moddr <- dredge(nw_pnat_medBS_mod)
nw_pnat_medBS_modavg <- model.avg(nw_pnat_medBS_moddr)
nw_pnat_medBS_modavg_summ <- summarizeRelImportance(nw_pnat_medBS_modavg)

oww_curr_medBS_mod <- lm(log(medianFruitLengthFilled) ~ scale(log(curr_medianBodySize)) + scale(regionalPC1) + scale(regionalPC2) + scale(regionalPC3) + scale(soilcount) + scale(lgm_ens_Tano) + scale(lgm_ens_Pano) + scale(curr_medRangeSize), data = tdwg_final_oww, na.action = "na.fail")
oww_pnat_medBS_mod <- update(oww_curr_medBS_mod, ~.-scale(log(curr_medianBodySize)) + scale(log(presNat_medianBodySize)))
vif(oww_curr_medBS_mod)
vif(oww_pnat_mod)

oww_curr_medBS_moddr <- dredge(oww_curr_medBS_mod)
oww_curr_medBS_modavg <- model.avg(oww_curr_medBS_moddr)
oww_curr_medBS_modavg_summ <- summarizeRelImportance(oww_curr_medBS_modavg)

oww_pnat_medBS_moddr <- dredge(oww_pnat_medBS_mod)
oww_pnat_medBS_modavg <- model.avg(oww_pnat_medBS_moddr)
oww_pnat_medBS_modavg_summ <- summarizeRelImportance(oww_pnat_medBS_modavg)

owe_curr_medBS_mod <- lm(log(medianFruitLengthFilled) ~ scale(log(curr_medianBodySize)) + scale(regionalPC1) + scale(regionalPC2) + scale(regionalPC3) + scale(soilcount) + scale(lgm_ens_Tano) + scale(lgm_ens_Pano) + scale(curr_medRangeSize), data = tdwg_final_owe, na.action = "na.fail")
owe_pnat_medBS_mod <- update(owe_curr_medBS_mod, ~.-scale(log(curr_medianBodySize)) + scale(log(presNat_medianBodySize)) )
vif(owe_curr_mod); vif(owe_pnat_mod)

owe_curr_medBS_moddr <- dredge(owe_curr_medBS_mod)
owe_curr_medBS_modavg <- model.avg(owe_curr_medBS_moddr)
owe_curr_medBS_modavg_summ <- summarizeRelImportance(owe_curr_medBS_modavg)
confint(owe_curr_medBS_modavg)
summary(owe_curr_medBS_modavg)
owe_pnat_medBS_moddr <- dredge(owe_pnat_medBS_mod)
owe_pnat_medBS_modavg <- model.avg(owe_pnat_medBS_moddr)
owe_pnat_medBS_modavg_summ <- summarizeRelImportance(owe_pnat_medBS_modavg)

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
# # Generate output table
# medianModSummary <- rbind(
#   summarizeLM(glob_curr_step_mod, scale = "Global", scenario = "Current"),
#   summarizeLM(glob_pnat_step_mod, scale = "Global", scenario = "Present-natural"),
#   summarizeLM(nw_curr_step_mod, scale = "Neotropics", scenario = "Current"),
#   summarizeLM(nw_pnat_step_mod, scale = "Neotropics", scenario = "Present-natural"),
#   summarizeLM(oww_curr_step_mod, scale = "Afrotropics", scenario = "Current"),
#   summarizeLM(oww_pnat_step_mod, scale = "Afrotropics", scenario = "Present-natural"),
#   summarizeLM(owe_curr_step_mod, scale = "Indotropics", scenario = "Current"),
#   summarizeLM(owe_pnat_step_mod, scale = "Indotropics", scenario = "Present-natural")
# )
# write.csv(medianModSummary, file.path(res.dir, "medianModSummary.csv"), fileEncoding="UTF-8",
#           row.names = F)

## MAXIMUM BODY SIZE ============
glob_curr_maxBS_mod <- lm(log(max95FruitLengthFilled) ~ scale(log(curr_max95BodySize)) + scale(globalPC1) + scale(globalPC2) + scale(globalPC3) + scale(soilcount) + scale(lgm_ens_Tano) + scale(lgm_ens_Pano), data =tdwg_final_glob, na.action = "na.fail")
glob_pnat_maxBS_mod <- update(glob_curr_maxBS_mod, ~.- scale(log(curr_max95BodySize))+ scale(log(presNat_max95BodySize)))
vif(glob_curr_maxBS_mod); vif(glob_pnat_maxBS_mod)

glob_curr_maxBS_moddr <- dredge(glob_curr_maxBS_mod)
glob_curr_maxBS_modavg <- model.avg(glob_curr_maxBS_moddr)
glob_curr_maxBS_modavg_summ <- summarizeRelImportance(glob_curr_maxBS_modavg)

glob_pnat_maxBS_moddr <- dredge(glob_pnat_maxBS_mod)
glob_pnat_maxBS_modavg <- model.avg(glob_pnat_maxBS_moddr)
glob_pnat_maxBS_modavg_summ <- summarizeRelImportance(glob_pnat_maxBS_modavg)

nw_curr_maxBS_mod <- lm(log(max95FruitLengthFilled) ~ scale(log(curr_max95BodySize)) + scale(regionalPC1) + scale(regionalPC2) + scale(regionalPC3) + scale(soilcount) + scale(lgm_ens_Tano) + scale(lgm_ens_Pano), data = tdwg_final_nw, na.action = "na.fail")
nw_pnat_maxBS_mod <- update(nw_curr_maxBS_mod, ~.- scale(log(curr_max95BodySize))+ scale(log(presNat_max95BodySize)))
vif(nw_curr_maxBS_mod); vif(nw_pnat_maxBS_mod)

nw_curr_maxBS_moddr <- dredge(nw_curr_maxBS_mod)
nw_curr_maxBS_modavg <- model.avg(nw_curr_maxBS_moddr)
nw_curr_maxBS_modavg_summ <- summarizeRelImportance(nw_curr_maxBS_modavg)

nw_pnat_maxBS_moddr <- dredge(nw_pnat_maxBS_mod)
nw_pnat_maxBS_modavg <- model.avg(nw_pnat_maxBS_moddr)
nw_pnat_maxBS_modavg_summ <- summarizeRelImportance(nw_pnat_maxBS_modavg)

oww_curr_maxBS_mod <- lm(log(max95FruitLengthFilled) ~ scale(log(curr_max95BodySize)) + scale(regionalPC1) + scale(regionalPC2) + scale(regionalPC3) + scale(soilcount) + scale(lgm_ens_Tano) + scale(lgm_ens_Pano), data = tdwg_final_oww, na.action = "na.fail")
oww_pnat_maxBS_mod <- update(oww_curr_maxBS_mod, ~.-scale(log(curr_max95BodySize)) + scale(log(presNat_max95BodySize)))
vif(oww_curr_maxBS_mod); vif(oww_pnat_maxBS_mod)

oww_curr_maxBS_moddr <- dredge(oww_curr_maxBS_mod)
oww_curr_maxBS_modavg <- model.avg(oww_curr_maxBS_moddr)
oww_curr_maxBS_modavg_summ <- summarizeRelImportance(oww_curr_maxBS_modavg)

oww_pnat_maxBS_moddr <- dredge(oww_pnat_maxBS_mod)
oww_pnat_maxBS_modavg <- model.avg(oww_pnat_maxBS_moddr)
oww_pnat_maxBS_modavg_summ <- summarizeRelImportance(oww_pnat_maxBS_modavg)

owe_curr_maxBS_mod <- lm(log(max95FruitLengthFilled) ~ scale(log(curr_max95BodySize)) + scale(regionalPC1) + scale(regionalPC2) + scale(regionalPC3) + scale(soilcount) + scale(lgm_ens_Tano) + scale(lgm_ens_Pano), data = tdwg_final_owe, na.action = "na.fail")
owe_pnat_maxBS_mod <- update(owe_curr_maxBS_mod, ~.-scale(log(curr_max95BodySize)) + scale(log(presNat_max95BodySize)))
vif(owe_curr_maxBS_mod); vif(owe_pnat_maxBS_mod)

owe_curr_maxBS_moddr <- dredge(owe_curr_maxBS_mod)
owe_curr_maxBS_modavg <- model.avg(owe_curr_maxBS_moddr)
owe_curr_maxBS_modavg_summ <- summarizeRelImportance(owe_curr_maxBS_modavg)

owe_pnat_maxBS_moddr <- dredge(owe_pnat_maxBS_mod)
owe_pnat_maxBS_modavg <- model.avg(owe_pnat_maxBS_moddr)
owe_pnat_maxBS_modavg_summ <- summarizeRelImportance(owe_pnat_maxBS_modavg)

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

# maxModSummary <- rbind(
#   summarizeLM(glob_curr_step_mod, scale = "Global", scenario = "Current"),
#   summarizeLM(glob_pnat_step_mod, scale = "Global", scenario = "Present-natural"),
#   summarizeLM(nw_curr_step_mod, scale = "Neotropics", scenario = "Current"),
#   summarizeLM(nw_pnat_step_mod, scale = "Neotropics", scenario = "Present-natural"),
#   summarizeLM(oww_curr_step_mod, scale = "Afrotropics", scenario = "Current"),
#   summarizeLM(oww_pnat_step_mod, scale = "Afrotropics", scenario = "Present-natural"),
#   summarizeLM(owe_curr_step_mod, scale = "Indotropics", scenario = "Current"),
#   summarizeLM(owe_pnat_step_mod, scale = "Indotropics", scenario = "Present-natural")
# )
# write.csv(maxModSummary, file.path(res.dir, "maxModSummary.csv"), fileEncoding="UTF-8",
#           row.names = F)


## DISPERSION IN BODY SIZE ============
tdwg_final_glob_disp <- subset(tdwg_final_glob, palm_nSp > 1 & curr_nSp > 1 & presNat_nSp > 1)
tdwg_final_nw_disp <- subset(tdwg_final_glob_disp, REALM_LONG == "Neotropics")
tdwg_final_oww_disp <- subset(tdwg_final_glob_disp, REALM_LONG == "Afrotropics")
tdwg_final_owe_disp <- subset(tdwg_final_glob_disp, REALM_LONG %in% c("Australasia","IndoMalay"))

glob_curr_dispBS_mod <- lm(dispFruitLengthFilled ~ scale(curr_dispBodySize) + scale(globalPC1) + scale(globalPC2) + scale(globalPC3) + scale(soilcount) + scale(lgm_ens_Tano) + scale(lgm_ens_Pano), data = tdwg_final_glob_disp, na.action = "na.fail")
glob_pnat_dispBS_mod <- update(glob_curr_dispBS_mod, ~.- scale(curr_dispBodySize)+ scale(presNat_dispBodySize))
vif(glob_curr_dispBS_mod); vif(glob_pnat_dispBS_mod)

glob_curr_dispBS_moddr <- dredge(glob_curr_dispBS_mod)
glob_curr_dispBS_modavg <- model.avg(glob_curr_dispBS_moddr)
glob_curr_dispBS_modavg_summ <- summarizeRelImportance(glob_curr_dispBS_modavg)

glob_pnat_dispBS_moddr <- dredge(glob_pnat_dispBS_mod)
glob_pnat_dispBS_modavg <- model.avg(glob_pnat_dispBS_moddr)
glob_pnat_dispBS_modavg_summ <- summarizeRelImportance(glob_pnat_dispBS_modavg)

nw_curr_dispBS_mod <- lm(dispFruitLengthFilled ~ scale(curr_dispBodySize) + scale(regionalPC1) + scale(regionalPC2) + scale(regionalPC3) + scale(soilcount) + scale(lgm_ens_Tano) + scale(lgm_ens_Pano), data =tdwg_final_nw_disp, na.action = "na.fail")
nw_pnat_dispBS_mod <- update(nw_curr_dispBS_mod, ~.- scale(curr_dispBodySize)+ scale(presNat_dispBodySize))
vif(nw_curr_dispBS_mod); vif(nw_pnat_dispBS_mod)

nw_curr_dispBS_moddr <- dredge(nw_curr_dispBS_mod)
nw_curr_dispBS_modavg <- model.avg(nw_curr_dispBS_moddr)
nw_curr_dispBS_modavg_summ <- summarizeRelImportance(nw_curr_dispBS_modavg)

nw_pnat_dispBS_moddr <- dredge(nw_pnat_dispBS_mod)
nw_pnat_dispBS_modavg <- model.avg(nw_pnat_dispBS_moddr)
nw_pnat_dispBS_modavg_summ <- summarizeRelImportance(nw_pnat_dispBS_modavg)

oww_curr_dispBS_mod <- lm(dispFruitLengthFilled ~ scale(curr_dispBodySize) + scale(regionalPC1) + scale(regionalPC2) + scale(regionalPC3) + scale(soilcount) + scale(lgm_ens_Tano) + scale(lgm_ens_Pano), data =tdwg_final_oww_disp, na.action = "na.fail")
oww_pnat_dispBS_mod <- update(oww_curr_dispBS_mod, ~.- scale(curr_dispBodySize)+ scale(presNat_dispBodySize))
vif(oww_curr_dispBS_mod); vif(oww_pnat_dispBS_mod)

oww_curr_dispBS_moddr <- dredge(oww_curr_dispBS_mod)
oww_curr_dispBS_modavg <- model.avg(oww_curr_dispBS_moddr)
oww_curr_dispBS_modavg_summ <- summarizeRelImportance(oww_curr_dispBS_modavg)

oww_pnat_dispBS_moddr <- dredge(oww_pnat_dispBS_mod)
oww_pnat_dispBS_modavg <- model.avg(oww_pnat_dispBS_moddr)
oww_pnat_dispBS_modavg_summ <- summarizeRelImportance(oww_pnat_dispBS_modavg)

owe_curr_dispBS_mod <- lm(dispFruitLengthFilled ~ scale(curr_dispBodySize) + scale(regionalPC1) + scale(regionalPC2) + scale(regionalPC3) + scale(soilcount) + scale(lgm_ens_Tano) + scale(lgm_ens_Pano), data =tdwg_final_owe_disp, na.action = "na.fail")
owe_pnat_dispBS_mod <- update(owe_curr_dispBS_mod, ~.- scale(curr_dispBodySize)+ scale(presNat_dispBodySize))
vif(owe_curr_dispBS_mod); vif(owe_pnat_dispBS_mod)

owe_curr_dispBS_moddr <- dredge(owe_curr_dispBS_mod)
owe_curr_dispBS_modavg <- model.avg(owe_curr_dispBS_moddr)
owe_curr_dispBS_modavg_summ <- summarizeRelImportance(owe_curr_dispBS_modavg)

owe_pnat_dispBS_moddr <- dredge(owe_pnat_dispBS_mod)
owe_pnat_dispBS_modavg <- model.avg(owe_pnat_dispBS_moddr)
owe_pnat_dispBS_modavg_summ <- summarizeRelImportance(owe_pnat_dispBS_modavg)

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