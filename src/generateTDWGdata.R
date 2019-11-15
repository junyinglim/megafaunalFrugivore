## Generate body size and fruit size at the scale of botanical countries
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
library(sp); library(rgdal)
library(car); library(MuMIn)
library(readxl)
library(picante)
source(file.path(src.dir, "ggmodavg.R"))

## Import palm dataset ========================
palm_occ <- read.csv(file.path(data.dir, "palms_in_tdwg3.csv"))
palm_trait <- read.csv(file.path(data.dir, "PalmTraits_10.csv"))
palm_trait$SpecName <- gsub(palm_trait$SpecName, pattern= " ", replacement = "_")

## Import tdwg dataset ========================
#tdwg_env <- read.csv(file.path(data.dir, "TDWG_Environment_AllData_2019Jan.csv"))
tdwg_env <- read.csv(file.path(data.dir, "TDWG_Environment_AllData_2019Feb.csv"))

## Import diet dataset ========================
# Import MammalDiet v1.0 dataset
mammalDiet <- read.delim(file.path(frug.dir, "mammalDiet_v1.0.txt"))
mammalDiet$SpecName <- paste(mammalDiet$Genus, mammalDiet$Species, sep = "_")

# Import phylacine dataset
phylacine.dir <- file.path(frug.dir, "PHYLACINE")
phylacine_trait <- read.csv(file.path(phylacine.dir, "Trait_data.csv"))
spatial_metadata <- read.csv(file.path(phylacine.dir, "Spatial_metadata.csv"))

phylacine_trait <-  merge(phylacine_trait, spatial_metadata[c("Binomial.1.2", "Number.Cells.Current.Range", "Number.Cells.Present.Natural.Range")], by = "Binomial.1.2")

extinct_herb_splist <- subset(phylacine_trait, IUCN.Status.1.2 == "EP" & Diet.Plant == 100)$Binomial.1.2
extinct_herb_all_splist <- subset(phylacine_trait, IUCN.Status.1.2 == "EP" & Diet.Plant<= 100 & Diet.Plant >= 50)$Binomial.1.2

## Import mammal occurrence datasets ========================
# Import mammal occ datasets
mammal_presnat_occ <- readRDS(file.path(data.dir, "mammal_presnat_occ_raw.rds"))
mammal_presnat_facherb_occ <- readRDS(file.path(data.dir, "mammal_presnat_facherb_occ_raw.rds"))
mammal_nochange_occ <- readRDS(file.path(data.dir, "mammal_nochange_occ_raw.rds"))
mammal_curr_occ <- readRDS(file.path(data.dir, "mammal_current_occ_raw.rds"))

# Import data with false positives flagged
mammal_flagged_occ <- read.csv(file.path(data.dir, "flaggedPresentNaturalRanges_Faurby.csv"))
mammal_facherb_flagged_occ <- read.csv(file.path(data.dir, "flaggedPresentNaturalFacHerbRanges_Faurby.csv"))

# Combine datasets
mammal_presnat_occ2 <- Reduce(f = "rbind", x = list(mammal_presnat_occ[,-1], mammal_presnat_facherb_occ[,-1], mammal_nochange_occ[,-1]))
mammal_curr_occ2 <- rbind(mammal_curr_occ[,-1], mammal_nochange_occ[,-1])

# Remove false positives from present-natural occurrences
falsePos <- subset(rbind(mammal_flagged_occ,mammal_facherb_flagged_occ), Keep == 0)[c("LEVEL_3_CO", "SpecName")]  # define which taxon occurrences should be removed 

flagFalsePositives <- function(x,y, col.target = c("LEVEL_3_CO", "SpecName"), value.name = "falsePos"){
  # x: dataframe of interest with duplicates to be flagged
  # y: dataframe containing a set of forbidden combinations
  # col.target = columns that specify unique combinations
  # value.name = new column name that flags duplicates that appear in y
  col.target = c("LEVEL_3_CO", "SpecName")
  forbComb <- paste0(y[[col.target[1]]], y[[col.target[2]]])
  estComb <- paste0(x[[col.target[1]]], x[[col.target[2]]])
  # if there is no match (i.e., NA), flag as "no", otherwise "yes"
  x[[value.name]] <- ifelse(is.na(match(estComb, forbComb)), "No", "Yes") 
  return(x)
}

mammal_presnat_occ3 <- flagFalsePositives(x = mammal_presnat_occ2, y = falsePos)
mammal_curr_occ3 <- flagFalsePositives(x = mammal_curr_occ2, y = falsePos)

# Flag tertiary frugivores
tertFrug_spList <- subset(mammalDiet, Fruit == 3)$SpecName

# Remove tertiary frugivores and duplicate rows
mammal_presnat_occ4 <- subset(mammal_presnat_occ3, falsePos == "No" & (!SpecName %in% tertFrug_spList))
mammal_presnat_occ_final <- mammal_presnat_occ4[!duplicated(mammal_presnat_occ4),]

mammal_curr_occ4 <- subset(mammal_curr_occ3, falsePos == "No" & (!SpecName %in% tertFrug_spList))
mammal_curr_occ_final <- mammal_curr_occ4[!duplicated(mammal_curr_occ4),]

# Only include areas that contain palms
palm_tdwg_list <- unique(palm_occ$Area_code_L3)
mammal_presnat_comb_occ <- subset(mammal_presnat_occ_final, LEVEL_3_CO %in% palm_tdwg_list)
mammal_curr_comb_occ <- subset(mammal_curr_occ_final, LEVEL_3_CO %in% palm_tdwg_list)

# Export data
write.csv(mammal_presnat_comb_occ, file.path(data.dir, "mammal_presnat_occ.csv"), row.names = FALSE)
write.csv(mammal_curr_comb_occ, file.path(data.dir, "mammal_curr_occ.csv"), row.names = FALSE)

## Explore sensitivity of extinct mammals
mammalDiet$SpecName <- paste(mammalDiet$Genus, mammalDiet$Species, sep = "_")

mammalDietFamilyList <- unique(mammalDiet$Family)
mammalDietFam_summary <- ddply(.data = mammalDiet, .variables = .(Family),
                               summarize,
                               Order = Order[1],
                               famPropFrug =  sum(Fruit %in% 1:2, na.rm = T)/length(Fruit),
                               famPropOblgFrug = sum(Fruit == 1, na.rm = T)/ length(Fruit))

mammalDietOrder_summary <- ddply(.data = mammalDiet, .variables = .(Order),
                                 summarize,
                                 orderPropFrug =  sum(Fruit %in% 1:2, na.rm = T)/length(Fruit),
                                 orderPropOblgFrug = sum(Fruit == 1, na.rm = T)/ length(Fruit))

propFrugHigherTaxa <- merge(mammalDietFam_summary, mammalDietOrder_summary)

phylacine_trait$Family <- toupper(phylacine_trait$Family.1.2)
phylacine_trait$Order <- toupper(phylacine_trait$Order.1.2)
ep_frug <- merge(subset(phylacine_trait, IUCN.Status.1.2 == "EP" & Diet.Plant >= 50), y = propFrugHigherTaxa,
      by =c("Order", "Family"), all.x = TRUE)

classifyFrug <- function(x){
  frugCons <- NA
  if( is.na(x$famPropFrug) & is.na(x$orderPropFrug) ){
    frugCons <- 0
  } else if( is.na(x$famPropFrug) & x$orderPropFrug >= 0.75 ){
    frugCons <- 1
  } else if( x$famPropFrug >= 0.75  ){
    frugCons <- 1
  } else {
    frugCons <- 0
  }
  data.frame(frugCons)
}

ep_frug_class <- ddply(.data = ep_frug, .variables=.(Binomial.1.2), .fun = classifyFrug)


## Genus-level gap filling for palm fruit sizes ======================== 
# Define function for gap filling
gapFill <- function(x){
  if( sum(is.na(x$AverageFruitLength_cm)) == 0 ){
    # If the genus is complete in data, then just return values
    return(x)  
  } else {
    # Else, assign the empty values with the genus mean
    genusMean <- mean(x$AverageFruitLength_cm, na.rm = TRUE)
    x$AverageFruitLength_cm[is.na(x$AverageFruitLength_cm)] <- genusMean
    return(x)
  }
}

# Implement gap filling
palm_trait_subset <- palm_trait[c("PalmTribe","accGenus", "SpecName","AverageFruitLength_cm")]
palm_trait_gapfill  <- ddply(.data = palm_trait_subset,
                             .variables = .(accGenus),
                             .fun = gapFill)
names(palm_trait_gapfill)[names(palm_trait_gapfill) == "AverageFruitLength_cm"] <- "AverageFruitLength_cm_filled"
# Note: Only one species from a monotypic genus not filled this way "Butyagrus nabonnandii". However, species is a sterile hybrid between queen palm (Syagrus romanzoffiana) and pindo palm (Butia odora).

# Merge the gap filled data with the main dataset
palm_trait <- merge(palm_trait,
                    palm_trait_gapfill[c("SpecName", "AverageFruitLength_cm_filled")],
                    by = "SpecName")

## Calculate tdwg level metrics ======================== 
# Merge trait data to palm occurrence checklist
palm_occ_trait <- merge(palm_occ, palm_trait[c("SpecName", "AverageFruitLength_cm", "AverageFruitLength_cm_filled")])

# Exclude species that are not animal dispersed
palm_occ_trait <- palm_occ_trait[!palm_occ_trait$SpecName %in% c("Lodoicea_maldivica", "Cocos_nucifera", "Nypa_fruticans"),]
palm_occ_trait <- subset(palm_occ_trait, !is.na(AverageFruitLength_cm_filled)) # removes the two occurrences of Butyagrus nabonnandii

# Find list of countries that the mammals occur in
mammal_countrylist <- Reduce(union, list(unique(mammal_presnat_occ$LEVEL_3_CO),
                                         unique(mammal_nochange_occ$LEVEL_3_CO),
                                         unique(mammal_curr_occ$LEVEL_3_CO)))

# Find list of countries where palms and mammals overlap
mammal_palm_intersect <- intersect(unique(palm_occ_trait$Area_code_L3), mammal_countrylist) # 190 countries in total

# Calculate maximum and median fruit size
tdwg_fruit_summary <- ddply(.data = subset(palm_occ_trait, Area_code_L3 %in% mammal_palm_intersect),
                        .variables = .(Area_code_L3),
                        .fun = summarise,
                        # No gap filling
                        meanFruitLength = mean(AverageFruitLength_cm, na.rm = T),
                        medianFruitLength = median(AverageFruitLength_cm, na.rm = T),
                        max95FruitLength = quantile(AverageFruitLength_cm,
                                                    probs = 0.95, na.rm = T, type = 8 ),
                        # Gap filling
                        medianFruitLengthFilled = median(AverageFruitLength_cm_filled, na.rm = T),
                        meanFruitLengthFilled = mean(AverageFruitLength_cm_filled, na.rm = T),
                        max95FruitLengthFilled = quantile(AverageFruitLength_cm_filled,
                                                          probs = 0.95, na.rm = T, type = 8 ),
                        sdLogFruitLengthFilled = sd(log(AverageFruitLength_cm_filled), na.rm = T),
                        megapalm_nsp = sum(AverageFruitLength_cm_filled > 4, na.rm = T),
                        palm_nSp = length(AverageFruitLength_cm))
names(tdwg_fruit_summary)[names(tdwg_fruit_summary) == "Area_code_L3"] <- "LEVEL_3_CO"
write.csv(subset(palm_occ_trait, Area_code_L3 %in% mammal_palm_intersect), file = file.path(res.dir, "tdwg_palm_occ_trait.csv"), row.names = F)

# Calculate maximum and median body sizes of present natural mammal assemblages
mammal_presnat_occ_trait <- merge(mammal_presnat_comb_occ, phylacine_trait, by.x = "SpecName", by.y = "Binomial.1.2", all.x = TRUE)
mammal_presnat_occ_trait <- merge(mammal_presnat_occ_trait, ep_frug_class, by.x = "SpecName", by.y = "Binomial.1.2", all.x = TRUE)
mammal_presnat_occ_trait$frugCons[is.na(mammal_presnat_occ_trait$frugCons)] <- 1

tdwg_presnat_summary <-
  ddply(.data = subset(mammal_presnat_occ_trait, LEVEL_3_CO %in% mammal_palm_intersect),
        .variables = .(LEVEL_3_CO),
        .fun = summarize,
        presNat_meanBodySize = mean(Mass.g, na.rm = T),
        presNat_medianBodySize = quantile(Mass.g, probs = 0.5, na.rm = T),
        presNat_max95BodySize = quantile(Mass.g, probs = 0.95, type = 8, na.rm = T ),
        presNat_medianBodySize_cons = quantile(Mass.g[frugCons == 1], probs = 0.5, na.rm = T),
        presNat_max95BodySize_cons = quantile(Mass.g[frugCons == 1], probs = 0.95, na.rm = T),
        presNat_sdBodySize = sd(log(Mass.g), na.rm = T),
        presNat_nSp = length(unique(SpecName)),
        presNat_meso_nSp = length(unique(SpecName[Mass.g > 10000])),
        presNat_mega_nSp = length(unique(SpecName[Mass.g > 44000])))
write.csv(subset(mammal_presnat_occ_trait, LEVEL_3_CO %in% mammal_palm_intersect), file.path(res.dir, "mammal_presnat_occ_trait.csv"), row.names = F)

# Current the maximum and median body sizes of current mammal assemblages
mammal_curr_occ_trait <- merge(mammal_curr_comb_occ, phylacine_trait, by.x = "SpecName", by.y = "Binomial.1.2", all.x = TRUE)

tdwg_curr_summary <- 
  ddply(.data = subset(mammal_curr_occ_trait, LEVEL_3_CO %in% mammal_palm_intersect),
        .variables = .(LEVEL_3_CO),
        .fun = summarize,
        curr_meanBodySize = mean(Mass.g, na.rm = T),
        curr_medianBodySize = quantile(Mass.g, probs = 0.5, type = 8, na.rm = T),
        curr_max95BodySize = quantile(Mass.g, probs = 0.95, type = 8, na.rm = T),
        curr_sdBodySize = sd(log(Mass.g), na.rm = T),
        curr_nSp = length(unique(SpecName)),
        curr_meso_nSp = length(unique(SpecName[Mass.g > 10000])),
        curr_mega_nSp = length(unique(SpecName[Mass.g > 44000]))
)
# NOTE: Pteropus niger (Mascarene fruit bat) is found on both Mauritius and Reunion but is only on Mauritius for the current dataset, as a result, there are no mammals on Reunion in the current case (Pteropus) but 2 in the present-natural dataset
write.csv(subset(mammal_curr_occ_trait, LEVEL_3_CO %in% mammal_palm_intersect), file.path(res.dir, "mammal_curr_occ_trait.csv"), row.names = F)

## Simulate extinction ====================s
mammal_curr_splist <- unique(mammal_curr_occ_trait$SpecName)
mammal_curr_sp <- subset(phylacine_trait, Binomial.1.2 %in% mammal_curr_splist)[c("Binomial.1.2", "IUCN.Status.1.2")]

ctmc_pExt <- read.csv(file.path(data.dir, "ctmc_pExt.csv"))

# extinctionProb <- data.frame(
#   IUCN.Status.1.2 = c("LC", "NT", "VU", "EN", "CR", "DD"),
#   #extP50 = c(0.00005, 0.004, 0.0513, 0.4276, 0.9688, 0.00005)
#   extP50 = c(0.0009, 0.0071, 0.0513, 0.4276, 0.9688, 0.0009)
# )
# Values from Moore are 0.00005, 0.004, 0.05, 0.42, 0.97, 0.00005
# Mooers uses the IUCN 2001 designations, i.e., P(ext_CR, 10 years) = 0.5; P(ext_EN, 20 years) = 0.2 and P(ext_VU, 100 years) = 0.1
# Rearranging the equations P(r, t) = 1 - exp(-rt) , the rates for CR = -log(0.5)/10, EN = -log(0.8)/20, VU = -log(0.9) / 100
# For LC and NT, they just assumed P(ext_LC, 100) of 0.0001, ext_LC = -log(0.9999)/100
# Using these same rates, you can standardize them for the same time horizons, we then get
# P(CR,50) = 0.96875, P(EN,50) = 0.427; P(VU,50) = 0.0513
# Ext prob after 50 years (Davis et al 2018) = 0.0009, 0.0071, 0.0513, 0.4276, 0.9688, 0.0009

mammal_curr_sp_status <- merge(x = mammal_curr_sp, y = ctmc_pExt, by.x = "IUCN.Status.1.2", by.y = "IUCN.Status")

simulateExtinction <- function(df, col){
  mammal_futr_splist <- vector()  
  for(i in 1:nrow(df)){
    mammal_futr_splist[i] <- sample(x = c(df$Binomial.1.2[i], NA),
                                    prob = c(1-df[[col]][i], df[[col]][i]),
                                    size = 1)
  } 
  return(mammal_futr_splist)
}

# Simulate 1000 extinction scenarios
mammal_dimarco_futr_splist <- list()
mammal_hoffmann_futr_splist <- list()
for(i in 1:1000){
  pb = txtProgressBar(min = 0., max = 1000, initial = 0, style = 3)
  setTxtProgressBar(pb, i)
  mammal_dimarco_futr_splist[[i]] <- simulateExtinction(mammal_curr_sp_status, col = "dimarco")
  mammal_hoffmann_futr_splist[[i]] <- simulateExtinction(mammal_curr_sp_status, col = "hoffmann")
}
close(pb)

# Calculate body size for the 1000 extinction simulations
tdwg_futr_meanBodySize_reps <- list()
for(i in 1:1000){
  pb = txtProgressBar(min = 0., max = 1000, initial = 0, style = 3)
  setTxtProgressBar(pb, i)
  tdwg_futr_meanBodySize_reps [[i]] <- 
    ddply(.data = subset(mammal_curr_occ_trait,
                         LEVEL_3_CO %in% mammal_palm_intersect),
        .variables = .(LEVEL_3_CO),
        .fun = summarize,
        futr_dimarco_medianBodySize = median(Mass.g[SpecName %in% mammal_dimarco_futr_splist[[i]] ] , na.rm = T),
        futr_dimarco_maxBodySize = quantile(Mass.g[SpecName %in% mammal_dimarco_futr_splist[[i]]],
                                    probs = 0.95, na.rm = T),
        futr_hoffmann_medianBodySize = median(Mass.g[SpecName %in% mammal_hoffmann_futr_splist[[i]] ] , na.rm = T),
        futr_hoffmann_maxBodySize = quantile(Mass.g[SpecName %in% mammal_hoffmann_futr_splist[[i]]],
                                            probs = 0.95, na.rm = T))
}
close(pb)

tdwg_futr_meanBodySize <- Reduce(tdwg_futr_meanBodySize_reps, f = "rbind")
saveRDS(tdwg_futr_meanBodySize, file = file.path(res.dir, "tdwg_futrBS_reps.rds"))

tdwg_futr_meanBodySize_summary <- 
  ddply(.data = tdwg_futr_meanBodySize,
        .variable = .(LEVEL_3_CO),
        .fun = summarize,
        futr_dimarco_medianBodySize = median(futr_dimarco_medianBodySize, na.rm = T),
        futr_dimarco_maxBodySize = median(futr_dimarco_maxBodySize, na.rm = T),
        futr_hoffmann_medianBodySize = median(futr_hoffmann_medianBodySize, na.rm = T),
        futr_hoffmann_maxBodySize = median(futr_hoffmann_maxBodySize, na.rm = T),
        nreps = length(LEVEL_3_CO))

tdwg_curr_summary <- merge(tdwg_curr_summary, tdwg_futr_meanBodySize_summary, by = "LEVEL_3_CO", all = TRUE)

# Merge present natural and current mammal assemblage summary statistics and fruit statistics
tdwg_mammal_all <- merge(tdwg_presnat_summary, tdwg_curr_summary, by = "LEVEL_3_CO", all = TRUE)
tdwg_res <- merge(tdwg_fruit_summary, tdwg_mammal_all, by = "LEVEL_3_CO")

# Remove Reunion from dataset
tdwg_res <- tdwg_res[!(is.na(tdwg_res$curr_medianBodySize) | is.na(tdwg_res$presNat_medianBodySize)),]

# Export dataset
write.csv(tdwg_res, file.path(res.dir, "tdwg_mammal.csv"), row.names = FALSE)


## Data handling and clean up ========================
# Calculate some derived summary statistics
tdwg_res$propMegaPalm <- tdwg_res$megapalm_nsp / tdwg_res$palm_nSp
tdwg_res$propMegaMam_curr <- tdwg_res$curr_mega_nSp / tdwg_res$curr_nSp
tdwg_res$propMegaMam_presnat <- tdwg_res$presNat_mega_nSp / tdwg_res$presNat_nSp
tdwg_res$deltaMedianBodySize = log(tdwg_res$presNat_medianBodySize) - log(tdwg_res$curr_medianBodySize)
tdwg_res$megaSpLoss <- tdwg_res$presNat_mega_nSp - tdwg_res$curr_mega_nSp 
tdwg_res$mesoSpLoss <- tdwg_res$presNat_meso_nSp - tdwg_res$curr_meso_nSp

# Merge with tdwg environmental data
tdwg_final <- merge(tdwg_res, tdwg_env, by = "LEVEL_3_CO", all.x = TRUE)
tdwg_final2 <- subset(tdwg_final, REALM_LONG %in% c("IndoMalay", "Australasia", "Neotropics", "Afrotropics"))

# Global level PCA
env.var <- paste0("bio", c(12, 15, 17, 1, 11, 4), "_mean")
  #paste0("bio", 1:19, "_mean")
#c("PREC_Sum", "PREC_CV", "P_drie_quart", "Tmean_mean", "T_cold_quart", "Temp_SD")
globalPCA <- prcomp(x = tdwg_final2[env.var], scale. = TRUE, center = TRUE)

tdwg_final2$globalPC1 <- globalPCA$x[,1]
tdwg_final2$globalPC2 <- globalPCA$x[,2]
tdwg_final2$globalPC3 <- globalPCA$x[,3]

globalPCA_varexp <- cumsum(globalPCA$sdev^2) /sum(globalPCA$sdev^2)
globalPCA_res <- rbind(globalPCA$rotation, globalPCA_varexp)
write.csv(roundNumbers(data.frame(globalPCA_res)), row.names = T, file = file.path(res.dir, "globalPCA_res.csv"))

# Regional PCA
tdwg_final_nw <- subset(tdwg_final2, REALM_LONG == "Neotropics")
tdwg_final_oww <- subset(tdwg_final2, REALM_LONG == "Afrotropics")
tdwg_final_owe <- subset(tdwg_final2, REALM_LONG %in% c("IndoMalay", "Australasia"))

nwPCA <- prcomp(x= tdwg_final_nw[env.var], scale. = TRUE, center = TRUE)
owwPCA <- prcomp(x= tdwg_final_oww[env.var], scale. = TRUE, center = TRUE)
owePCA <- prcomp(x= tdwg_final_owe[env.var], scale. = TRUE, center = TRUE)

nwPCA_varexp <- cumsum(nwPCA$sdev^2) / sum(nwPCA$sdev^2)
owwPCA_varexp <- cumsum(owwPCA$sdev^2) / sum(owwPCA$sdev^2)
owePCA_varexp <- cumsum(owePCA$sdev^2) / sum(owePCA$sdev^2)

nwPCA_res <- rbind(nwPCA$rotation, nwPCA_varexp)
write.csv(roundNumbers(data.frame(nwPCA_res)), file = file.path(res.dir, "nwPCA_res.csv"))
owwPCA_res <- rbind(owwPCA$rotation, owwPCA_varexp)
write.csv(roundNumbers(data.frame(owwPCA_res)), file.path(res.dir, "owwPCA_res.csv"))
owePCA_res <- rbind(owePCA$rotation, owePCA_varexp)
write.csv(roundNumbers(data.frame(owePCA_res)), file.path(res.dir, "owePCA_res.csv"))

tdwg_final_nw$regionalPC1 <- nwPCA$x[,1]
tdwg_final_nw$regionalPC2 <- nwPCA$x[,2]
tdwg_final_nw$regionalPC3 <- nwPCA$x[,3]

tdwg_final_oww$regionalPC1 <- owwPCA$x[,1]
tdwg_final_oww$regionalPC2 <- owwPCA$x[,2]
tdwg_final_oww$regionalPC3 <- owwPCA$x[,3]

tdwg_final_owe$regionalPC1 <- owePCA$x[,1]
tdwg_final_owe$regionalPC2 <- owePCA$x[,2]
tdwg_final_owe$regionalPC3 <- owePCA$x[,3]

pc.col <- c("LEVEL_3_CO", "regionalPC1", "regionalPC2", "regionalPC3")
regionalPCA <- Reduce(rbind, list(tdwg_final_nw[pc.col], tdwg_final_oww[pc.col], tdwg_final_owe[pc.col]))
tdwg_final2 <- merge(tdwg_final2, regionalPCA, by = "LEVEL_3_CO")

write.csv(tdwg_final2, file.path(data.dir,"tdwg_final.csv"))

## Generate some summary statistics
mammal_curr_occ_trait <- read.csv(file.path(res.dir, "mammal_curr_occ_trait.csv"))
mammal_presnat_occ_trait <- read.csv(file.path(res.dir, "mammal_presnat_occ_trait.csv"))
tdwg_final2 <- read.csv(file.path(data.dir,"tdwg_final.csv"))
phylacine.dir <- file.path(frug.dir, "PHYLACINE")
phylacine_trait <- read.csv(file.path(phylacine.dir, "Trait_data.csv"))

length(unique(subset(mammal_curr_occ_trait, LEVEL_3_CO %in% tdwg_final2$LEVEL_3_CO)$SpecName)) # 1604 taxa in the current 
length(unique(subset(mammal_presnat_occ_trait, LEVEL_3_CO %in% tdwg_final2$LEVEL_3_CO)$SpecName)) # 1811 in the present-natural

curr_status <- table(subset(phylacine_trait, Binomial.1.2 %in% unique(mammal_curr_occ_trait$SpecName))$IUCN.Status.1.2)
sum(curr_status[names(curr_status) %in% c("CR", "DD", "EN", "LC", "NT", "VU")])
pnat_status <- table(subset(phylacine_trait, Binomial.1.2 %in% unique(mammal_presnat_occ_trait$SpecName))$IUCN.Status.1.2)
sum(pnat_status[names(pnat_status) %in% c("CR", "DD", "EN", "LC", "NT", "VU")])

table(mammal_presnat_occ_trait$CONTINENT)
asd <- subset(mammal_presnat_occ_trait, IUCN.Status.1.2 == "EP" & CONTINENT == "ASIA-TROPICAL")
asdlsit <- unique(asd$SpecName)
test <- subset(phylacine_trait, Binomial.1.2 %in% asdlsit)
test[order(test$Mass.g), names(test) %in% c("Binomial.1.2", "Mass.g")]
