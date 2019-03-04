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
library(car); library(MuMIn)
library(readxl)
library(picante)

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
palm_occ_trait <- palm_occ_trait[!palm_occ_trait$SpecName %in% c("Lodoicea_maldivica", "Cocos_nucifera"),]
palm_occ_trait <- subset(palm_occ_trait, !is.na(AverageFruitLength_cm_filled)) # removes the two occurrences of Butyagrus nabonnandii

# Find list of countries that the mammals occur in
mammal_countrylist <- Reduce(union, list(unique(mammal_presnat_occ$LEVEL_3_CO),
                                         unique(mammal_nochange_occ$LEVEL_3_CO),
                                         unique(mammal_curr_occ$LEVEL_3_CO)))

# Find list of countries where palms and mammals overlap
mammal_palm_intersect <- intersect(unique(palm_occ_trait$Area_code_L3), mammal_countrylist) # 190 countries in total

# Calculate mean fruit size of TDWG units
tdwg_meanFruit <- ddply(.data = subset(palm_occ_trait, Area_code_L3 %in% mammal_palm_intersect),
                        .variables = .(Area_code_L3),
                        .fun = summarise,
                        meanFruitLength = mean(AverageFruitLength_cm, na.rm = T),
                        medianFruitLengthFilled = median(AverageFruitLength_cm_filled, na.rm = T),
                        meanFruitLengthFilled = mean(AverageFruitLength_cm_filled, na.rm = T),
                        maxFruitLengthFilled = max(AverageFruitLength_cm_filled, na.rm = T),
                        max95FruitLengthFilled = quantile(AverageFruitLength_cm_filled,
                                                          probs = 0.95, na.rm = T ),
                        minFruitLengthFilled  = min(AverageFruitLength_cm_filled, na.rm = T),
                        rangeFruitLengthFilled = log(maxFruitLengthFilled)-log(minFruitLengthFilled),
                        dispFruitLengthFilled = mean(log(AverageFruitLength_cm_filled) - 
                                                       log(medianFruitLengthFilled), na.rm = T) ,
                        sdLogFruitLengthFilled = sd(log(AverageFruitLength_cm_filled), na.rm = T),
                        megapalm_nsp = sum(AverageFruitLength_cm_filled > 4, na.rm = T),
                        palm_nSp = length(AverageFruitLength_cm))
names(tdwg_meanFruit)[names(tdwg_meanFruit) == "Area_code_L3"] <- "LEVEL_3_CO"

# # Calculate z-scores for fruit size range, regional source pools
# palm_occ_trait2 <- merge(y = subset(palm_occ_trait, Area_code_L3 %in% mammal_palm_intersect),
#                          x = tdwg_env[c("LEVEL_3_CO", "THREEREALM", "REALM_LONG")],
#                          by.y = "Area_code_L3", by.x = "LEVEL_3_CO")
# palm_occ_trait2 <- subset(palm_occ_trait2, REALM_LONG %in% c("Neotropics", "Afrotropics", "IndoMalay", "Australasia"))
# 
# ses_range <- function(x, prefix, value.var){
#   adj_mat <- acast(LEVEL_3_CO ~ SpecName, value.var = value.var, data = x, fill = 0) # only contains species from subset, country (rows), species (columns)
#   obs <- apply(adj_mat, MARGIN = 1, FUN = function(x){ diff(range(log(x[x > 0])))})
#   # apply(adj_mat, MARGIN = 1, FUN = function(x){ sum(x > 0)}) # checking to see how many species in each country
#   trait <- apply(adj_mat, MARGIN = 2, FUN = function(x) unique(x)[unique(x)> 0]) # zeros are not real values and all fruit lengths are non-zero
#   iterations = 1000
#   rand <- matrix(rep(NA, iterations * nrow(adj_mat)), nrow = nrow(adj_mat))
#   for(i in 1:iterations){
#     # print(i)
#     randMat <- randomizeMatrix( decostand(adj_mat, "pa"), null.model = "richness") 
#     # multiplies occurrences with trait values
#     randMat_trait <- mapply("*", as.data.frame(randMat), unlist(trait)) 
#     rand[,i] <- apply(randMat_trait, MARGIN = 1, FUN = function(x) {diff(range( log(x[x>0]) ) )})
#   }
#   randMean <- apply(rand, MARGIN = 1, FUN = mean)
#   randSd <- apply(rand, MARGIN = 1, FUN = sd)
#   z <- (obs-randMean) / randSd
#   res <- data.frame( z, "LEVEL_3_CO" = names(z))
#   names(res)[1] <- paste0(prefix, "range_z")
#   return(res)
# }
# 
# palm_range_z_realm <- ddply(.data = palm_occ_trait2,
#                             .variables = .(THREEREALM),
#                             .fun = ses_range,
#                             prefix = "fruit_realm", value.var = "AverageFruitLength_cm_filled")
# 
# palm_range_z_global <- ses_range(palm_occ_trait2, prefix = "fruit_global", value.var = "AverageFruitLength_cm_filled")
# palm_range_z_res <- merge(palm_range_z_realm, palm_range_z_global, by = "LEVEL_3_CO")
# tdwg_meanFruit <- merge(tdwg_meanFruit, palm_range_z_res, by = "LEVEL_3_CO", all.x = TRUE)

# Calculate mean and median body sizes of present natural mammal assemblages
mammal_presnat_occ_trait <- merge(mammal_presnat_comb_occ, phylacine_trait, by.x = "SpecName", by.y = "Binomial.1.2", all.x = TRUE)
tdwg_presnat_meanBodySize <-
  ddply(.data = subset(mammal_presnat_occ_trait, LEVEL_3_CO %in% mammal_palm_intersect),
        .variables = .(LEVEL_3_CO),
        .fun = summarize,
        presNat_meanBodySize = mean(Mass.g, na.rm = T),
        presNat_medianBodySize = median(Mass.g, na.rm = T),
        presNat_maxBodySize = max(Mass.g, na.rm = T),
        presNat_max95BodySize = quantile(Mass.g, probs = 0.95, na.rm = T ),
        presNat_minBodySize = min(Mass.g, na.rm = T),
        presNat_dispBodySize = mean(log(Mass.g) - log(presNat_medianBodySize)) ,
        presNat_rangeBodySize = log(presNat_maxBodySize) - log(presNat_minBodySize),
        presNat_sdBodySize = sd(log(Mass.g), na.rm = T),
        presNat_nSp = length(unique(SpecName)),
        presNat_meso_nSp = length(unique(SpecName[Mass.g > 10000])),
        presNat_mega_nSp = length(unique(SpecName[Mass.g > 44000])))
# mammal_presnat_occ_trait2 <- merge(x = subset(mammal_presnat_occ_trait, LEVEL_3_CO %in% mammal_palm_intersect), y = tdwg_env[c("LEVEL_3_CO", "REALM_LONG", "THREEREALM")], by = "LEVEL_3_CO")
# mammal_presnat_occ_trait2 <- subset(mammal_presnat_occ_trait2, REALM_LONG %in% c("Neotropics", "Afrotropics", "IndoMalay", "Australasia"))
# 
# mam_range_z_realm <- ddply(.data = mammal_presnat_occ_trait2,
#                            .variables = .(THREEREALM),
#                            .fun = ses_range,
#                            prefix = "pnat_realm", value.var = "Mass.g", .progress = "text")
# mam_range_z_global <- ses_range(mammal_presnat_occ_trait2, prefix = "pnat_global", value.var = "Mass.g")
# mam_range_z_res <- merge(mam_range_z_realm, mam_range_z_global, by = "LEVEL_3_CO")
# tdwg_presnat_meanBodySize <- merge(tdwg_presnat_meanBodySize, mam_range_z_res, by = "LEVEL_3_CO", all.x = TRUE)

# Current the mean and median body sizes of current mammal assemblages
mammal_curr_occ_trait <- merge(mammal_curr_comb_occ, phylacine_trait, by.x = "SpecName", by.y = "Binomial.1.2", all.x = TRUE)
tdwg_curr_meanBodySize <- 
  ddply(.data = subset(mammal_curr_occ_trait, LEVEL_3_CO %in% mammal_palm_intersect),
        .variables = .(LEVEL_3_CO),
        .fun = summarize,
        curr_meanBodySize = mean(Mass.g, na.rm = T),
        curr_medianBodySize = median(Mass.g, na.rm = T),
        curr_maxBodySize = max(Mass.g, na.rm = T),
        curr_max95BodySize = quantile(Mass.g, probs = 0.95, na.rm = T),
        curr_minBodySize = min(Mass.g, na.rm = T),
        curr_rangeBodySize = log(curr_maxBodySize) - log(curr_minBodySize),
        curr_dispBodySize = mean(log(Mass.g) - log(curr_medianBodySize)) ,
        curr_sdBodySize = sd(log(Mass.g), na.rm = T),
        curr_nSp = length(unique(SpecName)),
        curr_meso_nSp = length(unique(SpecName[Mass.g > 10000])),
        curr_mega_nSp = length(unique(SpecName[Mass.g > 44000])),
        futr_medianBodySize = median(Mass.g[!IUCN.Status.1.2 %in% c("CR","EW","EN","EX")], na.rm = T),
        futr_maxBodySize = quantile(Mass.g[!IUCN.Status.1.2 %in% c("CR", "EW", "EN", "EX")], probs = 0.95, na.rm = T))

# mammal_curr_occ_trait2 <- merge(x = subset(mammal_curr_occ_trait, LEVEL_3_CO %in% mammal_palm_intersect), y = tdwg_env[c("LEVEL_3_CO", "REALM_LONG", "THREEREALM")], by = "LEVEL_3_CO")
# mammal_curr_occ_trait2 <- subset(mammal_curr_occ_trait2, REALM_LONG %in% c("Neotropics", "Afrotropics", "IndoMalay", "Australasia"))
# 
# mam_curr_range_z_realm <- ddply(.data = mammal_curr_occ_trait2,
#                            .variables = .(THREEREALM),
#                            .fun = ses_range,
#                            prefix = "curr_realm", value.var = "Mass.g", .progress = "text")
# mam_curr_range_z_global <- ses_range(mammal_curr_occ_trait2, prefix = "curr_global", value.var = "Mass.g")
# mam_curr_range_z_res <- merge(mam_curr_range_z_realm, mam_curr_range_z_global, by = "LEVEL_3_CO")
# tdwg_curr_meanBodySize <- merge(tdwg_curr_meanBodySize, mam_curr_range_z_res, by = "LEVEL_3_CO", all.x = TRUE)

# NOTE: Pteropus niger (Mascarene fruit bat) is found on both Mauritius and Reunion but is only on Mauritius for the current dataset, as a result, there are no mammals on Reunion in the current case (Pteropus) but 2 in the present-natural dataset

# Merge present natural and current mammal assemblage summary statistics and fruit statistics
tdwg_mammal_all <- merge(tdwg_presnat_meanBodySize, tdwg_curr_meanBodySize, by = "LEVEL_3_CO", all = TRUE)
tdwg_res <- merge(tdwg_meanFruit, tdwg_mammal_all, by = "LEVEL_3_CO")

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

# Remove remote oceanic islands and countries with no environmental data
# remoteIslList <- subset(tdwg_final, ISISLAND == 1 & DistToCont_km > 2000) # List of islands to remove
# tdwg_final2 <- subset(tdwg_final, (!LEVEL_3_CO %in% union(remoteIslList$LEVEL_3_CO, c("MRQ", "TUA", "VNA"))) & Geology_3categ %in% c("continental", "mainland", "volcanic"))

tdwg_final2 <- tdwg_final
#tdwg_final2 <- subset(tdwg_final2, PALMSR>0)
# Islands that are greater than 2000 km from closest continent
# MRQ, TUA, VNA (Venuzuelan antilles removed due to lack of environmental information)
# In total 13 removed (not including REU which was removed earlier)

# Global level PCA
env.var <- paste0("bio", c(12, 15, 17, 1, 11, 4), "_mean")
  #paste0("bio", 1:19, "_mean")
#c("PREC_Sum", "PREC_CV", "P_drie_quart", "Tmean_mean", "T_cold_quart", "Temp_SD")
globalPCA <- prcomp(x = tdwg_final2[env.var], scale. = TRUE, center = TRUE)
cumsum(globalPCA$sdev) /sum(globalPCA$sdev)
tdwg_final2$globalPC1 <- globalPCA$x[,1]
tdwg_final2$globalPC2 <- globalPCA$x[,2]
tdwg_final2$globalPC3 <- globalPCA$x[,3]

# Regional PCA
tdwg_final_nw <- subset(tdwg_final2, REALM_LONG == "Neotropics")
tdwg_final_oww <- subset(tdwg_final2, REALM_LONG == "Afrotropics")
tdwg_final_owe <- subset(tdwg_final2, REALM_LONG %in% c("IndoMalay", "Australasia"))

nwPCA <- prcomp(x= tdwg_final_nw[env.var], scale. = TRUE, center = TRUE)
cumsum(nwPCA$sdev) / sum(nwPCA$sdev)
owwPCA <- prcomp(x= tdwg_final_oww[env.var], scale. = TRUE, center = TRUE)
cumsum(owwPCA$sdev) / sum(owwPCA$sdev)
owePCA <- prcomp(x= tdwg_final_owe[env.var], scale. = TRUE, center = TRUE)
cumsum(owePCA$sdev) / sum(owePCA$sdev)

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
