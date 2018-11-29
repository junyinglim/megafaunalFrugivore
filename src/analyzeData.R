# Analyze Data

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
library(plyr); library(dplyr)
library(ggplot2); library(sp); library(rgdal); library(viridis); library(gridExtra)
library(reshape2)

## Import palm data ========================
palm_occ <- read.csv(file.path(data.dir, "palms_in_tdwg3.csv"))
palm_trait <- read.csv(file.path(data.dir, "PalmTraits_10.csv"))
palm_trait$SpecName <- gsub(palm_trait$SpecName, pattern= " ", replacement = "_")

## Import tdwg data ========================
tdwg_env <- read.csv(file.path(data.dir, "TDWG_Environment_AllData_2014Dec.csv"))

## Import diet data ========================
# Import MammalDiet v1.0 dataset
mammalDiet <- read.delim(file.path(frug.dir, "mammalDiet_v1.0.txt"))
mammalDiet$SpecName <- paste(mammalDiet$Genus, mammalDiet$Species, sep = "_")

# Import phylacine dataset
phylacine.dir <- file.path(frug.dir, "PHYLACINE")
phylacine_trait <- read.csv(file.path(phylacine.dir, "Trait_data.csv"))
spatial_metadata <- read.csv(file.path(phylacine.dir, "Spatial_metadata.csv"))

extinct_herb_splist <- subset(phylacine_trait, IUCN.Status.1.2 == "EP" & Diet.Plant == 100)$Binomial.1.2


# Import mammal occ datasets
mammal_presnat_occ <- readRDS(file.path(data.dir, "mammal_presNat_occ_raw.rds"))
mammal_nochange_occ <- readRDS(file.path(data.dir, "mammal_nochange_occ_raw.rds"))
mammal_curr_occ <- readRDS(file.path(data.dir, "mammal_current_occ_raw.rds"))
mammal_flagged_occ <- read.csv(file.path(data.dir, "flaggedPresentNaturalRanges.csv"))

mammal_presnat_occ2 <- rbind(mammal_presnat_occ[,-1], mammal_nochange_occ[,-1])

# Remove false positives
falsePos <- subset(mammal_flagged_occ, Keep == 1)[c("LEVEL_3_CO", "SpecName")]

mammal_presnat_occ2$falsePos <- ifelse(is.na(match(paste0(mammal_presnat_occ2$LEVEL_3_CO, mammal_presnat_occ2$SpecName),paste0(falsePos$LEVEL_3_CO, falsePos$SpecName))),"No", "Yes")
sum(mammal_presnat_occ2$falsePos == "Yes")
mammal_presnat_comb_occ <- subset(mammal_presnat_occ2, falsePos == "No")

mammal_curr_occ2 <- rbind(mammal_curr_occ[,-1], mammal_nochange_occ[,-1])
mammal_curr_occ2$falsePos <- ifelse(is.na(match(paste0(mammal_curr_occ2$LEVEL_3_CO, mammal_curr_occ2$SpecName),paste0(falsePos$LEVEL_3_CO, falsePos$SpecName))),"No", "Yes")
sum(mammal_curr_occ2$falsePos == "Yes")
mammal_curr_comb_occ <- subset(mammal_curr_occ2, falsePos == "No")

# Only include areas that contain palms
palm_tdwg_list <- unique(palm_occ$Area_code_L3)
mammal_presnat_comb_occ <- subset(mammal_presnat_comb_occ, LEVEL_3_CO %in% palm_tdwg_list)
mammal_curr_comb_occ <- subset(mammal_curr_comb_occ, LEVEL_3_CO %in% palm_tdwg_list)

length(unique(mammal_curr_comb_occ$SpecName))
sum(!unique(mammal_presnat_comb_occ$SpecName) %in% extinct_herb_splist)

## Genus-level gap filling ======================== 
palm_trait_subset <- palm_trait[c("PalmTribe","accGenus", "SpecName","AverageFruitLength_cm")]

gapFill <- function(x){
  if( sum(is.na(x$AverageFruitLength_cm)) == 0 ){
    return(x)  
  } else {
    genusMean <- mean(x$AverageFruitLength_cm, na.rm = TRUE)
    x$AverageFruitLength_cm[is.na(x$AverageFruitLength_cm)] <- genusMean
    return(x)
  }
}

palm_trait_gapfill  <- ddply(.data = palm_trait_subset,
                             .variables = .(accGenus),
                             .fun = gapFill)
names(palm_trait_gapfill)[names(palm_trait_gapfill) == "AverageFruitLength_cm"] <- "AverageFruitLength_cm_filled"
# Only one gap not filled this way "Butyagrus nabonnandii". However, species is a sterile hybrid between queen palm (Syagrus romanzoffiana) and pindo palm (Butia odora)
palm_trait <- merge(palm_trait, palm_trait_gapfill[c("SpecName", "AverageFruitLength_cm_filled")])

## Calculate tdwg level metrics ======================== 
palm_occ_trait <- merge(palm_occ, palm_trait[c("SpecName", "AverageFruitLength_cm", "AverageFruitLength_cm_filled")])
palm_occ_trait <- palm_occ_trait[!palm_occ_trait$SpecName %in% c("Lodoicea_maldivica", "Cocos_nucifera"),]

mammal_countrylist <- Reduce(union, list(unique(mammal_presnat_occ$LEVEL_3_CO),
                                         unique(mammal_nochange_occ$LEVEL_3_CO),
                                         unique(mammal_curr_occ$LEVEL_3_CO)))

mammal_palm_intersect <- intersect(unique(palm_occ_trait$Area_code_L3), mammal_countrylist) # 190 countries in total

tdwg_meanFruit <- ddply(.data = subset(palm_occ_trait, Area_code_L3 %in% mammal_palm_intersect),
                        .variables = .(Area_code_L3),
                        .fun = summarise,
                        meanFruitLength = mean(AverageFruitLength_cm, na.rm = TRUE),
                        meanFruitLengthFilled = mean(AverageFruitLength_cm_filled, na.rm = TRUE),
                        palm_nSp = length(AverageFruitLength_cm))
names(tdwg_meanFruit)[names(tdwg_meanFruit) == "Area_code_L3"] <- "LEVEL_3_CO"

mammal_presnat_occ_trait <- merge(mammal_presnat_comb_occ, phylacine_trait, by.x = "SpecName", by.y = "Binomial.1.2", all.x = TRUE)
tdwg_presnat_meanBodySize <-
  ddply(.data = subset(mammal_presnat_occ_trait, LEVEL_3_CO %in% mammal_palm_intersect),
        .variables = .(LEVEL_3_CO),
        .fun = summarize,
        presNat_meanBodySize = mean(Mass.g, na.rm = TRUE),
        presNat_nSp = length(unique(SpecName)))

megaherbivores_splist <- subset(phylacine_trait, Mass.g >= 1000000)$Binomial.1.2
tdwg_presnat_mega_meanBodySize <- 
  ddply(.data = subset(mammal_presnat_occ_trait, SpecName %in% megaherbivores_splist & LEVEL_3_CO %in% mammal_palm_intersect),
        .variables = .(LEVEL_3_CO),
        .fun = summarize,
                presNat_mega_meanBodySize = mean(Mass.g, na.rm = TRUE),
                presNat_mega_nSp = length(unique(SpecName)))


mammal_curr_occ_trait <- merge(mammal_curr_comb_occ, phylacine_trait, by.x = "SpecName", by.y = "Binomial.1.2", all.x = TRUE)
tdwg_curr_meanBodySize <- 
  ddply(.data = subset(mammal_curr_occ_trait, LEVEL_3_CO %in% mammal_palm_intersect),
        .variables = .(LEVEL_3_CO),
        .fun = summarize,
        curr_meanBodySize = mean(Mass.g, na.rm = TRUE),
        curr_nSp = length(unique(SpecName)))

tdwg_curr_mega_meanBodySize <- 
  ddply(.data = subset(mammal_curr_occ_trait, SpecName %in% megaherbivores_splist & LEVEL_3_CO %in% mammal_palm_intersect),
        .variables = .(LEVEL_3_CO),
        .fun = summarize,
        curr_mega_meanBodySize = mean(Mass.g, na.rm = TRUE),
        curr_mega_nSp = length(unique(SpecName)))

tdwg_mammal_all <- Reduce(function(x, y) merge(x, y, by = "LEVEL_3_CO", all = T), list(tdwg_presnat_meanBodySize, tdwg_presnat_mega_meanBodySize, tdwg_curr_meanBodySize, tdwg_curr_mega_meanBodySize))
#tdwg_meanBodySize <- merge(tdwg_curr_meanBodySize, tdwg_presnat_meanBodySize, all = TRUE)
tdwg_res <- merge(tdwg_meanFruit, tdwg_mammal_all)#, all = TRUE)

write.csv(tdwg_res, file.path(res.dir, "tdwg_res.csv"), row.names = FALSE)

## Plot mammal body size and fruit sizes ========================
# Import tdwg shape files
tdwg_shape <- readOGR(file.path(rawdata.dir, "TDWG", "level3", "level3.shp"))
tdwg_shape@data$id <- rownames(tdwg_shape@data)
tdwg_shape_df <- fortify(tdwg_shape, region = "id")
tdwg_shape_df <- merge(tdwg_shape_df, tdwg_shape@data)

tdwg_plot_df <- merge(tdwg_res, tdwg_env, by = "LEVEL_3_CO")

fruitLengthPlot  <- ggplot() +
  geom_polygon(aes(y = lat, x = long, group = group), data = subset(tdwg_shape_df, !LEVEL_NAME == "Antarctica")) +
  geom_point(aes(x = LONG, y = LAT, fill = meanFruitLength, size = meanFruitLength), colour = "white", shape = 21, data = tdwg_plot_df, alpha = 0.7) +
  scale_fill_viridis() +
  theme(legend.position = "bottom")
ggsave(filename = file.path(fig.dir, "fruitSizePlot.pdf"), fruitLengthPlot, width = 14, height = 7)

gappedFruitLengthPlot  <- ggplot() +
  geom_polygon(aes(y = lat, x = long, group = group), data = subset(tdwg_shape_df, !LEVEL_NAME == "Antarctica")) +
  geom_point(aes(x = LONG, y = LAT, fill = meanFruitLengthFilled, size = meanFruitLengthFilled), colour = "white", shape = 21, data = tdwg_plot_df, alpha = 0.7) +
  scale_fill_viridis() + 
  theme(legend.position = "bottom")
ggsave(filename = file.path(fig.dir, "fruitSizeFilledPlot.pdf"), gappedFruitLengthPlot, width = 14, height = 7)

FruitLengthFilledvsNotPlot <- ggplot() + 
  geom_point(aes(y = log(meanFruitLength), x = log(meanFruitLengthFilled)), data = tdwg_plot_df) +
  geom_abline(aes(intercept = 0, slope = 1))
temp <- resid(lm(log(meanFruitLength) ~ log(meanFruitLengthFilled), data=tdwg_plot_df))

FruitLengthFilledvsNotPlotResid <- ggplot(data.frame("resid" = temp[which(tdwg_plot_df$meanFruitLength != tdwg_plot_df$meanFruitLengthFilled)])) +
  geom_histogram(aes(x = resid), bins = 15)

FruitLengthSensitivity<- grid.arrange(FruitLengthFilledvsNotPlot, FruitLengthFilledvsNotPlotResid, ncol =2)
ggsave(filename = file.path(fig.dir, "FruitLengthSensitivity.pdf"), FruitLengthSensitivity)


target.col <- c("LAT", "LONG", "LEVEL_3_CO", "curr_meanBodySize", "curr_nSp", "presNat_meanBodySize","presNat_nSp", "curr_mega_nSp", "presNat_mega_nSp", "presNat_mega_meanBodySize", "curr_mega_meanBodySize")
tdwg_plot_melt <- melt(tdwg_plot_df[,target.col], id.vars = c("LAT", "LONG", "LEVEL_3_CO"),
                      measure.vars = c("curr_meanBodySize",
                                       "presNat_meanBodySize",
                                       "curr_mega_meanBodySize",
                                       "presNat_mega_meanBodySize",
                                       "curr_nSp", "presNat_nSp",
                                       "curr_mega_nSp", "presNat_mega_nSp"))

mammalbodySizePlot <- ggplot() +
  geom_polygon(aes(y = lat, x = long, group = group), data = subset(tdwg_shape_df, !LEVEL_3_CO == "ANT")) +
  geom_point(aes(x = LONG, y = LAT, fill = value, size = value), colour = "white", shape = 21, data = subset(tdwg_plot_melt, variable %in% c("curr_meanBodySize", "presNat_meanBodySize", "curr_mega_meanBodySize", "presNat_mega_meanBodySize")), alpha = 0.7) +
  facet_wrap(~variable, nrow = 2) + 
  scale_fill_viridis()

ggsave(mammalbodySizePlot, filename = file.path(fig.dir, "meanBodySizeCombined.pdf"), width = 20, height = 10)


mammalSpRichPlot <- ggplot() +
  geom_polygon(aes(y = lat, x = long, group = group), data = subset(tdwg_shape_df, !LEVEL_3_CO == "ANT")) +
  geom_point(aes(x = LONG, y = LAT, fill = log(value), size = log(value)), colour = "white", shape = 21, data = subset(tdwg_plot_melt, variable %in% c("curr_nSp", "presNat_nSp", "curr_mega_nSp", "presNat_mega_nSp")), alpha = 0.7) +
  facet_wrap(~variable, nrow = 2) + 
  scale_fill_viridis()

ggsave(mammalSpRichPlot, filename = file.path(fig.dir, "mammalSpRichCombined.pdf"), width = 20, height = 10)

## Exploratory plots ========================
tdwg_final <- merge(tdwg_res, tdwg_env, by = "LEVEL_3_CO",all.x = TRUE)

ggplot(aes(y = log(meanFruitLengthFilled), x = log(curr_meanBodySize)), data = tdwg_final) +
  geom_point() +
  facet_wrap(~THREEREALM) +
  geom_smooth(method = "lm")
ggplot(aes(y = log(meanFruitLengthFilled), x = log(presNat_meanBodySize)), data = tdwg_final) +
  geom_point() +
  facet_wrap(~THREEREALM) +
  geom_smooth(method = "lm") 

test <- melt(tdwg_final, id.var = c("LEVEL_3_CO", "THREEREALM", "meanFruitLengthFilled"), measure.vars = c("presNat_meanBodySize", "curr_meanBodySize"))

fruitSize_mammalBodySizePlot  <- ggplot(aes(y = log(meanFruitLengthFilled), x = log(value)), data = subset(test, !LEVEL_3_CO %in% "VNA")) +
  geom_point() +
  facet_wrap(variable~THREEREALM) +
  geom_smooth(method = "lm") 
ggsave(fruitSize_mammalBodySizePlot, file = file.path(fig.dir, "fruitVsMammalPlot.pdf"))

summary(lm(log(meanFruitLengthFilled)~ log(curr_meanBodySize), data = subset(tdwg_final, THREEREALM == "NewWorld")))
summary(lm(log(meanFruitLengthFilled)~ log(presNat_meanBodySize), data = subset(tdwg_final, THREEREALM == "NewWorld")))

mammalCurrVsPresNatBodySize <- ggplot(data = tdwg_final) + geom_point(aes(y = log(presNat_meanBodySize), x = log(curr_meanBodySize), color = THREEREALM)) + geom_abline(aes(intercept = 0, slope = 1))
ggsave(mammalCurrVsPresNatBodySize, file = file.path(fig.dir, "mammalCurrVsPresNatBodySize.pdf"))

## SEM analysis ========================
library(lavaan)
library(semTools)
library(semPlot)

# PCA of climate globally
currClimVar <- c("PREC_Sum", "PREC_CV", "P_drie_quart", "Tmean_mean", "T_cold_quart", "Temp_SD")
tdwg_env_subset <- tdwg_env[complete.cases(tdwg_env[currClimVar]), ]
pcaClim <- prcomp(tdwg_env_subset[currClimVar], scale. = TRUE, center = TRUE)
summary(pcaClim)    # first 3 axes explain ~93.8%
tdwg_env_subset <- cbind(tdwg_env_subset, pcaClim$x[,1:3])

# Merge environmental variables 
tdwg_res_sem <- merge(tdwg_res, tdwg_env_subset, by = "LEVEL_3_CO", all.x = TRUE)
write.csv(tdwg_res_sem, file.path(res.dir, "semAnalysis.csv"), row.names = FALSE)

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


