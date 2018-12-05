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
extinct_herb_all_splist <- subset(phylacine_trait, IUCN.Status.1.2 == "EP" & Diet.Plant<= 100 & Diet.Plant >= 50)$Binomial.1.2

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

# Remove tertiary frugivores and duplicate rows
tertFrug_spList <- subset(mammalDiet, Fruit == 3)$SpecName

mammal_presnat_occ4 <- subset(mammal_presnat_occ3, falsePos == "No" & (!SpecName %in% tertFrug_spList))
mammal_presnat_occ_final <- mammal_presnat_occ4[!duplicated(mammal_presnat_occ4),]

mammal_curr_occ4 <- subset(mammal_curr_occ3, falsePos == "No" & (!SpecName %in% tertFrug_spList))
mammal_curr_occ_final <- mammal_curr_occ4[!duplicated(mammal_curr_occ4),]

# Only include areas that contain palms
palm_tdwg_list <- unique(palm_occ$Area_code_L3)
mammal_presnat_comb_occ <- subset(mammal_presnat_occ_final, LEVEL_3_CO %in% palm_tdwg_list)
mammal_curr_comb_occ <- subset(mammal_curr_occ_final, LEVEL_3_CO %in% palm_tdwg_list)

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
                        megapalm_nsp = sum(AverageFruitLength_cm_filled > 4),
                        palm_nSp = length(AverageFruitLength_cm))
names(tdwg_meanFruit)[names(tdwg_meanFruit) == "Area_code_L3"] <- "LEVEL_3_CO"

mammal_presnat_occ_trait <- merge(mammal_presnat_comb_occ, phylacine_trait, by.x = "SpecName", by.y = "Binomial.1.2", all.x = TRUE)
tdwg_presnat_meanBodySize <-
  ddply(.data = subset(mammal_presnat_occ_trait, LEVEL_3_CO %in% mammal_palm_intersect),
        .variables = .(LEVEL_3_CO),
        .fun = summarize,
        presNat_meanBodySize = mean(log(Mass.g), na.rm = TRUE),
        presNat_medianBodySize = median(log(Mass.g), na.rm = TRUE),
        presNat_nSp = length(unique(SpecName)))

# megaherbivores_splist <- subset(phylacine_trait, Mass.g >= 1000000)$Binomial.1.2
# tdwg_presnat_mega_meanBodySize <- 
#   ddply(.data = subset(mammal_presnat_occ_trait, SpecName %in% megaherbivores_splist & LEVEL_3_CO %in% mammal_palm_intersect),
#         .variables = .(LEVEL_3_CO),
#         .fun = summarize,
#                 presNat_mega_meanBodySize = mean(Mass.g, na.rm = TRUE),
#                 presNat_mega_medianBodySize = median(Mass.g, na.rm = TRUE),
#                 presNat_mega_nSp = length(unique(SpecName))
#           )

# Current the mean and median body sizes of current 
mammal_curr_occ_trait <- merge(mammal_curr_comb_occ, phylacine_trait, by.x = "SpecName", by.y = "Binomial.1.2", all.x = TRUE)
tdwg_curr_meanBodySize <- 
  ddply(.data = subset(mammal_curr_occ_trait, LEVEL_3_CO %in% mammal_palm_intersect),
        .variables = .(LEVEL_3_CO),
        .fun = summarize,
        curr_meanBodySize = mean(log(Mass.g), na.rm = TRUE),
        curr_medianBodySize = median(log(Mass.g), na.rm = TRUE),
        curr_nSp = length(unique(SpecName)))

# # Calculate the mean and median body sizes of mega herbivores for current
# tdwg_curr_mega_meanBodySize <- 
#   ddply(.data = subset(mammal_curr_occ_trait, SpecName %in% megaherbivores_splist & LEVEL_3_CO %in% mammal_palm_intersect),
#         .variables = .(LEVEL_3_CO),
#         .fun = summarize,
#         curr_mega_meanBodySize = mean(Mass.g, na.rm = TRUE),
#         curr_mega_medianBodySize = median(Mass.g, na.rm= TRUE),
#         curr_mega_nSp = length(unique(SpecName)))

#tdwg_mammal_all <- Reduce(function(x, y) merge(x, y, by = "LEVEL_3_CO", all = T), list(tdwg_presnat_meanBodySize, tdwg_presnat_mega_meanBodySize, tdwg_curr_meanBodySize, tdwg_curr_mega_meanBodySize))
tdwg_mammal_all <- merge(tdwg_presnat_meanBodySize, tdwg_curr_meanBodySize, by = "LEVEL_3_CO", all = TRUE)
tdwg_res <- merge(tdwg_meanFruit, tdwg_mammal_all)#, all = TRUE)

write.csv(tdwg_res, file.path(res.dir, "tdwg_res.csv"), row.names = FALSE)

# write.csv(merge(mammal_presnat_comb_occ, phylacine_trait, by.x = "SpecName", by.y = "Binomial.1.2", all.x = TRUE), file.path(data.dir, "mammal_presnat_comb_occ.csv"), row.names = FALSE)
# write.csv(merge(mammal_curr_comb_occ, phylacine_trait, by.x = "SpecName", by.y = "Binomial.1.2", all.x = TRUE), file.path(data.dir, "mammal_curr_comb_occ.csv"), row.names = FALSE)

## Plot mammal body size and fruit sizes ========================
# Import tdwg shape files
tdwg_shape <- readOGR(file.path(rawdata.dir, "TDWG", "level3", "level3.shp"))
tdwg_shape@data$id <- rownames(tdwg_shape@data)
tdwg_shape_df <- fortify(tdwg_shape, region = "id")
tdwg_shape_df <- merge(tdwg_shape_df, tdwg_shape@data)
tdwg_plot_df <- merge(tdwg_res, tdwg_env, by = "LEVEL_3_CO")

# fruitLengthPlot  <- ggplot() +
#   geom_polygon(aes(y = lat, x = long, group = group), data = subset(tdwg_shape_df, !LEVEL_NAME == "Antarctica")) +
#   geom_point(aes(x = LONG, y = LAT, fill = meanFruitLength, size = meanFruitLength), colour = "white", shape = 21, data = tdwg_plot_df, alpha = 0.7) +
#   scale_fill_viridis() +
#   theme(legend.position = "bottom")
# ggsave(filename = file.path(fig.dir, "fruitSizePlot.pdf"), fruitLengthPlot, width = 14, height = 7)

gappedFruitLengthPlot  <- ggplot() +
  geom_polygon(aes(y = lat, x = long, group = group), data = subset(tdwg_shape_df, !LEVEL_NAME == "Antarctica")) +
  geom_point(aes(x = LONG, y = LAT, fill = meanFruitLengthFilled, size = meanFruitLengthFilled), colour = "white", shape = 21, data = tdwg_plot_df, alpha = 0.7) +
  scale_fill_viridis() + 
  theme(legend.position = "bottom")
ggsave(filename = file.path(fig.dir, "fruitSizeFilledPlot.pdf"), gappedFruitLengthPlot, width = 14, height = 7)

# FruitLengthFilledvsNotPlot <- ggplot() + 
#   geom_point(aes(y = log(meanFruitLength), x = log(meanFruitLengthFilled)), data = tdwg_plot_df) +
#   geom_abline(aes(intercept = 0, slope = 1))
#temp <- resid(lm(log(meanFruitLength) ~ log(meanFruitLengthFilled), data=tdwg_plot_df))

#FruitLengthFilledvsNotPlotResid <- ggplot(data.frame("resid" = temp[which(tdwg_plot_df$meanFruitLength != tdwg_plot_df$meanFruitLengthFilled)])) +
#  geom_histogram(aes(x = resid), bins = 15)
#FruitLengthSensitivity<- grid.arrange(FruitLengthFilledvsNotPlot, FruitLengthFilledvsNotPlotResid, ncol =2)
#ggsave(filename = file.path(fig.dir, "FruitLengthSensitivity.pdf"), FruitLengthSensitivity)


target.col <- c("LAT", "LONG", "LEVEL_3_CO", "curr_meanBodySize", "curr_nSp", "presNat_meanBodySize","presNat_nSp")#, "curr_mega_nSp", "presNat_mega_nSp", "presNat_mega_meanBodySize", "curr_mega_meanBodySize")
tdwg_plot_melt <- melt(tdwg_plot_df[,target.col], id.vars = c("LAT", "LONG", "LEVEL_3_CO"),
                      measure.vars = c("curr_meanBodySize",
                                       "presNat_meanBodySize",
                                       #"curr_mega_meanBodySize",
                                       #"presNat_mega_meanBodySize",
                                       "curr_nSp", "presNat_nSp",
                                       "curr_mega_nSp", "presNat_mega_nSp"))


mammalbodySizePlot <- ggplot() +
  geom_polygon(aes(y = lat, x = long, group = group), data = subset(tdwg_shape_df, !LEVEL_3_CO == "ANT")) +
  geom_point(aes(x = LONG, y = LAT, fill = value, size = value), colour = "white", shape = 21, data = subset(tdwg_plot_melt, variable %in% c("curr_meanBodySize", "presNat_meanBodySize")), alpha = 0.7) +
  facet_wrap(~variable, nrow = 2) + 
  scale_fill_viridis()
#"curr_mega_meanBodySize", "presNat_mega_meanBodySize"

ggsave(mammalbodySizePlot, filename = file.path(fig.dir, "meanBodySizeCombined.pdf"), width = 10, height = 10)


mammalSpRichPlot <- ggplot() +
  geom_polygon(aes(y = lat, x = long, group = group), data = subset(tdwg_shape_df, !LEVEL_3_CO == "ANT")) +
  geom_point(aes(x = LONG, y = LAT, fill = value, size = value), colour = "white", shape = 21, data = subset(tdwg_plot_melt, variable %in% c("curr_nSp", "presNat_nSp", "curr_mega_nSp", "presNat_mega_nSp")), alpha = 0.7) +
  facet_wrap(~variable, nrow = 2) + 
  scale_fill_viridis()

ggsave(mammalSpRichPlot, filename = file.path(fig.dir, "mammalSpRichCombined.pdf"), width = 20, height = 10)

## Exploratory plots ========================
tdwg_final <- merge(tdwg_res, tdwg_env, by = "LEVEL_3_CO",all.x = TRUE)

test <- melt(tdwg_final, id.var = c("LEVEL_3_CO", "THREEREALM", "Geology_3categ","meanFruitLengthFilled"), measure.vars = c("presNat_meanBodySize", "curr_meanBodySize"))

fruitSize_mammalBodySizePlot  <- ggplot(aes(y = log(meanFruitLengthFilled), x = log(value)), data = subset(test, !LEVEL_3_CO %in% c("VNA"))) +
  geom_point() +
  facet_wrap(variable~THREEREALM) +
  geom_smooth(method = "lm") 
ggsave(fruitSize_mammalBodySizePlot, file = file.path(fig.dir, "fruitVsMammalPlot.pdf"))


fruitSize_mammalBodySizePlot_noisland  <- ggplot(aes(y = log(meanFruitLengthFilled), x = log(value)), data = subset(test, Geology_3categ %in% c("continental", "mainland"))) +
  geom_point() +
  facet_wrap(variable~THREEREALM) +
  geom_smooth(method = "lm") 

ggsave(fruitSize_mammalBodySizePlot_noisland, file = file.path(fig.dir, "fruitVsMammalPlot_noisland.pdf"))



summary(lm(log(meanFruitLengthFilled)~ log(curr_meanBodySize), data = subset(tdwg_final, THREEREALM == "NewWorld")))
summary(lm(log(meanFruitLengthFilled)~ log(presNat_meanBodySize), data = subset(tdwg_final, THREEREALM == "NewWorld" & ! LEVEL_3_CO %in% c("VNA", "MRQ"))))

mammalCurrVsPresNatBodySize <- ggplot(data = tdwg_final) + geom_point(aes(y = log(presNat_meanBodySize), x = log(curr_meanBodySize), color = THREEREALM)) + geom_abline(aes(intercept = 0, slope = 1))
ggsave(mammalCurrVsPresNatBodySize, file = file.path(fig.dir, "mammalCurrVsPresNatBodySize.pdf"))


ggplot(aes(x = log(curr_meanBodySize), y = meanFruitLengthFilled), data = tdwg_final) + geom_point() + geom_smooth(method = "lm") 
ggplot(aes(x = log(presNat_meanBodySize), y = log(meanFruitLengthFilled)), data = tdwg_final) + geom_point() + geom_smooth(method = "lm") 

hist(log(tdwg_final$meanFruitLengthFilled))
hist(tdwg_final$meanFruitLengthFilled)

# Calculate the proportion of mammals that are megafaunal, and the number of species of megafaunal fruits










## ISLAND TAXA
table(tdwg_env$GeologicalOrigin)
table(tdwg_env$Geology_3categ)
ggplot(aes(y = log(meanFruitLengthFilled), x=log(DistToCont_km), colour = Geology_3categ), data = subset(tdwg_final, Geology_3categ %in% c("continental", "mainland", "volcanic"))) + geom_point() + geom_smooth(method = "lm") + geom_label_repel(aes(y = log(meanFruitLengthFilled), x = log(DistToCont_km), label = LEVEL_3_CO))

library(ggrepel)

names(tdwg_env)
ggplot(aes(y = log(meanFruitLengthFilled), x=THREEREALM, color = Geology_3categ), data = tdwg_final) + geom_boxplot()
names(tdwg_env)

library(geosphere)

mat <- distm(tdwg_final[,c('LONG','LAT')], tdwg_final[,c('LONG','LAT')], fun=distGeo)
nearestTDWGlist <- tdwg_final$LEVEL_3_CO[unlist(apply(mat, MARGIN = 1, FUN = function(x){ which(min(x[x>0]) == x) }))]
distNearestTDWGlist <- unlist(apply(mat, MARGIN = 1, FUN = function(x){ min(x[x>0]) }))
tdwg_final$nearestTDWG <- nearestTDWGlist
tdwg_final$distNearestTDWG <- distNearestTDWGlist

ggplot(aes(y = log(meanFruitLengthFilled), x=distNearestTDWG, color = Geology_3categ), data = tdwg_final)
names(tdwg_env)

