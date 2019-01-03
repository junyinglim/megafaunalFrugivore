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
library(ggplot2); library(sp); library(rgdal); library(viridis); library(gridExtra)
library(car); library(MuMIn)
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
                        medianFruitLengthFilled = median(AverageFruitLength_cm_filled, na.rm = TRUE),
                        meanFruitLengthFilled = mean(AverageFruitLength_cm_filled, na.rm = TRUE),
                        megapalm_nsp = sum(AverageFruitLength_cm_filled > 4),
                        palm_nSp = length(AverageFruitLength_cm))
names(tdwg_meanFruit)[names(tdwg_meanFruit) == "Area_code_L3"] <- "LEVEL_3_CO"

mammal_presnat_occ_trait <- merge(mammal_presnat_comb_occ, phylacine_trait, by.x = "SpecName", by.y = "Binomial.1.2", all.x = TRUE)
tdwg_presnat_meanBodySize <-
  ddply(.data = subset(mammal_presnat_occ_trait, LEVEL_3_CO %in% mammal_palm_intersect),
        .variables = .(LEVEL_3_CO),
        .fun = summarize,
        presNat_meanBodySize = mean(Mass.g, na.rm = TRUE),
        presNat_medianBodySize = median(Mass.g, na.rm = TRUE),
        presNat_megaHerb_nSp = length(unique(SpecName[Mass.g > 44000])),
        presNat_nSp = length(unique(SpecName)))

# Current the mean and median body sizes of current 
mammal_curr_occ_trait <- merge(mammal_curr_comb_occ, phylacine_trait, by.x = "SpecName", by.y = "Binomial.1.2", all.x = TRUE)
tdwg_curr_meanBodySize <- 
  ddply(.data = subset(mammal_curr_occ_trait, LEVEL_3_CO %in% mammal_palm_intersect),
        .variables = .(LEVEL_3_CO),
        .fun = summarize,
        curr_meanBodySize = mean(Mass.g, na.rm = TRUE),
        curr_medianBodySize = median(Mass.g, na.rm = TRUE),
        curr_megaHerb_nSp = length(unique(SpecName[Mass.g > 44000])),
        curr_nSp = length(unique(SpecName)))

tdwg_mammal_all <- merge(tdwg_presnat_meanBodySize, tdwg_curr_meanBodySize, by = "LEVEL_3_CO", all = TRUE)
tdwg_res <- merge(tdwg_meanFruit, tdwg_mammal_all)#, all = TRUE)

write.csv(tdwg_res, file.path(res.dir, "tdwg_res.csv"), row.names = FALSE)

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
# 
gappedFruitLengthPlot  <- ggplot() +
  geom_polygon(aes(y = lat, x = long, group = group), data = subset(tdwg_shape_df, !LEVEL_3_CO == "ANT")) +
  geom_point(aes(x = LONG, y = LAT, fill = meanFruitLengthFilled, size = meanFruitLengthFilled), colour = "white", shape = 21, data = tdwg_plot_df, alpha = 0.7) +
  scale_fill_viridis() +
  theme(legend.position = "bottom")
ggsave(filename = file.path(fig.dir, "fruitSizeFilledPlot.pdf"), gappedFruitLengthPlot, width = 14, height = 7)

# FruitLengthFilledvsNotPlot <- ggplot() +
#   geom_point(aes(y = log(meanFruitLength), x = log(meanFruitLengthFilled)), data = tdwg_plot_df) +
#   geom_abline(aes(intercept = 0, slope = 1))
# temp <- resid(lm(log(meanFruitLength) ~ log(meanFruitLengthFilled), data=tdwg_plot_df))
# 
# FruitLengthFilledvsNotPlotResid <- ggplot(data.frame("resid" = temp[which(tdwg_plot_df$meanFruitLength != tdwg_plot_df$meanFruitLengthFilled)])) +
#  geom_histogram(aes(x = resid), bins = 15)
# FruitLengthSensitivity<- grid.arrange(FruitLengthFilledvsNotPlot, FruitLengthFilledvsNotPlotResid, ncol =2)
# ggsave(filename = file.path(fig.dir, "FruitLengthSensitivity.pdf"), FruitLengthSensitivity)


target.col <- c("LAT", "LONG", "LEVEL_3_CO", "curr_meanBodySize", "curr_nSp", "presNat_meanBodySize","presNat_nSp","curr_medianBodySize","presNat_medianBodySize")#, "curr_mega_nSp", "presNat_mega_nSp", "presNat_mega_meanBodySize", "curr_mega_meanBodySize")
tdwg_plot_melt <- melt(tdwg_plot_df[,target.col],
                       id.vars = c("LAT", "LONG", "LEVEL_3_CO"),
                       measure.vars = c("curr_meanBodySize",
                                       "presNat_meanBodySize",
                                       "curr_medianBodySize",
                                       "presNat_medianBodySize",
                                       "curr_nSp", "presNat_nSp"))

tdwg_plot_melt$variable <- factor(tdwg_plot_melt$variable, levels = c("curr_medianBodySize", "presNat_medianBodySize"), labels = c("Current", "Present-natural"))

#POSTER <- merge(x=tdwg_shape_df, y=tdwg_plot_df, by= "LEVEL_3_CO", all.x = TRUE)
#POSTER2 <- POSTER[order(POSTER$order),]
q_probs = seq(0.0, 1.0, 0.1)
q_labs <- levels(cut(log(tdwg_plot_melt$value), quantile(log(tdwg_plot_melt$value), probs = q_probs, na.rm = T)))
q_labs <- gsub(q_labs, pattern = "\\(|\\]", replacement = "")
q_labs <- gsub(q_labs, pattern = ",", replacement = "-")
#q_labs <- c("16.9 - 40.0", "40.0 - 64.1", "64.1 - 114.4", "114.4 - 254.7", "254.7 - 652.0", "652.0 - 3,498.2", "3,498.2 - 20,537.3", "20,537.3 - 120,571.7", "120,571.7 - 1,329,083")
q_labs <- c("0 - 17", "17 - 40", "40 - 64", "64 - 114", "114 - 254", "254 - 652", "652 - 3,498", "3,498 - 20,537", "20,537 - 120,572", "120,572 - 1,329,083")

tdwg_plot_melt$logvalue <- log(tdwg_plot_melt$value)

curr_mammalbodySizePlot <- ggplot() +
  geom_polygon(aes(y = lat, x = long, group = group), data = subset(tdwg_shape_df, !LEVEL_3_CO == "ANT")) +
  geom_point(aes(x = LONG, y = LAT, fill = cut(logvalue, quantile(logvalue, probs = q_probs), labels = q_labs), size = logvalue), colour = "white", shape = 21, data = subset(tdwg_plot_melt, variable %in% c("Current", "Present-natural") & !is.na(logvalue)), alpha = 0.9) +
  facet_wrap(~variable, nrow = 2, scales = "free")+
  theme(panel.background = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), legend.position = "bottom",legend.justification = "center", strip.background = element_blank(), strip.text = element_text(size = 18)) +
  guides(size=FALSE, fill = guide_legend(override.aes = list(size=10))) +
  scale_fill_manual(values = colorRamps::matlab.like(n = 10), name = "Median\nmammal body size (g)", breaks = q_labs)
ggsave(curr_mammalbodySizePlot, filename = file.path(fig.dir, "currmedianLogBodySizeCombined.pdf"), width = 10, height = 10)


library(RColorBrewer)
  scale_fill_viridis(limits = c(2,12),
                     name = "Log median\nmammal body size (g)", option = "A")
  #scale_fill_gradientn(limits = c(1,12), name = "Log median\nmammal body size (g)",
  #                     colors = wes_palette("Zissou1", 100, type = "continuous"))

# presnat_mammalbodySizePlot <- ggplot() +
#   geom_polygon(aes(y = lat, x = long, group = group), data = subset(POSTER2, !LEVEL_3_CO == "ANT")) +
#   geom_point(aes(x = LONG, y = LAT, fill = log(value), size = log(value)), colour = "white", shape = 21, data = subset(tdwg_plot_melt, variable %in% c("Current", "Present-natural")), alpha = 0.7) +
#   #facet_wrap(~variable, nrow = 2) + 
#   theme(panel.background = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), legend.position = "bottom",legend.justification = "center") +
#   guides(size=FALSE)+
#   scale_fill_viridis(limits = c(2,12), name = "Log median\nmammal body size (g)")
  #scale_fill_gradientn(limits = c(1,12), name = "Log median\nmammal body size (g)",
  #                     colors = wes_palette("Zissou1", 100, type = "continuous"))

#ggsave(mammalbodySizePlot, filename = file.path(fig.dir, "medianLogBodySizeCombined.pdf"), width = 10, height = 10)

#ggsave(presnat_mammalbodySizePlot, filename = file.path(fig.dir, "presnatmedianLogBodySizeCombined.pdf"), width = 10, height = 5)


mammalSpRichPlot <- ggplot() +
  geom_polygon(aes(y = lat, x = long, group = group), data = subset(tdwg_shape_df, !LEVEL_3_CO == "ANT")) +
  geom_point(aes(x = LONG, y = LAT, fill = value, size = value), colour = "white", shape = 21, data = subset(tdwg_plot_melt, variable %in% c("curr_nSp", "presNat_nSp")), alpha = 0.7) +
  facet_wrap(~variable, nrow = 2) + 
  scale_fill_viridis()

ggsave(mammalSpRichPlot, filename = file.path(fig.dir, "mammalSpRichCombined.pdf"), width = 10, height = 10)

## Exploratory plots ========================
tdwg_res$propMegaPalm <- tdwg_res$megapalm_nsp / tdwg_res$palm_nSp
tdwg_res$propMegaMam_curr <- tdwg_res$curr_megaHerb_nSp / tdwg_res$curr_nSp
tdwg_res$propMegaMam_presnat <- tdwg_res$presNat_megaHerb_nSp / tdwg_res$presNat_nSp
tdwg_final <- merge(tdwg_res, tdwg_env, by = "LEVEL_3_CO",all.x = TRUE)
tdwg_final$deltaMedianBodySize = log(tdwg_final$presNat_medianBodySize) - log(tdwg_final$curr_medianBodySize)
tdwg_final2 <- subset(tdwg_final, !LEVEL_3_CO %in% c("MRQ", "TUA", "VNA"))

# Global level PCA
env.var <- c("PREC_Sum", "PREC_CV", "P_drie_quart", "Tmean_mean", "T_cold_quart", "Temp_SD")
globalPCA <- prcomp(x = tdwg_final2[env.var], scale. = TRUE, center = TRUE)
tdwg_final2$PC1 <- globalPCA$x[,1]
tdwg_final2$PC2 <- globalPCA$x[,2]
tdwg_final2$PC3 <- globalPCA$x[,3]

library(MuMIn)
source(src.dir, "ggmodavg.R")

# Current + Global -------
target.col <- c("PC1", "PC2","PC3", "meanFruitLengthFilled", "curr_medianBodySize", "NPP_mean", "PC1", "PC2", "PC3", "ensLGM_Tano", "ensLGM_Pano", "THREEREALM")
tdwg_final3 <- tdwg_final2[complete.cases(tdwg_final2[target.col]),]
global_curr_mod <- lm(log(meanFruitLengthFilled) ~ scale(log(curr_medianBodySize)) + scale(NPP_mean) + scale(PC1) + scale(PC2) + scale(PC3) + scale(ensLGM_Tano) + scale(ensLGM_Pano), data = tdwg_final3, na.action = "na.fail")
global_curr_mod_dr <- dredge(global_curr_mod)
global_curr_mod_avg <- model.avg(global_curr_mod_dr)
ggsave(plotRelImportance(global_curr_mod_avg), filename = file.path(fig.dir, "global_curr_relimpt_plot.pdf"))

presid <- resid(global_curr_mod, type = "partial")
z <- data.frame(presid = presid[,1], var = scale(log(tdwg_final3$curr_medianBodySize))) 
global_curr_presid_plot <- ggplot(aes(y = presid, x =var),data = z)+ geom_point(shape = 1) + geom_smooth(method = "lm") + labs(y = "Mean fruit size\n(Partial residuals)", x = "Log median body size (Current)") + theme(panel.background = element_blank())
ggsave(global_curr_presid_plot, filename = file.path(fig.dir, "global_curr_presid_plot.pdf"))

# Current + New World -------
tdwg_final3_nw <- subset(tdwg_final3, THREEREALM == "NewWorld" & Geology_3categ %in% c("mainland", "continental"))
nw_curr_mod <- lm(log(meanFruitLengthFilled) ~ scale(log(curr_medianBodySize)) + scale(NPP_mean) + scale(PC1) + scale(PC2) + scale(PC3) + scale(ensLGM_Tano) + scale(ensLGM_Pano), data = tdwg_final3_nw, na.action = "na.fail")
nw_curr_mod_dr <- dredge(nw_curr_mod)
nw_curr_mod_avg <- model.avg(nw_curr_mod_dr)
nw_curr_mod_avg_stat <- plotRelImportance(nw_curr_mod_avg)
#ggsave(, filename = file.path(fig.dir, "nw_curr_relimpt_plot.pdf"))


presid <- resid(nw_curr_mod, type = "partial")
z <- data.frame(presid = presid[,1], var = scale(log(tdwg_final3_nw$curr_medianBodySize))) 
nw_curr_presid_plot <- ggplot(aes(y = presid, x =var),data = z)+ geom_point(shape = 1) + geom_smooth(method = "lm") + labs(y = "Mean fruit size\n(Partial residuals)", x = "Log median body size (Current)") + theme(panel.background = element_blank())
ggsave(nw_curr_presid_plot, filename = file.path(fig.dir, "nw_curr_presid_plot.pdf"))

# Current, OW East -------
tdwg_final3_owe <- subset(tdwg_final3, THREEREALM == "OWEast")
owe_curr_mod <- lm(log(meanFruitLengthFilled) ~ scale(log(curr_medianBodySize)) + scale(NPP_mean) + scale(PC1) + scale(PC2) + scale(PC3) + scale(ensLGM_Tano) + scale(ensLGM_Pano), data = tdwg_final3_owe, na.action = "na.fail")
owe_curr_mod_dr <- dredge(owe_curr_mod)
owe_curr_mod_avg <- model.avg(owe_curr_mod_dr)
ggsave(plotRelImportance(owe_curr_mod_avg), filename = file.path(fig.dir, "owe_curr_relimpt_plot.pdf"))

presid <- resid(owe_curr_mod, type = "partial")
z <- data.frame(presid = presid[,1], var = scale(log(tdwg_final3_owe$curr_medianBodySize))) 
owe_curr_presid_plot <- ggplot(aes(y = presid, x =var),data = z)+ geom_point(shape = 1) + geom_smooth(method = "lm") + labs(y = "Mean fruit size\n(Partial residuals)", x = "Log median body size (Current)") + theme(panel.background = element_blank())
ggsave(owe_curr_presid_plot, filename = file.path(fig.dir, "owe_curr_presid_plot.pdf"))

# Current, OW West -------
tdwg_final3_oww <- subset(tdwg_final3, THREEREALM == "OWWest")
oww_curr_mod <- lm(log(meanFruitLengthFilled) ~ scale(log(curr_medianBodySize)) + scale(NPP_mean) + scale(PC1) + scale(PC2) + scale(PC3) + scale(ensLGM_Tano) + scale(ensLGM_Pano), data = tdwg_final3_oww, na.action = "na.fail")
oww_curr_mod_dr <- dredge(oww_curr_mod)
oww_curr_mod_avg <- model.avg(oww_curr_mod_dr)
ggsave(plotRelImportance(oww_curr_mod_avg), filename = file.path(fig.dir, "oww_curr_relimpt_plot.pdf"))

presid <- resid(oww_curr_mod, type = "partial")
z <- data.frame(presid = presid[,1], var = scale(log(tdwg_final3_oww$curr_medianBodySize))) 
oww_curr_presid_plot <- ggplot(aes(y = presid, x =var),data = z)+ geom_point(shape = 1) + geom_smooth(method = "lm") + labs(y = "Mean fruit size\n(Partial residuals)", x = "Log median body size (Current)") + theme(panel.background = element_blank())
ggsave(oww_curr_presid_plot, filename = file.path(fig.dir, "oww_curr_presid_plot.pdf"))

# Present-natural, New World -------
target.col <- c("PC1", "PC2","PC3", "meanFruitLengthFilled", "presNat_medianBodySize", "NPP_mean", "PC1", "PC2", "PC3", "ensLGM_Tano", "ensLGM_Pano", "THREEREALM")
tdwg_final4 <- tdwg_final2[complete.cases(tdwg_final2[target.col]),]
tdwg_final4_nw <- subset(tdwg_final4, THREEREALM == "NewWorld" & Geology_3categ %in% c("mainland", "continental"))
nw_pnat_mod <- lm(log(meanFruitLengthFilled) ~ scale(log(presNat_medianBodySize)) + scale(NPP_mean) + scale(PC1) + scale(PC2) + scale(PC3) + scale(ensLGM_Tano) + scale(ensLGM_Pano), data = tdwg_final4_nw, na.action = "na.fail")
nw_pnat_mod_dr <- dredge(nw_pnat_mod)
nw_pnat_mod_avg <- model.avg(nw_pnat_mod_dr)
nw_pnat_mod_avg_stat <- plotRelImportance(nw_pnat_mod_avg)
#ggsave(plotRelImportance(nw_pnat_mod_avg), filename = file.path(fig.dir, "nw_pnat_relimpt_plot.pdf"))

presid <- resid(nw_pnat_mod, type = "partial")
z2 <- data.frame(presid = presid[,1], var = scale(log(tdwg_final4_nw$presNat_medianBodySize))) 
nw_pnat_presid_plot <- ggplot(aes(y = presid, x =var),data = z2)+ geom_point(shape = 1) + geom_smooth(method = "lm") + labs(y = "Mean fruit size\n(Partial residuals)", x = "Log median body size (Current)") + theme(panel.background = element_blank())
ggsave(nw_pnat_presid_plot, filename = file.path(fig.dir, "nw_pnat_presid_plot.pdf"))

z$type <- "Current"
z2$type<- "Present-natural"
z3<- rbind(z,z2)
z4 <- ggplot(aes(y = presid, x =var,color=type),data = z3)+ geom_point(shape = 1) + geom_smooth(method = "lm") + labs(y = "Mean fruit size\n(Partial residuals)", x = "Scaled log median body size") + theme(panel.background = element_blank(), legend.justification = "center", legend.title = element_blank(), legend.position = c(0.8, 0.1)) + scale_colour_manual(values = wes_palette("Cavalcanti1",2)) 
ggsave(z4, filename = file.path(fig.dir, "presid.pdf"), width = 5, height =4.5)

nw_curr_mod_avg_stat$treatment <- "Current"
nw_pnat_mod_avg_stat$treatment <- "Present-natural"
nw_all_mod_avg_stat <- rbind(nw_curr_mod_avg_stat,nw_pnat_mod_avg_stat)

newlabel <- data.frame(old = c("scale(ensLGM_Pano)","scale(ensLGM_Tano)", "scale(log(curr_medianBodySize))", "scale(NPP_mean)", "scale(PC1)", "scale(PC2)", "scale(PC3)", "scale(log(presNat_medianBodySize))"), new = c("LGM Prec anomaly", "LGM Temp anomaly", "Median Mammal Body Size", "NPP", "PC1", "PC2", "PC3", "Median Mammal Body Size"))
nw_all_mod_avg_stat$coefficient_label <- newlabel$new[match(nw_all_mod_avg_stat$coefficient, newlabel$old)]
nw_all_mod_avg_stat$coefficient_label <- factor(nw_all_mod_avg_stat$coefficient_label, level = c("Median Mammal Body Size", "NPP", "PC1", "PC2", "PC3", "LGM Prec anomaly", "LGM Temp anomaly"))

nw_relimpt_plot <- ggplot(data = nw_all_mod_avg_stat) + geom_point(aes(y = fullAvgCoef, x= coefficient_label, size = importance, color = treatment),position=position_dodge(width=0.5)) + geom_errorbar(aes(ymin = lower2.5, ymax = upper97.5, x = coefficient_label, color = treatment),position=position_dodge(width=0.5), width = 0.2) + labs(y = "Coefficients") + theme(axis.text.x = element_text(angle = 45, hjust = 0.9), axis.title.x = element_blank()) + scale_size_continuous(breaks = seq(0,1,0.25), limits = c(0,1), name = "Variable Importance") + scale_colour_manual(name = NULL, values = wes_palette("Cavalcanti1",2))

ggsave(nw_relimpt_plot, filename = file.path(fig.dir,"nw_relimpt_plot.pdf"), width = 7, height = 4)




# Present natural, global PCA
target.col <- c("PC1", "PC2","PC3", "meanFruitLengthFilled", "curr_medianBodySize", "NPP_mean", "PC1", "PC2", "PC3", "ensLGM_Tano", "ensLGM_Pano")
tdwg_final3 <- tdwg_final2[complete.cases(tdwg_final2[target.col]),]
fullmod <- lm(log(meanFruitLengthFilled) ~ log(curr_medianBodySize) + NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tano + ensLGM_Pano, data = tdwg_final3, na.action = "na.fail")
fullmod_dr <- dredge(fullmod)
fullmod_avg <- model.avg(fullmod_dr)
summary(fullmod_avg)











fullmod_step <- step(fullmod)
crPlot(fullmod, variable = "log(curr_medianBodySize)", smooth = FALSE, main = "Global ", ylab = "Partial Residuals (log mean fruit size)", xlab = "Log Current Median Body Size" )
summary(fullmod_step)

# Present natural global, global PCA
fullmod2 <- lm(log(medianFruitLengthFilled) ~ log(presNat_medianBodySize) + NPP_mean + PC1 + PC2 + PC3 + log(AREA_KM2) + ensLGM_Tano + ensLGM_Pano, data = tdwg_final2)
crPlot(fullmod2, variable = "log(presNat_medianBodySize)", smooth = FALSE, main = "Global ", ylab = "Partial Residuals (log median fruit size)", xlab = "LogCURRENT Median Body Size")
summary(fullmod2)

# New world only PCA
tdwg_final2_NW <- subset(tdwg_final2, THREEREALM == "NewWorld")
nwPCA <- prcomp(x = tdwg_final2_NW[env.var], scale. = T, center = T)
tdwg_final2_NW$PC1 <- nwPCA$x[,1]
tdwg_final2_NW$PC2 <- nwPCA$x[,2]
tdwg_final2_NW$PC3 <- nwPCA$x[,3]

# Current NW, NW PCA 
nwfullmod1 <- lm(log(meanFruitLengthFilled) ~ log(curr_medianBodySize) + NPP_mean + PC1 + PC2 + PC3+ ensLGM_Tano + ensLGM_Pano, data = tdwg_final2_NW)
step(nwfullmod1) # not significant after step-wise deletion
crPlot(nwfullmod1, variable = "log(curr_meanBodySize)", smooth = FALSE, main = "New World only ", ylab = "Partial Residuals (log mean fruit size)", xlab = "Log Current Mean Body Size")
nwmod1 <- lm(log(meanFruitLengthFilled) ~ log(curr_meanBodySize) + NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tano+ ensLGM_Pano, data = tdwg_final2_NW, na.action = 'na.fail')
nw_curr_modavg <- model.avg(dredge(nwmod1))
summary(nw_curr_modavg)

# Present natural NW, NW PCA
nwmod2 <- lm(log(meanFruitLengthFilled) ~ log(presNat_meanBodySize) + NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tano+ ensLGM_Pano, data = tdwg_final2_NW, na.action = 'na.fail')
crPlot(nwmod2, variable = "log(presNat_meanBodySize)", smooth = FALSE, main = "New World only ", ylab = "Partial Residuals (log mean fruit size)", xlab = "Log Present-natural Mean Body Size")
summary(nwmod2)
summary(model.avg(dredge(nwmod2)))

# Current NW, Global PCA
tdwg_final2_NWsubset <- subset(tdwg_final2, THREEREALM == "NewWorld")
#scale(log(medianFruitLengthFilled)) 
nwsubset_mod1 <- lm(propMegaPalm~ scale(log(presNat_medianBodySize)) + scale(NPP_mean) + scale(PC1) + scale(PC2) + scale(PC3)+ scale(ensLGM_Tano) + scale(ensLGM_Pano), data = subset(tdwg_final2, propMegaPalm >0 & THREEREALM == "NewWorld"))
summary(nwsubset_mod1)
crPlot(nwsubset_mod1, variable = "scale(log(presNat_medianBodySize))", smooth = FALSE, main = "New World only ", ylab = "Partial Residuals (log mean fruit size)", xlab = "Log Current Mean Body Size")

nwsubset_mod2 <- lm(scale(log(meanFruitLengthFilled)) ~ scale(log(presNat_medianBodySize)) + scale(NPP_mean) + scale(PC1) + scale(PC2) + scale(PC3)+ scale(ensLGM_Tano) + scale(ensLGM_Pano), data = tdwg_final2_NWsubset, na.action = "na.fail")
summary(nwsubset_mod2)
crPlot(nwsubset_mod1, variable = "log(curr_medianBodySize)", smooth = FALSE, main = "New World only ", ylab = "Partial Residuals (log mean fruit size)", xlab = "Log Current Mean Body Size")

summary(model.avg(dredge(nwsubset_mod2)))
confint(model.avg(dredge(nwsubset_mod2)))

nwsubset_mod2 <- lm(log(meanFruitLengthFilled) ~ log(presNat_meanBodySize) + NPP_mean + PC1 + PC2 + PC3+ ensLGM_Tano + ensLGM_Pano, data = tdwg_final2_NWsubset)
summary(nwsubset_mod2)
crPlot(nwsubset_mod2, variable = "log(presNat_meanBodySize)", smooth = FALSE, main = "New World only ", ylab = "Partial Residuals (log mean fruit size)", xlab = "Log Present-natural Mean Body Size")

# Change in median size
deltaSize_mod <- lm(log(meanFruitLengthFilled) ~ deltaMedianBodySize + NPP_mean + PC1 + PC2 + PC3+ ensLGM_Tano + ensLGM_Pano, data = tdwg_final2)
crPlot(deltaSize_mod, variable = "deltaMedianBodySize", smooth = FALSE, main = "New World only ", ylab = "Partial Residuals (log mean fruit size)", xlab = "Log Present-natural Mean Body Size")
summary(deltaSize_mod)

plot(deltaMedianBodySize ~ log(curr_medianBodySize), data = tdwg_final2)
plot(meanFruitLengthFilled ~ curr_meanBodySize, data = tdwg_final2)




test <- melt(tdwg_final, id.var = c("LEVEL_3_CO", "THREEREALM", "Geology_3categ","meanFruitLengthFilled", "megapalm_nsp","propMegaPalm", "AREA_KM2"), measure.vars = c("presNat_meanBodySize", "curr_meanBodySize", "presNat_megaHerb_nSp", "curr_megaHerb_nSp" , "propMegaMam_curr", "propMegaMam_presnat"))

ggplot(aes(y = meanFruitLengthFilled, x= presNat_meanBod)) + geom_point()




fruitSize_mammalBodySizePlot  <- ggplot(aes(y = propMegaPalm, x = value), data = subset(test, variable %in% c("propMegaMam_curr", "propMegaMam_presnat") & !LEVEL_3_CO %in% c("VNA", "MRQ") & Geology_3categ %in% c("continental", "mainland") & propMegaPalm > 0)) +
  geom_point() +
  #facet_wrap(~variable) +
  facet_wrap(variable~THREEREALM) +
  geom_smooth(method = "glm")
ggsave(fruitSize_mammalBodySizePlot, file = file.path(fig.dir, "fruitVsMammalPlot.pdf"))

summary(glm(megapalm_nsp ~ log(AREA_KM2), data = test, family = "poisson"))

tdwg_final$mamSdiff <- with(presNat_nSp-curr_nSp, data = tdwg_final)
tdwg_final$mamMdiff <- with(presNat_meanBodySize-curr_meanBodySize, data = tdwg_final)

summary(lm(meanFruitLengthFilled ~ curr_nSp, data = subset(tdwg_final, THREEREALM == "NewWorld")))

library(MuMIn)
tdwg_final_subset <- na.omit(tdwg_final[c("propMegaPalm", "megapalm_nsp", "meanFruitLengthFilled", "curr_meanBodySize", "presNat_meanBodySize", "NPP_mean", "PREC_Sum", "Tmean_mean", "soilcount", "Geology_3categ", "alt_range", "AREA_KM2", "mamMdiff", "mamSdiff")])
#tdwg_final_subset2 <- subset(tdwg_final_subset, Geology_3categ %in% c("mainland", "continental")) # 154 datapoints
fullmod <- lm(log(meanFruitLengthFilled) ~ curr_meanBodySize + presNat_meanBodySize + NPP_mean + PREC_Sum + Tmean_mean + soilcount + alt_range + mamMdiff + mamSdiff, data = tdwg_final_subset2, na.action = "na.fail")
moddredge <- dredge(fullmod, beta = "sd")
modavg <- model.avg(moddredge)
confint(modavg)
summary(modavg)


fullmod <- lm(propMegaPalm ~ curr_meanBodySize + presNat_meanBodySize + NPP_mean + PREC_Sum + Tmean_mean + soilcount + mamMdiff + mamSdiff, data = tdwg_final_subset, na.action = "na.fail")
moddredge <- dredge(fullmod, beta = "sd")
modavg <- model.avg(moddredge)
confint(modavg)
summary(modavg)

fullmod <- lm(log(megapalm_nsp) ~ curr_meanBodySize + presNat_meanBodySize + NPP_mean + PREC_Sum + Tmean_mean + soilcount + log(AREA_KM2) + mamMdiff + mamSdiff, data = subset(tdwg_final_subset2, megapalm_nsp >0), na.action = "na.fail", family = "poisson")
moddredge <- dredge(fullmod, beta = "sd")
modavg <- model.avg(moddredge)
confint(modavg)
summary(modavg)
ggplot() + geom_point(aes(y = log(meanFruitLengthFilled), x = mamMdiff, color = THREEREALM), data = tdwg_final)


hist(resid(fullmod))
plot(resid(fullmod)~ fitted(fullmod))
abline(0,0)
library(lme4); library(lmerTest)
summary(lm(FruitSizeStd~mamMdiff, data = tdwg_final))
plot(propMegaPalm~mamMdiff, data = subset(tdwg_final_subset2, megapalm_nsp >0))
plot(log(meanFruitLengthFilled)~mamSdiff, data = subset(tdwg_final_subset2))
boxplot(FruitSizeStd~THREEREALM, tdwg_final)

meanFruitByRealm <- tapply(tdwg_final$meanFruitLengthFilled, INDEX = tdwg_final$THREEREALM, FUN = mean)
meanFruitByRealm <- data.frame(meanFruitByRealm, THREEREALM = names(meanFruitByRealm))
tdwg_final <- merge(tdwg_final, meanFruitByRealm)
tdwg_final$FruitSizeStd <- (tdwg_final$meanFruitLengthFilled - tdwg_final$meanFruitByRealm) / sd(tdwg_final$meanFruitByRealm)


fruitSize_mammalBodySizePlot_noisland  <- ggplot(aes(y = log(meanFruitLengthFilled), x = log(value)), data = subset(test, Geology_3categ %in% c("continental", "mainland"))) +
  geom_point() +
  facet_wrap(variable~THREEREALM) +
  geom_smooth(method = "lm") 

ggsave(fruitSize_mammalBodySizePlot_noisland, file = file.path(fig.dir, "fruitVsMammalPlot_noisland.pdf"))



summary(lm(log(meanFruitLengthFilled)~ mamMdiff, data = subset(tdwg_final, Geology_3categ %in% c("continental", "mainland"))))

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

