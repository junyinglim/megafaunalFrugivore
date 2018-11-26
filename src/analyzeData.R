# Analyze Data

## Directories ========================
main.dir <- "/Users/junyinglim/Dropbox/Projects/2019/palms/projects/megafaunalFrugivore"
data.dir <- file.path(main.dir, "data")
res.dir <- file.path(main.dir, "results")
src.dir <- file.path(main.dir, "src")
fig.dir <- file.path(main.dir, "figs")
rawdata.dir <- "~/Dropbox/Projects/2019/palms/data"
frug.dir <- file.path(rawdata.dir, "frugivores")

## Packages ========================
library(plyr)
library(ggplot2); library(sp); library(rgdal); library(viridis)
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

# Import mammal occ datasets
mammal_presnat_occ <- read.table(file.path(data.dir, "mammal_presNat_nobuff_occ.txt"), header = TRUE)
mammal_curr_occ <- read.table(file.path(data.dir, "mammal_current_occ.txt"), header = TRUE)

## Calculate tdwg level metrics ======================== 
palm_occ_trait <- merge(palm_occ, palm_trait[c("SpecName", "AverageFruitLength_cm")])
palm_occ_trait <- palm_occ_trait[!palm_occ_trait$SpecName %in% c("Lodoicea_maldivica", "Cocos_nucifera"),]
tdwg_meanFruit <- ddply(.data = palm_occ_trait,
                        .variables = .(Area_code_L3),
                        .fun = summarise,
                        meanFruitLength = mean(AverageFruitLength_cm, na.rm = TRUE),
                        palm_nSp = length(AverageFruitLength_cm))
names(tdwg_meanFruit)[names(tdwg_meanFruit) == "Area_code_L3"] <- "LEVEL_3_CO"

mammal_presnat_occ_trait <- merge(mammal_presnat_occ, phylacine_trait, by.x = "SpecName", by.y = "Binomial.1.2", all.x = TRUE)
tdwg_presnat_meanBodySize <- ddply(.data = mammal_presnat_occ_trait,
                                   .variables = .(LEVEL_3_CO),
                                   .fun = summarize,
                                   presNat_meanBodySize = mean(Mass.g, na.rm = TRUE),
                                   presNat_nSp = length(Mass.g))

mammal_curr_occ_trait <- merge(mammal_curr_occ, phylacine_trait, by.x = "SpecName", by.y = "Binomial.1.2", all.x = TRUE)
tdwg_curr_meanBodySize <- ddply(.data = mammal_curr_occ_trait,
                                .variables = .(LEVEL_3_CO),
                                .fun = summarize,
                                curr_meanBodySize = mean(Mass.g, na.rm = TRUE),
                                curr_nSp = length(Mass.g))
tdwg_meanBodySize <- merge(tdwg_curr_meanBodySize, tdwg_presnat_meanBodySize, all = TRUE)
tdwg_res <- merge(tdwg_meanFruit, tdwg_meanBodySize, all.x = TRUE)

tdwg_res[which.max(tdwg_res$meanFruitLength),]
subset(palm_occ, Area_code_L3 %in% "MRQ")
subset(palm_trait, SpecName %in% c("Pelagodoxa_henryana", "Cocos_nucifera"))
subset(palm_trait, SpecName %in% c("Washingtonia_robusta"))
subset(palm_occ, Area_code_L3 %in% "MXI")

# Merge
write.csv(tdwg_res, file.path(res.dir, "tdwg_res.csv"), row.names = FALSE)

## Plot mammal body size and fruit sizes ========================
# Import tdwg shape files
tdwg_shape <- readOGR(file.path(rawdata.dir, "TDWG", "level3", "level3.shp"))
tdwg_shape@data$id <- rownames(tdwg_shape@data)
#tdwg_shape_reproj <- spTransform(tdwg_shape, CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
tdwg_shape_df <- fortify(tdwg_shape, region = "id")
tdwg_shape_df <- merge(tdwg_shape_df, tdwg_shape@data)

tdwg_plot_df <- merge(tdwg_res, tdwg_env, by = "LEVEL_3_CO")

fruitLengthPlot  <- ggplot() +
  geom_polygon(aes(y = lat, x = long, group = group), data = subset(tdwg_shape_df, !LEVEL_NAME == "Antarctica")) +
  geom_point(aes(x = LONG, y = LAT, fill = meanFruitLength, size = meanFruitLength), colour = "white", shape = 21, data = tdwg_plot_df, alpha = 0.7) +
  scale_fill_viridis()
ggsave(filename = file.path(fig.dir, "fruitSizePlot.pdf"), fruitLengthPlot)

target.col <- c("LAT", "LONG", "LEVEL_3_CO", "curr_meanBodySize", "curr_nSp", "presNat_meanBodySize","presNat_nSp")
tdwg_plot_melt <- melt(tdwg_plot_df[,target.col], id.vars = c("LAT", "LONG", "LEVEL_3_CO"), measure.vars = c("curr_meanBodySize", "curr_nSp", "presNat_meanBodySize","presNat_nSp"))

meanBodySizePlot <- ggplot() +
  geom_polygon(aes(y = lat, x = long, group = group), data = subset(tdwg_shape_df, !LEVEL_3_CO == "ANT")) +
  geom_point(aes(x = LONG, y = LAT, fill = value, size = value), colour = "white", shape = 21, data = subset(tdwg_plot_melt, variable %in% c("curr_meanBodySize", "presNat_meanBodySize")), alpha = 0.7) +
  facet_wrap(~variable, nrow = 2) + 
  scale_fill_viridis()

ggsave(meanBodySizePlot, filename = file.path(fig.dir, "meanBodySizeCombined.pdf"), width = 10, height = 10)


presNat_meanBodySizePlot  <- ggplot() +
  geom_polygon(aes(y = lat, x = long, group = group), data = subset(tdwg_shape_df, !LEVEL_3_CO == "ANT")) +
  geom_point(aes(x = LONG, y = LAT, fill = presNat_meanBodySize, size = presNat_meanBodySize), colour = "white", shape = 21, data = tdwg_plot_df, alpha = 0.7) +
  scale_fill_viridis()
ggsave(filename = file.path(fig.dir, "presNat_meanBodySizePlot.pdf"), presNat_meanBodySizePlot,width = 10, height = 4)

curr_meanBodySizePlot  <- ggplot() +
  geom_polygon(aes(y = lat, x = long, group = group), data = subset(tdwg_shape_df, !LEVEL_NAME == "Antarctica")) +
  geom_point(aes(x = LONG, y = LAT, fill = curr_meanBodySize, size = curr_meanBodySize), colour = "white", shape = 21, data = tdwg_plot_df, alpha = 0.7) +
  scale_fill_viridis()
ggsave(filename = file.path(fig.dir, "curr_meanBodySizePlot.pdf"), curr_meanBodySizePlot)



library(raster)
x <- raster(file.path(phylacine.dir, "Present_Natural", "Loxodonta_africana.tif"))
y <- raster(file.path(phylacine.dir, "Current", "Loxodonta_africana.tif"))
x_poly <-rasterToPolygons(x, fun = function(x){x>0})
x_fort <- fortify(x_poly)
mean(subset(phylacine_trait, Binomial.1.2 %in% subset(mammal_presnat_occ, LEVEL_3_CO == "CNY")$SpecName)$Mass.g)

ggplot() +
  geom_polygon(aes(y = lat, x = long, group = group), data = subset(tdwg_shape_df, !LEVEL_NAME == "Antarctica")) +
  geom_polygon(aes(y = lat, x = long, group = group), data = x_fort)
head(x_fort)
proj4string(x_poly)

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


