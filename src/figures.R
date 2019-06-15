## GENERATE FIGURES FOR MANUSCRIPT AND DATA EXPLORATION
# Author: Jun Ying Lim

## Packages ========================
rm(list = ls())
library(ggplot2)
library(reshape2)
library(cowplot)
library(wesanderson); library(RColorBrewer); library(viridis)
library(scales)
library(stringr)
library(png); library(grid)

## Directories ========================
main.dir <- "/Users/junyinglim/Dropbox/Projects/2019/palms/projects/megafaunalFrugivore"
data.dir <- file.path(main.dir, "data")
res.dir <- file.path(main.dir, "results")
src.dir <- file.path(main.dir, "src")
fig.dir <- file.path(main.dir, "figs")
rawdata.dir <- "~/Dropbox/Projects/2019/palms/data"
frug.dir <- file.path(rawdata.dir, "frugivores")
options(stringsAsFactors =FALSE)

## Utility functions ========================
cleanCuts <- function(x, round = TRUE){
  # Cleans up vector of factor levels generated using `cut` function for nicer plotting
  x <- gsub(",", " - ",x)
  x <- gsub("\\[|\\(|\\]|\\)", "", x)
  if(round){
    num <- str_split_fixed(x, pattern = " - ", n=2)
    x <- paste(round(as.numeric(num[,1]), 3), round(as.numeric(num[,2]), 3), sep = " \u2013 ")
  }
  return(x)
}

roundCorrectly <- function(x){
  if(! is.na(x)){
    return(as.character(as.expression(substitute(italic(R)^2~"="~x, list(x = formatC(x, format = "f", digits = 2))))))  
  } else {
    return(NA)
  }
}

## Import tdwg polygon data ========================
# tdwg_shp_raw <- readOGR("~/Dropbox/Projects/2019/palms/data/TDWG/level3/level3.shp")
# tdwg_shp@data$id <- rownames(tdwg_shp@data)
# tdwg_fort <- fortify(tdwg_shp, region = "id")
# tdwg_fort <- merge(tdwg_fort, tdwg_shp@data, by = "id")
# saveRDS(tdwg_fort, file.path(data.dir, "tdwg_shp_fort.rds") )
tdwg_shp <- readRDS(file.path(data.dir, "tdwg_shp_fort.rds"))
tdwg_shp2 <- subset(tdwg_shp, !LEVEL_3_CO == "ANT")

## Import tdwg environmental data ========================
tdwg_env <- read.csv(file.path(data.dir, "TDWG_Environment_AllData_2014Dec.csv"))

## Import tdwg mammal and palm data ========================
tdwg_res <- read.csv(file.path(res.dir, "tdwg_mammal.csv"))

## Merge data ========================
tdwg_final <- read.csv(file.path(data.dir,"tdwg_final.csv"))

## Plotting utilities ========================
# Define map theme
map_theme <- theme(panel.background = element_blank(),
                   axis.text = element_blank(),
                   axis.line = element_blank(),
                   axis.ticks = element_blank(),
                   axis.title = element_blank(),
                   legend.position = "bottom",
                   legend.justification = "center",
                   strip.background = element_blank(),
                   strip.text = element_text(size = 18))

point_theme <- theme(
                     axis.text = element_text(color = "grey20"),
                     legend.position = "bottom",
                     strip.background = element_blank(),
                     strip.text = element_text(size = 18))

sloth_grob <- readPNG(file.path(fig.dir, "symbols", "megasloth.png"))
rhino_grob <- readPNG(file.path(fig.dir, "symbols", "rhino.png"))

## Plot biogeography panels =========================
neot_countrylist <- subset(tdwg_final, THREEREALM == "NewWorld")$LEVEL_3_CO
afrot_countrylist <- subset(tdwg_final, THREEREALM == "OWWest")$LEVEL_3_CO
indot_countrylist <- subset(tdwg_final, THREEREALM == "OWEast")$LEVEL_3_CO

tdwg_shp2$Global <- ifelse(tdwg_shp2$LEVEL_3_CO %in% tdwg_final$LEVEL_3_CO, 1, 0)
tdwg_shp2$Neotropics <- ifelse(tdwg_shp2$LEVEL_3_CO %in% neot_countrylist, 1, 0)
tdwg_shp2$Afrotropics <- ifelse(tdwg_shp2$LEVEL_3_CO %in% afrot_countrylist, 1, 0)
tdwg_shp2$Indotropics <- ifelse(tdwg_shp2$LEVEL_3_CO %in% indot_countrylist, 1, 0)

map_globaltropics <- ggplot() +
  geom_polygon(aes(y = lat, x = long, group = group, fill = factor(Global)),
               data = tdwg_shp2) + 
  scale_fill_manual(values = c("grey60", "grey20")) +
  map_theme + theme(legend.position = "none")
ggsave(map_globaltropics,
       filename = file.path(fig.dir, "map_globaltropics.pdf"), height = 6, width = 14)

map_neotropics <- ggplot() +
  geom_polygon(aes(y = lat, x = long, group = group, fill = factor(Neotropics)),
               data = tdwg_shp2) + 
  scale_fill_manual(values = c("grey60", wes_palette("Cavalcanti1", n = 5)[1])) +
  map_theme + theme(legend.position = "none")
ggsave(map_neotropics, filename = file.path(fig.dir, "map_neotropics.pdf"), height = 6, width = 14)

map_afrotropics <- ggplot() +
  geom_polygon(aes(y = lat, x = long, group = group, fill = factor(Afrotropics)),
                                    data = tdwg_shp2) + 
  scale_fill_manual(values = c("grey60", wes_palette("Cavalcanti1", n = 5)[2])) +
  map_theme + theme(legend.position = "none")
ggsave(map_afrotropics, filename = file.path(fig.dir, "map_afrotropics.pdf"), height = 6, width = 14)

map_indotropics <- ggplot() +
  geom_polygon(aes(y = lat, x = long, group = group, fill = factor(Indotropics)), data = tdwg_shp2) +   scale_fill_manual(values = c("grey60", wes_palette("Cavalcanti1", n = 5)[5])) +
  map_theme + theme(legend.position = "none")
ggsave(map_indotropics, filename = file.path(fig.dir, "map_indotropics.pdf"), height = 6, width = 14)

## Fig 1: Fruit size and body size ========================
# Define some custom themes
title_theme <- theme(plot.title = element_blank())
  #theme(plot.title = element_text(size = 20, hjust = 0))

# Plot fruit sizes
target.col <- c("LAT", "LONG", "LEVEL_3_CO", "medianFruitLengthFilled")
tdwg_medFS <- melt(tdwg_final[,target.col], id.vars = c("LAT", "LONG", "LEVEL_3_CO"), measure.vars = "medianFruitLengthFilled")

q_probs = seq(0.0, 1.0, 0.1)
tdwg_medFS$medFS_q <- cut(tdwg_medFS$value, quantile(tdwg_medFS$value, probs = q_probs, na.rm = T), include.lowest = T)
levels(tdwg_medFS$medFS_q) <- gsub(levels(tdwg_medFS$medFS_q), pattern = ",", replacement = " - ")
levels(tdwg_medFS$medFS_q) <- gsub(levels(tdwg_medFS$medFS_q), pattern = "\\(|\\[|\\]", replacement = "")

medFS_p <- ggplot() +
  geom_polygon(aes(y = lat, x = long, group = group), fill = "grey60",
              data = tdwg_shp2) +
  geom_point(aes(x = LONG, y = LAT,
                 colour = medFS_q, size = value), data = tdwg_medFS, alpha = 0.9) +
  coord_fixed() +
  ggtitle("Median fruit size") +
  map_theme + title_theme +
  guides(size=FALSE, colour = guide_legend(title = "Fruit length (cm)", override.aes = list(size=5))) +
  scale_color_viridis(discrete = TRUE)
ggsave(file.path(fig.dir, "medFS.pdf"), medFS_p, width = 9, height = 5)

target.col <- c("LAT", "LONG", "LEVEL_3_CO", "max95FruitLengthFilled")
tdwg_maxFS <- melt(tdwg_final[,target.col], id.vars = c("LAT", "LONG", "LEVEL_3_CO"), measure.vars = "max95FruitLengthFilled")
q_probs = seq(0.0, 1.0, 0.1)
tdwg_maxFS$maxFS_q <- cut(tdwg_maxFS$value, quantile(tdwg_maxFS$value, probs = q_probs, na.rm = T), include.lowest = T)
levels(tdwg_maxFS$maxFS_q) <- gsub(levels(tdwg_maxFS$maxFS_q), pattern = ",", replacement = " - ")
levels(tdwg_maxFS$maxFS_q) <- gsub(levels(tdwg_maxFS$maxFS_q), pattern = "\\(|\\[|\\]", replacement = "")

maxFS_p <- ggplot() +
  geom_polygon(aes(y = lat, x = long, group = group), fill = "grey60",
               data = tdwg_shp2) +
  geom_point(aes(x = LONG, y = LAT,
                 colour = maxFS_q, size = value),
             data =tdwg_maxFS, alpha = 0.9) +
  coord_fixed() +
  ggtitle("Maximum (95th percentile) fruit size") +
  map_theme + title_theme +
  guides(size=FALSE, colour = guide_legend(title = "Fruit length (cm)", override.aes = list(size=5))) +
  scale_color_viridis(discrete = TRUE)

ggsave(file.path(fig.dir, "maxFS.pdf"), maxFS_p, width = 9, height = 5)

# Plotting mammal data on map
target.col <- c("LAT", "LONG", "LEVEL_3_CO", "curr_medianBodySize", "presNat_medianBodySize","curr_max95BodySize","presNat_max95BodySize")

tdwg_medBS_melt <- melt(tdwg_final[,target.col],
                        id.vars = c("LAT", "LONG", "LEVEL_3_CO"),
                        measure.vars = c("curr_medianBodySize",
                                         "presNat_medianBodySize"))
tdwg_medBS_melt$variable <- factor(tdwg_medBS_melt$variable, levels = c("curr_medianBodySize", "presNat_medianBodySize"), labels = c("Current", "Present-natural"))

medBS_arthm_min <- min(tdwg_medBS_melt$value)/1000
medBS_arthm_max <- max(tdwg_medBS_melt$value)/1000

medBS_p <- ggplot() + 
  geom_polygon(aes(y = lat, x = long, group = group), fill = "grey60", data = tdwg_shp2) +
  geom_point(aes(y = LAT, x= LONG, colour = (value/1000), size = log(value)), data = tdwg_medBS_melt) +
  facet_wrap(~variable, nrow = 2)+
  # some padding so the legend takes up same space as the maximum BS
  scale_colour_viridis(discrete = FALSE,
                       breaks = c(medBS_arthm_min, 0.1, 1, 10, medBS_arthm_max),
                       trans = log10_trans(),
                       labels = c(round(medBS_arthm_min,2), 0.1, 1, 10, round(medBS_arthm_max) )) +
  coord_fixed() +
  ggtitle("Median body size") +
  guides(size = FALSE, colour = guide_colorbar(title = "Body mass (kg)",
                                               label.theme = element_text(angle = 45, vjust = 0.5))) +
  map_theme +
  title_theme +
  theme(legend.key.width = unit(1, "cm"))
ggsave(file.path(fig.dir, "medBS.pdf"), medBS_p, width = 9, height = 10)

tdwg_maxBS_melt <- melt(tdwg_final[,target.col],
                        id.vars = c("LAT", "LONG", "LEVEL_3_CO"),
                        measure.vars = c("curr_max95BodySize",
                                         "presNat_max95BodySize"))
tdwg_maxBS_melt$variable <- factor(tdwg_maxBS_melt$variable, levels = c("curr_max95BodySize", "presNat_max95BodySize"), labels = c("Current", "Present-natural"))

maxBS_arthm_min <- min(tdwg_maxBS_melt$value)/1000
maxBS_arthm_max <- max(tdwg_maxBS_melt$value)/1000

maxBS_p <- ggplot() + 
  geom_polygon(aes(y = lat, x = long, group = group), fill = "grey60", data = tdwg_shp2) +
  geom_point(aes(y = LAT, x= LONG, colour = value/1000, size = log(value)), data = tdwg_maxBS_melt) +
  facet_wrap(~variable, nrow = 2)+
  coord_fixed() +
  scale_colour_viridis(discrete = FALSE,
                       breaks = c(maxBS_arthm_min, 0.1, 1, 10, 100, 1000, maxBS_arthm_max),
                       trans = log10_trans(),
                       labels = c(round(maxBS_arthm_min, 2), 0.1, 1, 10, 100, 1000, round(maxBS_arthm_max))) +  
  ggtitle("Maximum (95th percentile) body size") +
  guides(size = FALSE,
         colour = guide_colorbar(title = "Body mass (kg)",
                                 label.theme = element_text(angle = 45, vjust = 0.5))) +
  map_theme +
  title_theme +
  theme(legend.key.width = unit(1, "cm"))

ggsave(file.path(fig.dir, "maxBS.pdf"), maxBS_p, width = 9, height = 10)

# Combine fruit size and body size 
fig1_label_size = 20
FS_comb <- plot_grid(medFS_p + theme(legend.position = "none"),
                     maxFS_p + theme(legend.position = "none"),
                     nrow = 1, labels = c("a","b"), label_size = fig1_label_size)

FS_comb_leg <- plot_grid(get_legend(medFS_p), get_legend(maxFS_p), nrow = 1)
FS_comb_wleg <- plot_grid(FS_comb, FS_comb_leg, rel_heights = c(1, 0.2), nrow = 2)

BS_comb <- plot_grid(medBS_p + theme(legend.position = "none"),
                     maxBS_p + theme(legend.position = "none"), 
                     nrow = 1, labels= c("c", "d"), label_size = fig1_label_size)

BS_comb_leg <- plot_grid(get_legend(medBS_p), get_legend(maxBS_p), nrow = 1)

fig1_BSFS <- plot_grid(FS_comb, FS_comb_leg, BS_comb, BS_comb_leg,
                       nrow = 4, rel_heights = c(1,0.2,2,0.2))

ggsave(file.path(fig.dir, "fig1_BSFS.pdf"), fig1_BSFS, width = 18, height = 12)


# Plot species richness ===========
misc_melt <- melt(tdwg_final,
                    id.vars = c("LAT", "LONG", "LEVEL_3_CO", "THREEREALM"),
                    measure.vars = c("curr_nSp", "presNat_nSp", "palm_nSp",
                                     "curr_mega_nSp", "presNat_mega_nSp", "megapalm_nsp",
                                     "deltaMedianBodySize", "futr_maxBodySize", "futr_medianBodySize"))
misc_maps <- list()
for(i in 1:length(unique(misc_melt$variable))){
  if(i %in% 1:7){
    misc_maps[[i]] <- ggplot() + geom_polygon(aes(y = lat, x = long, group = group), data = tdwg_shp2) + geom_point(aes(y = LAT, x = LONG, fill = scale(value), size = scale(value)), shape = 21, colour = "white", data = subset(misc_melt, variable == unique(misc_melt$variable)[i] )) + map_theme + scale_fill_viridis() + labs(title = unique(misc_melt$variable)[i])+ guides(size = FALSE)
  } else {
    misc_maps[[i]] <- ggplot() + geom_polygon(aes(y = lat, x = long, group = group), data = tdwg_shp2) + geom_point(aes(y = LAT, x = LONG, fill = scale(log(value)), size = scale(log(value))), shape = 21, colour = "white", data = subset(misc_melt, variable == unique(misc_melt$variable)[i] )) + map_theme + scale_fill_viridis() + labs(title = unique(misc_melt$variable)[i]) + guides(size = FALSE)
  }
}
misc_map_comb <- plot_grid(plotlist = misc_maps, nrow = 3, ncol = 3)

ggsave(file.path(fig.dir, "miscMaps.pdf"), misc_map_comb, width = 26, height = 16.5)


mammalSpRichPlot <- ggplot() +
  geom_polygon(aes(y = lat, x = long, group = group),
               data = subset(tdwg_shape_df, !LEVEL_3_CO == "ANT")) +
  geom_point(aes(x = LONG, y = LAT, fill = value, size = value), colour = "white", shape = 21,
             data = subset(tdwg_plot_melt, variable %in% c("curr_nSp", "presNat_nSp")), alpha = 0.7) +
  facet_wrap(~variable, nrow = 2) + 
  map_theme +
  theme(panel.background = element_blank()) +
  scale_fill_viridis()
  
ggsave(mammalSpRichPlot, filename = file.path(fig.dir, "mammalSpRichCombined.pdf"), width = 10, height = 10)

# Plot mammal sp loss on map ==========
megaSpLossPlot <- ggplot() + 
  geom_polygon(aes(y = lat, x = long, group = group),
               data = subset(tdwg_shp, !LEVEL_3_CO == "ANT")) +
  geom_point(aes(y = LAT, x = LONG, fill = megaSpLoss, size = megaSpLoss), data = tdwg_final,
             colour = "white", shape = 21) +
  mammal_theme +
  scale_fill_viridis()
ggsave(file.path(fig.dir, "megaSpLoss.pdf"), megaSpLossPlot, width = 10, height = 5)

megaSpLossCategPlot <- ggplot() + 
  geom_polygon(aes(y = lat, x = long, group = group),
               data = subset(tdwg_shp, !LEVEL_3_CO == "ANT")) +
  geom_point(aes(y = LAT, x = LONG, fill = ifelse(megaSpLoss > median(megaSpLoss), "high", "low")), data = tdwg_final,
             colour = "white", shape = 21, size = 3) +
  mammal_theme 
ggsave(file.path(fig.dir, "megaSpLossCateg.pdf"), megaSpLossCategPlot, width = 10, height = 5)

mesoSpLossPlot <- ggplot() + 
  geom_polygon(aes(y = lat, x = long, group = group),
               data = subset(tdwg_shp, !LEVEL_3_CO == "ANT")) +
  geom_point(aes(y = LAT, x = LONG, fill = mesoSpLoss, size = mesoSpLoss), data = tdwg_final,
             colour = "white", shape = 21) +
  mammal_theme +
  scale_fill_viridis() 
ggsave(file.path(fig.dir, "mesoSpLoss.pdf"), mesoSpLossPlot, width = 10, height = 5)

mesoSpLossCategPlot <- ggplot() + 
  geom_polygon(aes(y = lat, x = long, group = group),
               data = subset(tdwg_shp, !LEVEL_3_CO == "ANT")) +
  geom_point(aes(y = LAT, x = LONG, fill = ifelse(mesoSpLoss > median(mesoSpLoss), "high", "low")), data = tdwg_final,
             colour = "white", shape = 21, size = 3) +
  mammal_theme 
ggsave(file.path(fig.dir, "mesoSpLossCateg.pdf"), mesoSpLossCategPlot, width = 10, height = 5)

par(mfrow = c(1,2))
hist(tdwg_final$megaSpLoss, main = "Absolute loss of diversity")
abline(v = median(tdwg_final$megaSpLoss), col = "red")
z <- tdwg_final$megaSpLoss / tdwg_final$presNat_mega_nSp
z[tdwg_final$megaSpLoss == 0] <- 0
hist(z, main = "Proportional loss of diversity")

# Plot palm species richness ===========
tdwg_final <- merge(tdwg_res, tdwg_env, by = "LEVEL_3_CO",all.x = TRUE)
remoteIslList <- subset(tdwg_final, ISISLAND == 1 & DistToCont_km > 2000) # List of islands to remove
tdwg_final2 <- subset(tdwg_final, !LEVEL_3_CO %in% union(remoteIslList$LEVEL_3_CO, c("MRQ", "TUA", "VNA")))

palmSpRichPlot <- ggplot() +
  geom_polygon(aes(y = lat, x = long, group = group),
               data = subset(tdwg_shape_df, !LEVEL_3_CO == "ANT")) +
  geom_point(aes(x = LONG, y = LAT, fill = PALMSR, size = PALMSR), colour = "white", shape = 21,
             data = tdwg_final2, alpha = 0.7) +
  map_theme +
  guides(size = FALSE) +
  theme(panel.background = element_blank(), legend.position = "bottom") +
  scale_fill_viridis()
ggsave(file.path(fig.dir, "palmSpRichPlot.pdf"), palmSpRichPlot, width = 10, height = 5, device = cairo_pdf)

palmDivPlot <- ggplot() +
  geom_polygon(aes(y = lat, x = long, group = group),
               data = subset(tdwg_shape_df, !LEVEL_3_CO == "ANT")) +
  geom_point(aes(x = LONG, y = LAT, fill = PALMSR>5, size = PALMSR), colour = "white", shape = 21,
             data = tdwg_final2, alpha = 0.7) +
  map_theme +
  guides(size = FALSE) +
  theme(panel.background = element_blank(), legend.position = "bottom")

ggsave(file.path(fig.dir, "palmDiv.pdf"), palmDivPlot, width = 10, height = 5, device = cairo_pdf)



## Plotting all abiotic variables on map  ===============

biotic_melt <- melt(tdwg_final,
                    id.vars = c("LAT", "LONG", "LEVEL_3_CO", "THREEREALM"),
                    measure.vars = c("curr_medianBodySize", "presNat_medianBodySize", "medianFruitLengthFilled",
                                     "curr_max95BodySize", "presNat_max95BodySize", "max95FruitLengthFilled",
                                     "curr_dispBodySize", "presNat_dispBodySize", "dispFruitLengthFilled"))

biotic_maps <- list()
for(i in 1:length(unique(biotic_melt$variable))){
  if(i %in% 7:9){
    biotic_maps[[i]] <- ggplot() + geom_polygon(aes(y = lat, x = long, group = group), data = tdwg_shp2) + geom_point(aes(y = LAT, x = LONG, fill = scale(value), size = scale(value)), shape = 21, colour = "white", data = subset(biotic_melt, variable == unique(biotic_melt$variable)[i] )) + map_theme + scale_fill_viridis() + labs(title = unique(biotic_melt$variable)[i])+ guides(size = FALSE)
  } else {
    biotic_maps[[i]] <- ggplot() + geom_polygon(aes(y = lat, x = long, group = group), data = tdwg_shp2) + geom_point(aes(y = LAT, x = LONG, fill = scale(log(value)), size = scale(log(value))), shape = 21, colour = "white", data = subset(biotic_melt, variable == unique(biotic_melt$variable)[i] )) + map_theme + scale_fill_viridis() + labs(title = unique(biotic_melt$variable)[i]) + guides(size = FALSE)
  }
}
biotic_map_comb <- plot_grid(plotlist = biotic_maps, nrow = 3, ncol = 3)

ggsave(file.path(fig.dir, "bioticMaps.pdf"), biotic_map_comb, width = 26, height = 16.5)

abiotic_melt <- melt(tdwg_final,
                     id.vars = c("LAT", "LONG", "LEVEL_3_CO", "THREEREALM"),
                     measure.vars = c("globalPC1", "globalPC2", "globalPC3",
                                      "regionalPC1", "regionalPC2", "regionalPC3",
                                      "soilcount", "lgm_ens_Tano", "lgm_ens_Pano"))
abiotic_maps <- list()
for(i in 1:length(unique(abiotic_melt$variable))){
  abiotic_maps[[i]] <- ggplot() + geom_polygon(aes(y = lat, x = long, group = group), data = tdwg_shp2) + geom_point(aes(y = LAT, x = LONG, fill = scale(value), size = scale(value)), shape = 21, colour = "white", data = subset(abiotic_melt, variable == unique(abiotic_melt$variable)[i]))  + map_theme + scale_fill_viridis() + guides(size = FALSE) + labs(title = unique(abiotic_melt$variable)[i])
}

abiotic_map_comb <- plot_grid(plotlist = abiotic_maps, nrow = 3, ncol = 3)

ggsave(file.path(fig.dir, "abioticMaps.pdf"), abiotic_map_comb, width = 26, height = 16.5)

## Plot model average relative importance of maximum body size ===============
maxBS_modavg_res <- read.csv(file.path(res.dir, "maxBS_modavg_res.csv"), stringsAsFactors = TRUE)
maxBS_modavg_res$Variable <- factor(maxBS_modavg_res$Variable,
                                       levels = c("curr_logMax95BS_scl", "pnat_logMax95BS_scl",
                                                  "globalPC1_scl", "globalPC2_scl", "globalPC3_scl",
                                                  "regionalPC1_scl", "regionalPC2_scl", "regionalPC3_scl",
                                                  "lgm_ens_Pano_scl", "lgm_ens_Tano_scl"),
                                       labels = c("Max. body size", "Max. body size",
                                                  "PC1", "PC2", "PC3",
                                                  "PC1", "PC2", "PC3",
                                                  "LGM Prec. anomaly", "LGM Temp. anomaly"))

maxBS_modavg_res$GeographicScale <- factor(maxBS_modavg_res$GeographicScale, levels = c("Global", "Afrotropics", "Neotropics", "Indotropics"))

maxBS_modavg_plot <- ggplot(data = maxBS_modavg_res) + 
  geom_hline(aes(yintercept = 0), linetype = "dashed", colour = "grey50") +
  geom_point(aes(y = avgStdCoef, x = Variable, colour = Scenario, size = var_impt), position = position_dodge(0.8)) + 
  geom_errorbar(aes(ymin = avgStdCIlower, ymax = avgStdCIupper, x = Variable, colour = Scenario), position = position_dodge(0.8), width = 0.1 ) +
  scale_size_continuous(limits = c(0, 1), name = "Variable\nimportance") +
  labs(title = "Body size 95th percentile\n(scaled)", y = "Standardized\nmodel averaged coefficients", x= "Variables") +
  facet_wrap(~ GeographicScale, drop = TRUE, nrow = 4) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
ggsave(file.path(fig.dir, "maxBS_modavg.pdf"), maxBS_modavg_plot, height= 9, width = 8)

## Fig 3: Plot variance explained in full BS models ===============
varexp_theme <- theme(axis.text.x = element_text(angle = 46, hjust = 1),
                      strip.background = element_blank(),
                      strip.text = element_text(color = "grey20", size = 18))
r2_labelsize <- 5
medBS_col <- wes_palette("GrandBudapest1", n = 4)[2]

medBS_levels <- c("curr_logMedBS_scl", "pnat_logMedBS_scl", "globalPC1_scl", "globalPC2_scl", "globalPC3_scl", "regionalPC1_scl", "regionalPC2_scl", "regionalPC3_scl","lgm_ens_Pano_scl", "lgm_ens_Tano_scl")
medBS_labels <- c("Median body size", "Median body size","Climate PC1", "Climate PC2", "Climate PC3",
                  "Climate PC1", "Climate PC2", "Climate PC3","LGM Prec. Anom.", "LGM Temp. Anom.")

medBS_ols_modavg_res <- read.csv(file.path(res.dir, "medBS_ols_modavg.csv"), stringsAsFactors = TRUE)
medBS_ols_modavg_res$Variable <- factor(medBS_ols_modavg_res$coefficient,
                                    levels = medBS_levels, labels = medBS_labels)

medBS_ols_modavg_res$Geographic.Scale <- factor(medBS_ols_modavg_res$Geographic.Scale, levels = c("Global", "Afrotropics", "Neotropics", "Indotropics"))
medBS_ols_modavg_res$totalR2[duplicated(medBS_ols_modavg_res$totalR2)] <- NA
medBS_ols_modavg_res$totalR2plot <- sapply(X = medBS_ols_modavg_res$totalR2, FUN = roundCorrectly)
medBS_ols_modavg_plot <-ggplot(data = subset(medBS_ols_modavg_res, !is.na(Variable) )) + 
  geom_bar(aes(y = varexp,  x = Variable, fill = Variable), stat = "identity") +
  facet_grid(Geographic.Scale ~ Scenario) +
  labs(y = "Variance explained", x = "") +
  geom_text(aes(y = Inf , x = Inf, label = ifelse(is.na(totalR2),"",totalR2plot)), hjust = 1, vjust = 1, parse = T, size = r2_labelsize) +
  scale_y_continuous(breaks = c(0,0.1, 0.2, 0.3)) +
  guides(fill = FALSE) +
  scale_fill_manual(values = c(medBS_col, rep("grey20",5))) +
  varexp_theme
ggsave(file.path(fig.dir, "medBS_ols_modavg_relaimpo.pdf"),
       medBS_ols_modavg_plot, height= 12, width = 8)

medBS_sar_modavg_res <- read.csv(file.path(res.dir, "medBS_sar_modavg.csv"), stringsAsFactors = TRUE)
medBS_sar_modavg_res$Scenario <- factor(medBS_sar_modavg_res$Scenario, levels = 
                                          c("Current", "Present-Natural"), labels = c("Current", "Present-"))
medBS_sar_modavg_res$Variable <- factor(medBS_sar_modavg_res$coefficient,
                                        levels = medBS_levels, labels = medBS_labels)

medBS_sar_modavg_res$Geographic.Scale <- factor(medBS_sar_modavg_res$Geographic.Scale, levels = c("Global", "Afrotropics", "Neotropics", "Indotropics"))
medBS_sar_modavg_res$totalR2[duplicated(medBS_sar_modavg_res$totalR2)] <- NA
medBS_sar_modavg_res$totalR2plot <- sapply(X = medBS_sar_modavg_res$totalR2, FUN = roundCorrectly)
medBS_sar_modavg_plot <-ggplot(data = subset(medBS_sar_modavg_res, !is.na(Variable) )) + 
  geom_bar(aes(y = varexp,  x = Variable, fill = Variable), stat = "identity") +
  facet_grid(Geographic.Scale ~ Scenario) +
  labs(y = "Variance explained", x = "") +
  geom_text(aes(y = Inf , x = Inf, label = ifelse(is.na(totalR2),"",totalR2plot)), hjust = 1, vjust = 1, parse = T, size = r2_labelsize) +
  guides(fill = FALSE) +
  scale_y_continuous(breaks = c(0,0.1, 0.2, 0.3)) +
  scale_fill_manual(values = c(medBS_col, rep("grey20",5))) +
  varexp_theme
ggsave(file.path(fig.dir, "medBS_sar_modavg_relaimpo.pdf"),
       medBS_sar_modavg_plot, height= 12, width = 8)

maxBS_ols_modavg_res <- read.csv(file.path(res.dir, "maxBS_ols_modavg.csv"), stringsAsFactors = TRUE)
maxBS_levels <- c("curr_logMax95BS_scl", "pnat_logMax95BS_scl", "globalPC1_scl", "globalPC2_scl", "globalPC3_scl", "regionalPC1_scl", "regionalPC2_scl", "regionalPC3_scl","lgm_ens_Pano_scl", "lgm_ens_Tano_scl")
maxBS_labels <- c("Maximum body size", "Maximum body size","Climate PC1", "Climate PC2", "Climate PC3",
                  "Climate PC1", "Climate PC2", "Climate PC3","LGM Prec. Anom.", "LGM Temp. Anom.")

maxBS_ols_modavg_res$coefficient <- factor(maxBS_ols_modavg_res$coefficient,
                                           levels = maxBS_levels,
                                           labels = maxBS_labels)

maxBS_ols_modavg_res$Geographic.Scale <- factor(maxBS_ols_modavg_res$Geographic.Scale, levels = c("Global", "Afrotropics", "Neotropics", "Indotropics"))
maxBS_ols_modavg_res$totalR2[duplicated(maxBS_ols_modavg_res$totalR2)] <- NA
maxBS_ols_modavg_res$totalR2plot <- sapply(X = maxBS_ols_modavg_res$totalR2, FUN = roundCorrectly)

maxBS_col <- wes_palette("GrandBudapest1", n = 4)[4]
maxBS_ols_modavg_plot <- ggplot(data = subset(maxBS_ols_modavg_res, !is.na(coefficient))) + 
  geom_bar(aes(y = varexp,  x = coefficient, fill = coefficient), stat = "identity") +
  facet_grid(Geographic.Scale ~ Scenario) +
  labs(y = "Variance explained", x = "") +
  geom_text(aes(y = Inf , x = Inf, label = ifelse(is.na(totalR2),"",totalR2plot)), hjust = 1, vjust = 1, parse = T, size = r2_labelsize) +
  guides(fill = FALSE) +
  scale_y_continuous(breaks = c(0,0.1, 0.2, 0.3)) +
  scale_fill_manual(values = c(maxBS_col, rep("grey20",5))) +
  varexp_theme
ggsave(file.path(fig.dir, "maxBS_ols_modavg_relaimpo.pdf"), maxBS_ols_modavg_plot, height= 12, width = 8)

maxBS_sar_modavg_res <- read.csv(file.path(res.dir, "maxBS_sar_modavg.csv"))
maxBS_sar_modavg_res$coefficient <- factor(maxBS_sar_modavg_res$coefficient,
                                           levels = maxBS_levels,
                                           labels = maxBS_labels)

maxBS_sar_modavg_res$Geographic.Scale <- factor(maxBS_sar_modavg_res$Geographic.Scale, levels = c("Global", "Afrotropics", "Neotropics", "Indotropics"))
maxBS_sar_modavg_res$totalR2[duplicated(maxBS_sar_modavg_res$totalR2)] <- NA
maxBS_sar_modavg_res$totalR2plot <- sapply(X = maxBS_sar_modavg_res$totalR2, FUN = roundCorrectly)

maxBS_sar_modavg_plot <- ggplot(data = subset(maxBS_sar_modavg_res, !is.na(coefficient))) + 
  geom_bar(aes(y = varexp,  x = coefficient, fill = coefficient), stat = "identity") +
  facet_grid(Geographic.Scale ~ Scenario) +
  labs(y = "Variance explained", x = "") +
  geom_text(aes(y = Inf , x = Inf, label = ifelse(is.na(totalR2),"",totalR2plot)), hjust = 1, vjust = 1, parse = T, size = r2_labelsize) +
  guides(fill = FALSE) +
  scale_y_continuous(breaks = c(0,0.1, 0.2, 0.3)) +
  scale_fill_manual(values = c(maxBS_col, rep("grey20",5))) +
  varexp_theme

ggsave(file.path(fig.dir, "maxBS_sar_modavg_relaimpo.pdf"), maxBS_sar_modavg_plot, height= 12, width = 8)



## Fig 2: Plotting partial residuals =========
maxBS_pnat_presid <- readRDS(file.path(res.dir, "maxBS_pnat_presid.rds"))
maxBS_curr_presid <- readRDS(file.path(res.dir, "maxBS_curr_presid.rds"))
medBS_pnat_presid <- readRDS(file.path(res.dir, "medBS_pnat_presid.rds"))
medBS_curr_presid <- readRDS(file.path(res.dir, "medBS_curr_presid.rds"))

library(wesanderson)
curr_col <- wes_palette("Cavalcanti1",n = 4)[1]
pnat_col <- wes_palette("Cavalcanti1",n = 4)[2]
presid_theme <- theme(panel.background = element_blank(), panel.border = element_rect(fill = "transparent"))

medBS_curr_presid_p <- ggplot(data = medBS_curr_presid$points) + geom_point(aes(y = presid, x = curr_logMedBS_scl),shape= 21, color = curr_col) + geom_abline(slope = medBS_curr_presid$slope, intercept = medBS_curr_presid$intercept, color = curr_col) + labs(x = "Current log median body size\n(standard deviations)", y = "Log median fruit size \n(partial residuals)") + presid_theme

medBS_pnat_presid_p <- ggplot(data = medBS_pnat_presid$points) + geom_point(aes(y = presid, x = pnat_logMedBS_scl),shape= 21, color = pnat_col) + geom_abline(slope = medBS_pnat_presid$slope, intercept = medBS_pnat_presid$intercept, color = pnat_col) + labs(x = "Present-natural log median body size\n(standard deviations)", y = "Log median fruit size \n(partial residuals)") + presid_theme

maxBS_curr_presid_p <- ggplot(data = maxBS_curr_presid$points) + geom_point(aes(y = presid, x = curr_logMax95BS_scl),shape= 21, color = curr_col) + geom_abline(slope = maxBS_curr_presid$slope, intercept = maxBS_curr_presid$intercept, color = curr_col) + labs(x = "Current log maximum body size\n(standard deviations)", y = "Log maximum fruit size \n(partial residuals)") + presid_theme

maxBS_pnat_presid_p <- ggplot(data = maxBS_pnat_presid$points) + geom_point(aes(y = presid, x = pnat_logMax95BS_scl),shape= 21, color = pnat_col) + geom_abline(slope = maxBS_pnat_presid$slope, intercept = maxBS_pnat_presid$intercept, color = pnat_col) + labs(x = "Present-natural log maximum body size\n(standard deviations)", y = "Log maximum fruit size \n(partial residuals)") + presid_theme

presid_comb_p <- plot_grid(medBS_curr_presid_p, medBS_pnat_presid_p, maxBS_curr_presid_p, maxBS_pnat_presid_p, nrow = 2, labels = "AUTO", rel_widths = c(1,1))
ggsave(file.path(fig.dir, "presid_comb.pdf"), presid_comb_p, width = 9, height = 8)

## Plotting partial residuals (alternate combined plot)
names(maxBS_curr_presid$points) <- c("presid", "x")
names(maxBS_pnat_presid$points) <- c("presid", "x")
maxBS_comb_presid <- rbind(maxBS_curr_presid$points, maxBS_pnat_presid$points)
maxBS_comb_presid$Scenario <- rep(c("Current", "Present-natural"), each = 129)

comb_theme <- theme(legend.position = c(0.05,0.9),
      legend.background = element_rect(color = "grey90"),
      panel.background = element_blank(),
      legend.box.background = element_rect(fill = "white"),
      legend.title = element_blank(),
      panel.border = element_rect(fill = NA, color = "grey90"))

max_comb_presid_p <- ggplot(data= maxBS_comb_presid) +
  geom_point(aes(y = presid, x = x, color = Scenario), shape = 21)+ labs(x = "Log maximum body size\n(standard deviations)", y = "Log maximum fruit size \n(partial residuals)") + 
  geom_abline(slope = maxBS_curr_presid$slope, intercept = maxBS_curr_presid$intercept, colour = curr_col) +
  geom_abline(slope = maxBS_pnat_presid$slope, intercept = maxBS_pnat_presid$intercept, colour = pnat_col) +
  comb_theme +
  scale_color_manual(values = c(curr_col, pnat_col))

names(medBS_curr_presid$points) <- c("presid", "x")
names(medBS_pnat_presid$points) <- c("presid", "x")
medBS_comb_presid <- rbind(medBS_curr_presid$points, medBS_pnat_presid$points)
medBS_comb_presid$Scenario <- rep(c("Current", "Present-Natural"), each = 129)

med_comb_presid_p <- ggplot(data= medBS_comb_presid) +
  geom_point(aes(y = presid, x = x, color = Scenario), shape = 21) + labs(x = "Log median body size\n(standard deviations)", y = "Log median fruit size \n(partial residuals)") + 
  geom_abline(slope = medBS_curr_presid$slope, intercept = medBS_curr_presid$intercept, colour = curr_col) +
  geom_abline(slope = medBS_pnat_presid$slope, intercept = medBS_pnat_presid$intercept, colour = pnat_col) +
  comb_theme +
  scale_color_manual(values = c(curr_col, pnat_col))

presid_comb_p2 <- plot_grid(med_comb_presid_p, max_comb_presid_p, labels = "auto", nrow = 1)
ggsave(file.path(fig.dir, "fig2_presid_comb2.pdf"), presid_comb_p2, width = 8, height = 4)

# Fig 4: change in fruit size =========
# Generate histograms of fruit size changee
fruitsizechange <- read.csv(file.path(res.dir, "tdwgFruitSizeChange.csv")) # units are in cm since that is the original fruit length units
fruitsizechange <- subset(fruitsizechange, !is.na(changeInMedFruitSize))
medFSchange_hist <- ggplot(data = fruitsizechange) + geom_histogram(aes(changeInMedFruitSize), bins = 10) +
  labs(x = "Projected difference in\nmedian fruit length(cm)", y = "Count") +
  geom_vline(xintercept = mean(fruitsizechange$changeInMedFruitSize, na.rm = T), size = 1, linetype = "dashed", color = "coral")

maxFSchange_hist <- ggplot(data = fruitsizechange) + geom_histogram(aes(changeInMaxFruitSize), bins = 10) +
  labs(x = "Projected difference in\nmaximum fruit length (cm)", y = "Count") +
  geom_vline(xintercept = mean(fruitsizechange$changeInMaxFruitSize, na.rm = T),
             size = 1, linetype = "dashed", color = wes_palette("Cavalcanti1", 5)[4])

## Import photos
library(png)
smallfruitphoto <- readPNG(file.path(data.dir,"Arecaceae; Calyptrocalyx sp.; Baker 571_William J. Baker.png"))
largefruitphoto <- readPNG(file.path(data.dir,"Lemurophoenix halleuxii.png"))
smallfruitGrob <- grid::rasterGrob(smallfruitphoto)
largefruitGrob <- grid::rasterGrob(largefruitphoto)

## Geographic distribution of change in frugivore mammal body size
medFSchange_zCuts <- quantile(fruitsizechange$changeInMedFruitSize, probs = seq(0,1,length.out = 6), na.rm = TRUE)
fruitsizechange$medFSchange_categ <-  cut(fruitsizechange$changeInMedFruitSize, breaks = medFSchange_zCuts, include.lowest = T, ordered_result = T)
levels(fruitsizechange$medFSchange_categ) <- cleanCuts(levels(fruitsizechange$medFSchange_categ))
medFSchangemap <- ggplot() +
  geom_polygon(aes(y = lat, x = long, group = group), data = tdwg_shp2) +
  geom_point(aes(y = LAT, x= LONG, size = changeInMedFruitSize/2, color = medFSchange_categ), data = fruitsizechange) + 
  map_theme +
  coord_fixed() +
  #ggtitle(label = "Projected change in median fruit size") +
  guides(size = FALSE, colour = guide_legend(title = "Change in fruit length (cm)")) +
  theme(legend.box = 'vertical', legend.title = element_text(size = 10), legend.text = element_text(size = 10)) +
  annotation_custom(smallfruitGrob, xmin = -Inf, xmax = Inf, ymin = Inf, ymax = Inf) +
  scale_color_manual(values = brewer.pal(9,"RdBu")[c(9,7,6,3,1)])
ggsave(file.path(fig.dir, "medFSchange_map.pdf"), medFSchangemap, width = 9, height = 6, device = cairo_pdf)

maxFSchange_zCuts <- quantile(fruitsizechange$changeInMaxFruitSize, probs = seq(0,1,length.out = 6), na.rm = TRUE)
fruitsizechange$maxFSchange_categ <-  cut(fruitsizechange$changeInMaxFruitSize, breaks = maxFSchange_zCuts, include.lowest = T, ordered_result = T)
levels(fruitsizechange$maxFSchange_categ) <- cleanCuts(levels(fruitsizechange$maxFSchange_categ))

maxFSchangemap <- ggplot() +
  geom_polygon(aes(y = lat, x = long, group = group), data = tdwg_shp2) +
  geom_point(aes(y = LAT, x= LONG, size = changeInMaxFruitSize/2, color = maxFSchange_categ), data = fruitsizechange) + 
  #ggtitle(label = "Projected change in maximum (95th percentile) fruit size") +
  guides(size = FALSE, colour = guide_legend(title = "Change in fruit length (cm)")) +
  coord_fixed() +
  map_theme +
  theme(legend.box = 'vertical', legend.title = element_text(size = 10), legend.text = element_text(size = 10)) +
  scale_color_manual(values = brewer.pal(9,"RdBu")[c(9,7,6,3,1)])
ggsave(file.path(fig.dir, "maxFSchange_map.pdf"), maxFSchangemap, width = 9, height = 6, device = cairo_pdf)

medFSchange_comb <- ggdraw() +
  draw_plot(medFSchangemap + theme(legend.position = "none"), 0,0, 1,1) +
  draw_plot(smallfruitGrob, x = 0.14, y = 0.3, 0.12, 0.18) +
  draw_plot(medFSchange_hist + theme(text = element_text(size = 8),
                                     axis.text = element_text(size = 6)), 0.05, 0.1, 0.2, 0.4)
maxFSchange_comb <- ggdraw() +
  draw_plot(maxFSchangemap + theme(legend.position = "none"), 0,0, 1,1) +
  draw_plot(largefruitGrob, x = 0.14, y = 0.3, 0.12, 0.18) +
  draw_plot(maxFSchange_hist + theme(text = element_text(size = 8),
                                     axis.text = element_text(size = 6)), 0.05, 0.1, 0.2, 0.4)

medFSchange_comb_wleg <- plot_grid(medFSchange_comb, get_legend(medFSchangemap), nrow = 2, rel_heights = c(1,0.1))
maxFSchange_comb_wleg <- plot_grid(maxFSchange_comb, get_legend(maxFSchangemap), nrow = 2, rel_heights = c(1,0.1))
FSchange_comb_wleg <- plot_grid(medFSchange_comb_wleg, maxFSchange_comb_wleg, nrow = 2, labels = "auto")

ggsave(file.path(fig.dir, "fig4_FSchange.pdf"), FSchange_comb_wleg, device = cairo_pdf, width = 9, height = 10)


## Supplementary Figure 1: Plot sensitivity of gap filling =========
medFS_cor <- cor(tdwg_final$medianFruitLengthFilled, tdwg_final$medianFruitLength, method = "pearson")
medFS_gap <- ggplot() + 
  geom_point(aes(y = medianFruitLengthFilled, x = medianFruitLength), pch = 21, col = "grey40", data = tdwg_final) +
  geom_abline(aes(intercept = 0, slope = 1), col = "grey20") +
  labs(y = "Median fruit length\nwith gap filling (cm)", x = "Median fruit length\nwith gaps excluded (cm)") +
  annotate("text", x = -Inf, y = Inf, label = paste0("Pearson correlation = ", round(medFS_cor,3)), hjust = -0.1, vjust = 1.5)

maxFS_cor <- cor(tdwg_final$max95FruitLengthFilled, tdwg_final$max95FruitLength, method = "pearson")
maxFS_gap <- ggplot() +
  geom_point(aes(y = max95FruitLengthFilled, x = max95FruitLength), pch = 21, col = "grey40", data = tdwg_final) +
  geom_abline(aes(intercept = 0, slope = 1), col = "grey20") + 
  labs(y = "Maximum fruit length\nwith gap filling (cm)", x = "Maximum fruit length\nwith gaps excluded (cm)") + 
  annotate("text", x = -Inf, y = Inf, label = paste0("Pearson correlation = ", round(maxFS_cor,3)), hjust = -0.1, vjust = 1.5)

supp_fig1 <- plot_grid(medFS_gap, maxFS_gap, labels = "auto")
ggsave(supp_fig1, filename = file.path(fig.dir, "suppfig1_gapfill.pdf"), height = 4, width = 8)


## Test Fig 1.2: histogram of 
tdwg_BS_hist <- melt(tdwg_final, id.vars = "LEVEL_3_CO",
     measure.vars = c("presNat_medianBodySize", "curr_medianBodySize", "futr_medianBodySize"))

ggplot(tdwg_BS_hist, aes(log10(value), fill = variable)) + 
  geom_density(alpha = 0.5, aes(y = ..density..), position = 'identity')

tdwg_final$
tdwg_final$presNat_medianBodySize

target_col <- c("SpecName", "Mass.g", "IUCN.Status.1.2")
presnat_mamm <- read.csv(file.path(res.dir,"mammal_presNat_occ_trait.csv"))
test <- subset(presnat_mamm, LEVEL_3_CO %in% tdwg_final$LEVEL_3_CO)[target_col]
test2 <- test[!duplicated(test),]

test2$Status <- NA
test2$Status[test2$IUCN.Status.1.2 %in% c("DD", "LC", "NT")] <- "Survivors"
test2$Status[test2$IUCN.Status.1.2 %in% c("EW","EX","EP")] <- "Pleistocene + Recent extinctions"
test2$Status[test2$IUCN.Status.1.2 %in% c("EN", "VU", "CR")] <- "Possible future"

test3 <- ggplot(aes(fill = Status, x = log10(Mass.g)), data = test2) + geom_histogram(binwidth = 0.3) #+ facet_wrap(~Scenario, nrow = 3)
ggsave(test3, filename = file.path(fig.dir, "figXX_BS_time.pdf"))
                                           