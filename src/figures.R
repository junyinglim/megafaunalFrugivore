## Palm-frugivores
# Generate figures for manuscript and for data exploration
# Author: Jun Ying Lim

## Packages ========================
library(ggplot2)
library(reshape2)
library(cowplot)
library(wesanderson); library(RColorBrewer); library(viridis)
library(scales)
library(stringr)

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



## Plotting fruit data ========================

# Plot median fruit sizes
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
  map_theme + theme(plot.title = element_text(size = 20)) +
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
  map_theme + theme(plot.title = element_text(size = 20))+
  guides(size=FALSE, colour = guide_legend(title = "Fruit length (cm)", override.aes = list(size=5))) +
  scale_color_viridis(discrete = TRUE)

ggsave(file.path(fig.dir, "maxFS.pdf"), maxFS_p, width = 9, height = 5)

## Plotting biogeography data ========================
threerealmplot <- ggplot() +
  geom_polygon(aes(y = lat, x= long, group = group), data = subset(tdwg_shape_df, !LEVEL_3_CO == "ANT")) +
  geom_point(aes(x = LONG, y = LAT, fill = THREEREALM), data = tdwg_plot_df, shape = 21, size = 3, color = "white")+
  scale_fill_viridis(discrete = TRUE) +
  theme(legend.position = "bottom")
ggsave(threerealmplot, filename = file.path(fig.dir, "threerealmPlot.pdf"), width = 14, height = 7 )

islandplot <- ggplot() +
    geom_polygon(aes(y = lat, x= long, group = group), data = subset(tdwg_shape_df, !LEVEL_3_CO == "ANT")) +
    geom_point(aes(x = LONG, y = LAT, fill = Geology_3categ), data = subset(tdwg_plot_df, !is.na(Geology_3categ)), shape = 21, size = 3, color = "white")+
  scale_fill_viridis(discrete = TRUE) +
  theme(legend.position = "bottom")
ggsave(islandplot, filename = file.path(fig.dir, "islandPlot.pdf"), width = 14, height = 7 )

## Plotting mammal data on map ========================
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
  geom_point(aes(y = LAT, x= LONG, colour = (value/1000), size = (value/1000)), data = tdwg_medBS_melt) +
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
  theme(plot.title = element_text(size = 20), legend.key.width = unit(1, "cm"))
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
  geom_point(aes(y = LAT, x= LONG, colour = value/1000, size = value/1000), data = tdwg_maxBS_melt) +
  facet_wrap(~variable, nrow = 2)+
  coord_fixed() +
  scale_colour_viridis(discrete = FALSE,
                       breaks = c(maxBS_arthm_min, 0.1, 1, 10, 100, 1000, maxBS_arthm_max),
                       trans = log10_trans(),
                       labels = c(round(maxBS_arthm_min, 2), 0.1, 1, 10, 100, 1000, round(maxBS_arthm_max))) +  
  ggtitle("Maximum (95th percentile) body size") +
  guides(size = FALSE, colour = guide_colorbar(title = "Body mass (kg)",
                                               label.theme = element_text(angle = 45, vjust = 0.5))) +
  map_theme +
  theme(plot.title = element_text(size = 20), legend.key.width = unit(1, "cm"))

ggsave(file.path(fig.dir, "maxBS.pdf"), maxBS_p, width = 9, height = 10)

# Fig 1: Combine fruit size and body size ===========
fig1_label_size = 20
FS_comb <- plot_grid(medFS_p + theme(legend.position = "none"),
                     maxFS_p + theme(legend.position = "none"),
                     nrow = 1, labels = c("A","B"), label_size = fig1_label_size)

FS_comb_leg <- plot_grid(get_legend(medFS_p), get_legend(maxFS_p), nrow = 1)
FS_comb_wleg <- plot_grid(FS_comb, FS_comb_leg, rel_heights = c(1, 0.2), nrow = 2)

BS_comb <- plot_grid(medBS_p + theme(legend.position = "none"),
                     maxBS_p + theme(legend.position = "none"), 
                     nrow = 1, labels= c("C", "D"), label_size = fig1_label_size)

BS_comb_leg <- plot_grid(get_legend(medBS_p), get_legend(maxBS_p), nrow = 1)

fig1_BSFS <- plot_grid(FS_comb, FS_comb_leg, BS_comb, BS_comb_leg,
                       nrow = 4, rel_heights = c(1,0.2,2,0.2))

ggsave(file.path(fig.dir, "fig1_prelim.pdf"), fig1_BSFS, width = 18, height = 15)


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



## Plot median fruit size and porportion of megafauna
#propMegaMam_curr
medBS_v_FS_plot <- ggplot(aes(y = log(medianFruitLengthFilled), x = propMegaMam_presnat, color = THREEREALM),
                          data = tdwg_final) + 
  geom_point() + 
#  facet_wrap(~variable, nrow = 2) +
  geom_smooth(method = "lm") +
  geom_text_repel(aes(label= LEVEL_3_CO), size = 1, colour = "grey20",segment.color = "grey20", segment.size = 0.3 ) +
  labs(y = "Median fruit size\n(scaled)", x = "Proportion of megafaunal") +
  theme(legend.position = "bottom")


## Plot mammal body size against  fruit length ===============
medBS_v_FS_melt <- melt(tdwg_final,
                        id.vars = c("LAT", "LONG", "LEVEL_3_CO", "THREEREALM", "medianFruitLengthFilled", "megaSpLoss"),
                        measure.vars = c("curr_medianBodySize", "presNat_medianBodySize"))
medBS_v_FS_melt$variable <- factor(medBS_v_FS_melt$variable,
                                   levels = c("curr_medianBodySize", "presNat_medianBodySize"),
                                   labels = c("Current", "Present-natural"))
medBS_v_FS_plot <- ggplot(aes(y = log(medianFruitLengthFilled), x = scale(log(value)), color = THREEREALM, size = megaSpLoss),
                          data = medBS_v_FS_melt) + 
  geom_point() + 
  facet_wrap(~variable, nrow = 2) +
  geom_smooth(method = "lm") +
  geom_text_repel(aes(label= LEVEL_3_CO), size = 1, colour = "grey20",segment.color = "grey20", segment.size = 0.3 ) +
  labs(x = "Median body size\n(scaled)", y = "Median fruit size\n(scaled)", title = "Median body size\nscaled") +
  theme(legend.position = "bottom")
ggsave(file.path(fig.dir, "medBS_v_FS.pdf"), medBS_v_FS_plot, width = 6, height = 8)

maxBS_v_FS_melt <- melt(tdwg_final,
                        id.vars = c("LAT", "LONG", "LEVEL_3_CO", "THREEREALM", "max95FruitLengthFilled", "megaSpLoss"),
                        measure.vars = c("curr_max95BodySize", "presNat_max95BodySize"))
maxBS_v_FS_melt$variable <- factor(maxBS_v_FS_melt$variable,
                                   levels = c("curr_max95BodySize", "presNat_max95BodySize"),
                                   labels = c("Current", "Present-natural"))
maxBS_v_FS_plot <- ggplot(aes(y = log(max95FruitLengthFilled), x = scale(log(value)), color = THREEREALM, size = megaSpLoss),
                          data = maxBS_v_FS_melt) + 
  labs(x = "Body size 95th percentile\n(scaled)", y = "Fruit size 95th percentile\n(scaled)", title = "Body size 95th percentile\nscaled") +
  geom_point() + 
  geom_text_repel(aes(label= LEVEL_3_CO), size = 1, colour = "grey20",segment.color = "grey20", segment.size = 0.3 ) +
  facet_wrap(~variable, nrow = 2) +
  geom_smooth(method = "lm") +
  theme(legend.position = "bottom")
ggsave(file.path(fig.dir, "maxBS_v_FS.pdf"), maxBS_v_FS_plot, width = 6, height = 8)

bs_v_fs_comb_plot <- plot_grid(plotlist = list(medBS_v_FS_plot, maxBS_v_FS_plot), ncol = 2)
ggsave(bs_v_fs_comb_plot, filename = file.path(fig.dir, "bs_v_fs_comb.pdf"), height = 8, width = 12)

# dispBS_v_FS_melt <- melt(tdwg_final,
#                         id.vars = c("LAT", "LONG", "LEVEL_3_CO", "THREEREALM", "dispFruitLengthFilled"),
#                         measure.vars = c("curr_dispBodySize", "presNat_dispBodySize"))
# dispBS_v_FS_melt$variable <- factor(dispBS_v_FS_melt$variable,
#                                     levels = c("curr_dispBodySize", "presNat_dispBodySize"),
#                                     labels = c("Current", "Present-natural"))
# dispBS_v_FS_plot <- ggplot(aes(y = dispFruitLengthFilled, x = scale(value), color = THREEREALM),
#                           data = dispBS_v_FS_melt) + 
#   labs(x = "Body size dispersion\n(scaled)", y = "Fruit size dispersion\n(scaled)") +
#   geom_point() + 
#   facet_wrap(THREEREALM~variable, nrow = 3) +
#   geom_text_repel(aes(label= LEVEL_3_CO), size = 1, colour = "grey20",segment.color = "grey20", segment.size = 0.3 ) +
#   geom_smooth(method = "lm") +
#   point_theme
# ggsave(file.path(fig.dir, "dispBS_v_FS.pdf"), dispBS_v_FS_plot, width = 6, height = 8)



## Plot relative importance of maximum body size ===============
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

## Plot OLS model averaging results for median body size ===============
medBS_modavg_res <- read.csv(file.path(res.dir, "medBS_modavg_res.csv"), stringsAsFactors = TRUE)
medBS_modavg_res$Variable <- factor(medBS_modavg_res$Variable,
                                       levels = c("curr_logMedBS_scl", "pnat_logMedBS_scl",
                                                  "globalPC1_scl", "globalPC2_scl", "globalPC3_scl",
                                                  "regionalPC1_scl", "regionalPC2_scl", "regionalPC3_scl",
                                                  "lgm_ens_Pano_scl", "lgm_ens_Tano_scl", "soilcount_scl"),
                                       labels = c("Median body size", "Median body size",
                                                  "PC1", "PC2", "PC3",
                                                  "PC1", "PC2", "PC3",
                                                  "LGM Prec. anomaly", "LGM Temp. anomaly", "Soil diversity"))

medBS_modavg_res$GeographicScale <- factor(medBS_modavg_res$GeographicScale, levels = c("Global", "Afrotropics", "Neotropics", "Indotropics"))

medBS_modavg_plot <- ggplot(data = medBS_modavg_res) + 
  geom_hline(aes(yintercept = 0), linetype = "dashed", colour = "grey50") +
  geom_point(aes(y = avgStdCoef, x = Variable, colour = Scenario, size = var_impt), position = position_dodge(0.8)) + 
  geom_errorbar(aes(ymin = avgStdCIlower, ymax = avgStdCIupper, x = Variable, colour = Scenario), position = position_dodge(0.8), width = 0.1 ) +
  scale_size_continuous(limits = c(0, 1), name = "Variable\nimportance") +
  labs(title = "Median body size\n(scaled)", y = "Standardized\nmodel averaged coefficients", x= "Variables") +
  facet_wrap(~ GeographicScale, drop = TRUE, nrow = 4) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(fig.dir, "medBS_modavg.pdf"), medBS_modavg_plot, height= 9, width = 8)


# Combined OLS OLS model averaging results
modavg_comb_plot <- plot_grid(plotlist = list(medBS_modavg_plot, maxBS_modavg_plot), ncol = 2)
ggsave(modavg_comb_plot, filename = file.path(fig.dir, "modavg_comb.pdf"), height = 10, width = 15 )





## Plot SAR model averaging results for median body size
sar_knear_medBS_modavg_res <- read.csv(file.path(res.dir, "sar_knear_medBS_modavg.csv"))
sar_soi_medBS_modavg_res <- read.csv(file.path(res.dir, "sar_soi_medBS_modavg.csv"))
sar_medBS_modavg_res <- rbind(sar_knear_medBS_modavg_res, sar_soi_medBS_modavg_res)

sar_medBS_modavg_res$coefficient <- factor(sar_medBS_modavg_res$coefficient,
                                    levels = c("curr_logMedBS_scl", "pnat_logMedBS_scl",
                                               "globalPC1_scl", "globalPC2_scl", "globalPC3_scl",
                                               "regionalPC1_scl", "regionalPC2_scl", "regionalPC3_scl",
                                               "lgm_ens_Pano_scl", "lgm_ens_Tano_scl"),
                                    labels = c("Log median body size", "Log median body size",
                                               "PC1", "PC2", "PC3",
                                               "PC1", "PC2", "PC3",
                                               "LGM Prec. anomaly", "LGM Temp. anomaly"))





medBS_sarmodavg_plot <- ggplot(data = sar_medBS_modavg_res) + 
  geom_hline(aes(yintercept = 0), linetype = "dashed", colour = "grey50") +
  geom_point(aes(y = fullAvgCoef, x = coefficient, colour = Scenario, size = importance), position = position_dodge(0.8)) + 
  geom_errorbar(aes(ymin = lower2.5, ymax = upper97.5, x = coefficient, colour = Scenario), position = position_dodge(0.8), width = 0.1 ) +
  scale_size_continuous(limits = c(0, 1), name = "Variable\nimportance") +
  facet_wrap(Method~Weighting) +
  labs(title = "Median body size", y = "Standardized\nmodel averaged coefficients", x= "Variables") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(medBS_sarmodavg_plot, filename = file.path(fig.dir, "sar_medBS_modavg.pdf"), height = 8, width = 14)


## Plotting partial residuals of FS against median and maximum BS =========
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
medBS_comb_presid$Scenario <- rep(c("Current", "Present-natural"), each = 129)

med_comb_presid_p <- ggplot(data= medBS_comb_presid) +
  geom_point(aes(y = presid, x = x, color = Scenario), shape = 21) + labs(x = "Log median body size\n(standard deviations)", y = "Log median fruit size \n(partial residuals)") + 
  geom_abline(slope = medBS_curr_presid$slope, intercept = medBS_curr_presid$intercept, colour = curr_col) +
  geom_abline(slope = medBS_pnat_presid$slope, intercept = medBS_pnat_presid$intercept, colour = pnat_col) +
  comb_theme +
  scale_color_manual(values = c(curr_col, pnat_col))

presid_comb_p2 <- plot_grid(med_comb_presid_p, max_comb_presid_p, labels = "AUTO", nrow = 1)
ggsave(file.path(fig.dir, "presid_comb2.pdf"), presid_comb_p2, width = 10, height = 5)

## Histograms of changes in frugivore mammal body size
# Notes: One country ROD (Rodriguez) only has one frugivore that is also threatened and so won't have any estimates 
fruitsizechange <- read.csv(file.path(res.dir, "tdwgFruitSizeChange.csv")) # units are in cm since that is the original fruit length units
fruitsizechange <- subset(fruitsizechange, !is.na(changeInMedFruitSize))
medFSchange_hist <- ggplot(data = fruitsizechange) + geom_histogram(aes(changeInMedFruitSize), bins = 15) +
  labs(x = "Projected difference in\nmedian fruit length(cm)", y = "Count") +
  geom_vline(xintercept = mean(fruitsizechange$changeInMedFruitSize, na.rm = T), size = 1, linetype = "dashed", color = "coral")

maxFSchange_hist <- ggplot(data = fruitsizechange) + geom_histogram(aes(changeInMaxFruitSize), bins = 10) +
  labs(x = "Projected difference in\nmaximum fruit length (cm)", y = "Count") +
  geom_vline(xintercept = mean(fruitsizechange$changeInMaxFruitSize, na.rm = T),
             size = 1, linetype = "dashed", color = wes_palette("Cavalcanti1", 5)[4])

#FSchange_histcomb <- plot_grid(medFSchange_hist, maxFSchange_hist, labels = "AUTO", nrow = 2)
#ggsave(file.path(fig.dir, "FSchange_histcomb.pdf"), FSchange_histcomb, width = 5, height = 8)

## Geographic distribution of change in frugivore mammal body size
medFSchange_zCuts <- quantile(fruitsizechange$changeInMedFruitSize, probs = seq(0,1,length.out = 6), na.rm = TRUE)
fruitsizechange$medFSchange_categ <-  cut(fruitsizechange$changeInMedFruitSize, breaks = medFSchange_zCuts, include.lowest = T, ordered_result = T)
levels(fruitsizechange$medFSchange_categ) <- cleanCuts(levels(fruitsizechange$medFSchange_categ))
medFSchangemap <- ggplot() +
  geom_polygon(aes(y = lat, x = long, group = group), data = tdwg_shp2) +
  geom_point(aes(y = LAT, x= LONG, size = changeInMedFruitSize/2, color = medFSchange_categ), data = fruitsizechange) + 
  map_theme +
  coord_fixed() +
  ggtitle(label = "Projected change in median fruit size") +
  guides(size = FALSE, colour = guide_legend(title = "Change in fruit length (cm)")) +
  theme(legend.box = 'vertical', legend.title = element_text(size = 10), legend.text = element_text(size = 10)) +
  scale_color_manual(values = brewer.pal(9,"RdBu")[c(9,7,6,3,1)])
ggsave(file.path(fig.dir, "medFSchange_map.pdf"), medFSchangemap, width = 9, height = 6, device = cairo_pdf)

maxFSchange_zCuts <- quantile(fruitsizechange$changeInMaxFruitSize, probs = seq(0,1,length.out = 6), na.rm = TRUE)
fruitsizechange$maxFSchange_categ <-  cut(fruitsizechange$changeInMaxFruitSize, breaks = maxFSchange_zCuts, include.lowest = T, ordered_result = T)
levels(fruitsizechange$maxFSchange_categ) <- cleanCuts(levels(fruitsizechange$maxFSchange_categ))
maxFSchangemap <- ggplot() +
  geom_polygon(aes(y = lat, x = long, group = group), data = tdwg_shp2) +
  geom_point(aes(y = LAT, x= LONG, size = changeInMaxFruitSize/2, color = maxFSchange_categ), data = fruitsizechange) + 
  ggtitle(label = "Projected change in maximum (95th percentile) fruit size") +
  guides(size = FALSE, colour = guide_legend(title = "Change in fruit length (cm)")) +
  coord_fixed() +
  map_theme +
  theme(legend.box = 'vertical', legend.title = element_text(size = 10), legend.text = element_text(size = 10)) +
  scale_color_manual(values = brewer.pal(9,"RdBu")[c(9,7,6,3,1)])
ggsave(file.path(fig.dir, "maxFSchange_map.pdf"), maxFSchangemap, width = 9, height = 6, device = cairo_pdf)


# Fig 3: change in fruit size =========
medFSchange_comb <- ggdraw() +
  draw_plot(medFSchangemap + theme(legend.position = "none"), 0,0, 1,1) +
  draw_plot(medFSchange_hist + theme(text = element_text(size = 8),
                                     axis.text = element_text(size = 6)), 0.05, 0.1, 0.2, 0.4)
maxFSchange_comb <- ggdraw() +
  draw_plot(maxFSchangemap + theme(legend.position = "none"), 0,0, 1,1) +
  draw_plot(maxFSchange_hist + theme(text = element_text(size = 8),
                                     axis.text = element_text(size = 6)), 0.05, 0.1, 0.2, 0.4)

medFSchange_comb_wleg <- plot_grid(medFSchange_comb, get_legend(medFSchangemap), nrow = 2, rel_heights = c(1,0.1))
maxFSchange_comb_wleg <- plot_grid(maxFSchange_comb, get_legend(maxFSchangemap), nrow = 2, rel_heights = c(1,0.1))
FSchange_comb_wleg <- plot_grid(medFSchange_comb_wleg, maxFSchange_comb_wleg, nrow = 2, labels = "AUTO")

ggsave(file.path(fig.dir, "fig3_FSchange.pdf"), FSchange_comb_wleg, device = cairo_pdf, width = 9, height = 10)



levels(fruitsizechange$medFSchange_categ)
library(RColorBrewer)
fruitsizechange$changeInMedFruitSize[109]
subset(fruitsizechange, THREEREALM == "NewWorld" & changeInMedFruitSize > 0.2)

mammal_curr_occ_trait <- read.csv(file.path(res.dir, "mammal_curr_occ_trait.csv"))
mammal_pnat_occ_trait <- read.csv(file.path(res.dir, "mammal_presnat_occ_trait.csv"))
subset(mammal_curr_occ_trait, CONTINENT == "AUSTRALASIA")
table(mammal_curr_occ_trait$CON)

subset(tdwg_final, CONTINENT == "AUSTRALASIA")
