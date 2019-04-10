## Packages ========================
library(ggplot2)
library(viridis)
library(reshape2)
library(cowplot)

## Directories ========================
main.dir <- "/Users/junyinglim/Dropbox/Projects/2019/palms/projects/megafaunalFrugivore"
data.dir <- file.path(main.dir, "data")
res.dir <- file.path(main.dir, "results")
src.dir <- file.path(main.dir, "src")
fig.dir <- file.path(main.dir, "figs")
rawdata.dir <- "~/Dropbox/Projects/2019/palms/data"
frug.dir <- file.path(rawdata.dir, "frugivores")
options(stringsAsFactors =FALSE)

## Import tdwg polygon data ========================
tdwg_shp_raw <- readOGR("~/Dropbox/Projects/2019/palms/data/TDWG/level3/level3.shp")
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

# Plot fruit sizes
target.col <- c("LAT", "LONG", "LEVEL_3_CO", "medianFruitLengthFilled")
tdwg_final_melt <- melt(tdwg_final[,target.col], id.vars = c("LAT", "LONG", "LEVEL_3_CO"), measure.vars = "medianFruitLengthFilled")

q_probs = seq(0.0, 1.0, 0.1)
q_labs <- levels(cut(tdwg_final_melt$value, quantile(tdwg_final_melt$value, probs = q_probs, na.rm = T)))
q_labs <- gsub(q_labs, pattern = "\\(|\\]", replacement = "")
q_labs <- gsub(q_labs, pattern = ",", replacement = " - ")

fruitSizePlot <- ggplot() +
  geom_polygon(aes(y = lat, x = long, group = group),
               data = tdwg_shp2) +
  geom_point(aes(x = LONG, y = LAT, fill = cut(value, quantile(value, probs = q_probs),
                                               labels = q_labs), size = value),
             colour = "white", shape = 21, data = tdwg_final_melt, alpha = 0.9) +
  mammal_theme() +
  guides(size=FALSE, fill = guide_legend(override.aes = list(size=10))) +
  scale_fill_manual(values = colorRamps::matlab.like(n = 10), name = "Median\nfruit length (cm)", breaks = q_labs)

ggsave("~/Desktop/fruitSizes.pdf", fruitSizePlot, width = 9, height = 5)

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
target.col <- c("LAT", "LONG", "LEVEL_3_CO", "curr_medianBodySize", "presNat_meanBodySize","curr_medianBodySize","presNat_medianBodySize", "curr_nSp", "presNat_nSp")
tdwg_final_melt <- melt(tdwg_final[,target.col],
                        id.vars = c("LAT", "LONG", "LEVEL_3_CO"),
                        measure.vars = c("curr_nSp",
                                         "presNat_nSp"))
tdwg_final_melt$variable <- factor(tdwg_final_melt$variable, levels = c("curr_nSp", "presNat_nSp"), labels = c("Current", "Present-natural"))

q_probs = seq(0.0, 1.0, 0.1)
q_labs <- levels(cut(log10(tdwg_final_melt$value), quantile(log10(tdwg_final_melt$value), probs = q_probs, na.rm = T)))
q_labs <- gsub(q_labs, pattern = "\\(|\\]", replacement = "")
q_labs <- gsub(q_labs, pattern = ",", replacement = "-")
library(stringr)
q_labs_label <- str_split_fixed(q_labs, "-", n = 2)
q_labs_label2 <- paste0(round(10^as.numeric(q_labs_label[,1]), 0), " - ", round(10^as.numeric(q_labs_label[,2]), 0))

tdwg_final_melt$logvalue <- log10(tdwg_final_melt$value)

mammalbodySizePlot <- ggplot() +
  geom_polygon(aes(y = lat, x = long, group = group),
               data = subset(tdwg_shp, !LEVEL_3_CO == "ANT")) +
  geom_point(aes(x = LONG, y = LAT, fill = cut(logvalue, quantile(logvalue, probs = q_probs),
                                               labels = q_labs_label2), size = logvalue),
             colour = "white", shape = 21, data = subset(tdwg_final_melt, variable %in% c("Current", "Present-natural")), alpha = 0.9) +
  facet_wrap(~variable, nrow = 2, scales = "free")+
  mammal_theme() +
  guides(size=FALSE, fill = guide_legend(override.aes = list(size=10))) +
  scale_fill_manual(values = colorRamps::matlab.like(n = 10), name = "Species richness", breaks = q_labs_label2)
#"Mean\nmammal body size (g)"
ggsave(mammalbodySizePlot, filename = "~/Desktop/mammalNsp.pdf", width = 10, height = 10)
deltabodySizePlot <- ggplot() +
  geom_polygon(aes(y = lat, x = long, group = group),
               data = subset(tdwg_shape_df, !LEVEL_3_CO == "ANT")) +
  geom_point(aes(x = LONG, y = LAT, fill = log10(presNat_medianBodySize)-log10(curr_medianBodySize),
                 size = log10(presNat_medianBodySize)-log10(curr_medianBodySize)), colour = "white", shape = 21, #size =3,
             data = tdwg_plot_df, alpha = 0.9) +
  theme(panel.background = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), legend.position = "bottom",legend.justification = "center", strip.background = element_blank(), strip.text = element_text(size = 18)) +
  guides(size=FALSE) +
  scale_fill_viridis(name = "Difference in\nlog median\nmammal body size (g)")
ggsave(deltabodySizePlot, filename = file.path(fig.dir, "deltaBodySize.pdf"), width = 10, height = 5)

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


## Plotting median fruit size and area
asd <- ggplot(data = tdwg_final) + geom_point(aes(y = log(medianFruitLengthFilled), x = log(AREA_KM2), colour = factor(ISISLAND)))
ggsave(asd, filename = "~/Desktop/asd.pdf")



summary(lm(log(medianFruitLengthFilled) ~ log(AREA_KM2) * factor(ISISLAND), data = tdwg_final))

rangesize_melt <- melt(tdwg_final,
                    id.vars = c("LAT", "LONG", "LEVEL_3_CO", "THREEREALM"),
                    measure.vars = c("curr_medRangeSize", "presNat_medRangeSize"))

rangesizePlot <- ggplot() + geom_polygon(aes(y = lat, x = long, group = group), data = tdwg_shp2) + geom_point(aes(y = LAT, x = LONG, fill = scale(value), size = scale(value)), shape = 21, colour = "white", data = rangesize_melt) + facet_wrap(~variable) + map_theme + scale_fill_viridis() +  guides(size = FALSE)
ggsave("~/Desktop/asda.pdf", rangesizePlot, width= 14, height = 7)

ggplot(data = tdwg_final) + geom_point(aes(y = log(medianFruitLengthFilled), x = presNat_medRangeSize, color = THREEREALM))

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

## Plot relative importance of body size dispersion ===============
dispBS_modavg_res <- read.csv(file.path(res.dir, "dispBS_modavg_res.csv"), stringsAsFactors = TRUE)

dispBS_modavg_res$Variable <- factor(dispBS_modavg_res$Variable,
                                        levels = c("curr_dispBodySize_scl", "presNat_dispBodySize_scl",
                                                   "globalPC1_scl", "globalPC2_scl", "globalPC3_scl",
                                                   "regionalPC1_scl", "regionalPC2_scl", "regionalPC3_scl",
                                                   "lgm_ens_Pano_scl", "lgm_ens_Tano_scl", "soilcount_scl"),
                                        labels = c("Body size dispersion", "Body size dispersion",
                                                   "PC1", "PC2", "PC3",
                                                   "PC1", "PC2", "PC3",
                                                   "LGM Prec. anomaly", "LGM Temp. anomaly", "Soil diversity"))

dispBS_modavg_res$GeographicScale <- factor(dispBS_modavg_res$GeographicScale, levels = c("Global", "Afrotropics", "Neotropics", "Indotropics"))

dispBS_modavg_plot <- ggplot(data = dispBS_modavg_res) + 
  geom_hline(aes(yintercept = 0), linetype = "dashed", colour = "grey50") +
  geom_point(aes(y = avgStdCoef, x = Variable, colour = Scenario, size = var_impt), position = position_dodge(0.8)) + 
  geom_errorbar(aes(ymin = avgStdCIlower, ymax = avgStdCIupper, x = Variable, colour = Scenario), position = position_dodge(0.8), width = 0.1 ) +
  scale_size_continuous(limits = c(0, 1), name = "Variable\nimportance") +
  labs(title = "Body size dispersion\n(scaled)", y = "Standardized\nmodel averaged coefficients", x= "Variables") +
  facet_wrap(~ GeographicScale, drop = TRUE, nrow = 4) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(fig.dir, "dispBS_modavg.pdf"), dispBS_modavg_plot, height= 9, width = 8)


#modavg_comb_plot <- plot_grid(plotlist = list(medBS_modavg_plot, maxBS_modavg_plot, dispBS_modavg_plot), ncol = 3)
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


## Plotting partial residuals
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
fruitsizechange <- read.csv(file.path(res.dir, "tdwgFruitSizeChange.csv")) # units are in cm since that is the original fruit length units
medFSchange_hist <- ggplot(data = fruitsizechange) + geom_histogram(aes(changeInMedFruitSize), bins = 15) +
  labs(x = "Projected difference in\nmedian fruit length(cm)", y = "Count") +
  geom_vline(xintercept = mean(fruitsizechange$changeInMedFruitSize, na.rm = T), size = 2, linetype = "dashed", color = "coral")

maxFSchange_hist <- ggplot(data = fruitsizechange) + geom_histogram(aes(changeInMaxFruitSize), bins = 15) +
  labs(x = "Projected difference in\nmaximum fruit length (cm)", y = "Count") +
  geom_vline(xintercept = mean(fruitsizechange$changeInMaxFruitSize, na.rm = T), size = 2, linetype = "dashed", color = wes_palette("Cavalcanti1", 5)[4])

FSchange_histcomb <- plot_grid(medFSchange_hist, maxFSchange_hist, labels = "AUTO", nrow = 2)
ggsave(file.path(fig.dir, "FSchange_histcomb.pdf"), FSchange_histcomb, width = 5, height = 8)

## Geographic distribution of change in frugivore mammal body size
library(scales)
zCuts <- quantile(fruitsizechange$changeInMedFruitSize, probs = seq(0,1,length.out = 8), na.rm = TRUE)

fruitsizechange$medFSchange_categ <-  cut(fruitsizechange$changeInMedFruitSize, breaks = zCuts, include.lowest = T, ordered_result = T)
FSchangemap <- ggplot() +
  geom_polygon(aes(y = lat, x = long, group = group), data = tdwg_shp2) +
  geom_point(aes(y = LAT, x= LONG, size = changeInMedFruitSize, color = medFSchange_categ), data = fruitsizechange) + 
  map_theme +
  theme(legend.box = 'vertical') +
  scale_color_manual(values = brewer.pal(9,"RdBu")[c(9,8,7,6,4,3,2,1)])

ggsave(file.path(fig.dir, "FSchangemap.pdf"), FSchangemap, width = 9, height = 6)


range(fruitsizechange$changeInMedFruitSize, na.rm = T)
levels(fruitsizechange$medFSchange_categ)
library(RColorBrewer)
fruitsizechange$changeInMedFruitSize[109]
subset(fruitsizechange, THREEREALM == "NewWorld" & changeInMedFruitSize > 0.2)

mammal_curr_occ_trait <- read.csv(file.path(res.dir, "mammal_curr_occ_trait.csv"))
subset(mammal_curr_occ_trait, LEVEL_3_CO == "SWC")
