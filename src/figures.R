
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
tdwg_shape <- readOGR(file.path(rawdata.dir, "TDWG", "level3", "level3.shp"))
tdwg_shape@data$id <- rownames(tdwg_shape@data)
tdwg_shape_df <- fortify(tdwg_shape, region = "id")
tdwg_shape_df <- merge(tdwg_shape_df, tdwg_shape@data)

## Import tdwg environmental data ========================
tdwg_env <- read.csv(file.path(data.dir, "TDWG_Environment_AllData_2014Dec.csv"))

## Import tdwg mammal and palm data ========================
tdwg_res <- read.csv(file.path(res.dir, "tdwg_mammal.csv"))

## Merge data ========================
tdwg_plot_df <- merge(tdwg_res, tdwg_env, by = "LEVEL_3_CO", all.x = TRUE)

## Plotting utilities ========================
# Define map theme
map_theme <- theme(axis.ticks = element_blank(), axis.title = element_blank(), axis.text = element_blank())


## Plotting fruit data ========================
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
  map_theme + 
  theme(legend.position = "bottom")
ggsave(filename = file.path(fig.dir, "fruitSizeFilledPlot.pdf"), gappedFruitLengthPlot,
       width = 14, height = 7)

# FruitLengthFilledvsNotPlot <- ggplot() +
#   geom_point(aes(y = log(meanFruitLength), x = log(meanFruitLengthFilled)), data = tdwg_plot_df) +
#   geom_abline(aes(intercept = 0, slope = 1))
# temp <- resid(lm(log(meanFruitLength) ~ log(meanFruitLengthFilled), data=tdwg_plot_df))
# 
# FruitLengthFilledvsNotPlotResid <- ggplot(data.frame("resid" = temp[which(tdwg_plot_df$meanFruitLength != tdwg_plot_df$meanFruitLengthFilled)])) +
#  geom_histogram(aes(x = resid), bins = 15)
# FruitLengthSensitivity<- grid.arrange(FruitLengthFilledvsNotPlot, FruitLengthFilledvsNotPlotResid, ncol =2)
# ggsave(filename = file.path(fig.dir, "FruitLengthSensitivity.pdf"), FruitLengthSensitivity)

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

## Plotting mammal data  ========================
target.col <- c("LAT", "LONG", "LEVEL_3_CO", "curr_meanBodySize", "curr_nSp", "presNat_meanBodySize","presNat_nSp","curr_medianBodySize","presNat_medianBodySize")
tdwg_plot_melt <- melt(tdwg_plot_df[,target.col],
                       id.vars = c("LAT", "LONG", "LEVEL_3_CO"),
                       measure.vars = c("curr_meanBodySize",
                                        "presNat_meanBodySize",
                                        "curr_medianBodySize",
                                        "presNat_medianBodySize",
                                        "curr_nSp", "presNat_nSp"))
#tdwg_plot_melt$variable <- factor(tdwg_plot_melt$variable, levels = c("curr_medianBodySize", "presNat_medianBodySize"), labels = c("Current", "Present-natural"))

# Plot mammal body size world wide
q_probs = seq(0.0, 1.0, 0.1)
q_labs <- levels(cut(log(tdwg_plot_melt$value), quantile(log(tdwg_plot_melt$value), probs = q_probs, na.rm = T)))
q_labs <- gsub(q_labs, pattern = "\\(|\\]", replacement = "")
q_labs <- gsub(q_labs, pattern = ",", replacement = "-")
q_labs <- c("0 - 17", "17 - 40", "40 - 64", "64 - 114", "114 - 254", "254 - 652", "652 - 3,498", "3,498 - 20,537", "20,537 - 120,572", "120,572 - 1,329,083")

tdwg_plot_melt$logvalue <- log(tdwg_plot_melt$value)

mammalbodySizePlot <- ggplot() +
  geom_polygon(aes(y = lat, x = long, group = group),
               data = subset(tdwg_shape_df, !LEVEL_3_CO == "ANT")) +
  geom_point(aes(x = LONG, y = LAT, fill = cut(logvalue, quantile(logvalue, probs = q_probs), labels = q_labs), size = logvalue), colour = "white", shape = 21,
             data = subset(tdwg_plot_melt, variable %in% c("Current", "Present-natural") & !is.na(logvalue)), alpha = 0.9) +
  facet_wrap(~variable, nrow = 2, scales = "free")+
  theme(panel.background = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), legend.position = "bottom",legend.justification = "center", strip.background = element_blank(), strip.text = element_text(size = 18)) +
  guides(size=FALSE, fill = guide_legend(override.aes = list(size=10))) +
  scale_fill_manual(values = colorRamps::matlab.like(n = 10), name = "Median\nmammal body size (g)", breaks = q_labs)
ggsave(mammalbodySizePlot, filename = file.path(fig.dir, "currmedianLogBodySizeCombined.pdf"), width = 10, height = 10)

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



 # Plot mammal species richness
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

# Plot palm species richness
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
