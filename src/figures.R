## Generate figures
# Author: Jun Ying Lim
# Reference: Lim, J.Y., Svenning, J.-C., Göldel, B., Faurby, S. & Kissling, W.D. Past and future extinctions shape the body size - fruit size relationship between palms and mammalian frugivores.

## Packages ========================
rm(list = ls())
library(ggplot2)
library(reshape2)
library(cowplot)
library(wesanderson); library(RColorBrewer); library(viridis)
library(scales)
library(stringr); library(scales)
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
changeSciNot <- function(n) {
  # courtesy of https://stackoverflow.com/questions/29785555/in-r-using-scientific-notation-10-rather-than-e
  output <- format(n, scientific = TRUE) #Transforms the number into scientific notation even if small
  output <- sub("e", "*10^", output) #Replace e with 10^
  output <- sub("\\+0?", "", output) #Remove + symbol and leading zeros on expoent, if > 1
  output <- sub("-0?", "-", output) #Leaves - symbol but removes leading zeros on expoent, if < 1
  output
}

cleanCuts <- function(x){
  # Cleans up vector of factor levels generated using `cut` function for nicer plotting
  x <- gsub(",", " - ",x)
  x <- gsub("\\[|\\(|\\]|\\)", "", x)
  num <- str_split_fixed(x, pattern = " - ", n=2)
  new <- rep(NA, length(x))
  for(i in 1:length(x)){
    if(abs(as.numeric(num[i,1])) < 0.001){
      #start <- changeSciNot(scales::scientific(as.numeric(num[i,1]), 0))
      start <- 0 #scales::scientific(as.numeric(num[i,1]), 1)
    } else {
      start <- round(as.numeric(num[i,1]), 3 )
    }  
    if(abs(as.numeric(num[i,2])) < 0.001){
      #end <- changeSciNot(scales::scientific(as.numeric(num[i,2]), 0))
      end <- 0 #scales::scientific(as.numeric(num[i,2]), 1)
    } else {
      end <- round(as.numeric(num[i,2]), 3 )
    }
    new[i] <- paste(start,end, sep = " - ")# em-dash
    #new[i] <- paste(start,end, sep = " \u2013 ")# em-dash
  }
  return(new)
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
#tdwg_env <- read.csv(file.path(data.dir, "TDWG_Environment_AllData_2014Dec.csv"))
tdwg_env <- read.csv(file.path(data.dir, "TDWG_Environment_AllData_2019Feb.csv"))

## Import tdwg mammal and palm data ========================
tdwg_res <- read.csv(file.path(res.dir, "tdwg_mammal.csv"))

## Merge data ========================
tdwg_final <- read.csv(file.path(data.dir,"tdwg_final.csv"))


## Plotting utilities ========================
# Define map theme
map_theme <- theme(panel.background = element_blank(),
                   panel.grid = element_blank(),
                   panel.border = element_blank(),
                   axis.text = element_blank(),
                   axis.line = element_blank(),
                   axis.ticks = element_blank(),
                   axis.title = element_blank(),
                   legend.position = "bottom",
                   legend.justification = "center",
                   legend.title.align = 0,
                   strip.background = element_blank(),
                   legend.key = element_blank(),
                   strip.text = element_text(size = 18))

## Plot biogeography panels =========================
tdwg_shp_raw <- rgdal::readOGR("~/Dropbox/Projects/2019/palms/data/TDWG/level3/level3.shp")
tdwg_shp_raw@data$id <- rownames(tdwg_shp_raw@data)
tdwg_shp_small <- rgeos::gSimplify(tdwg_shp_raw, tol = 0.1, topologyPreserve = T)

tdwg_shp_small_fort <- fortify(tdwg_shp_small, region = "id")
tdwg_shp_small_fort2 <- merge(tdwg_shp_small_fort, tdwg_shp_raw@data, by = "id")
tdwg_shp_small_fort2 <- subset(tdwg_shp_small_fort2, !LEVEL_3_CO == "ANT")

neot_countrylist <- subset(tdwg_final, THREEREALM == "NewWorld")$LEVEL_3_CO
afrot_countrylist <- subset(tdwg_final, THREEREALM == "OWWest")$LEVEL_3_CO
indot_countrylist <- subset(tdwg_final, THREEREALM == "OWEast")$LEVEL_3_CO

tdwg_shp_small_fort2$Global <- ifelse(tdwg_shp_small_fort2$LEVEL_3_CO %in% tdwg_final$LEVEL_3_CO, 1, 0)
tdwg_shp_small_fort2$Neotropics <- ifelse(tdwg_shp_small_fort2$LEVEL_3_CO %in% neot_countrylist, 1, 0)
tdwg_shp_small_fort2$Afrotropics <- ifelse(tdwg_shp_small_fort2$LEVEL_3_CO %in% afrot_countrylist, 1, 0)
tdwg_shp_small_fort2$Indotropics <- ifelse(tdwg_shp_small_fort2$LEVEL_3_CO %in% indot_countrylist, 1, 0)

map_globaltropics <- ggplot() +
  geom_polygon(aes(y = lat, x = long, group = group, fill = factor(Global)),
               data = tdwg_shp_small_fort2) + 
  scale_fill_manual(values = c("grey60", "grey20")) +
  map_theme + theme(legend.position = "none",
                    plot.background = element_rect(fill = "transparent", color = NA))
ggsave(map_globaltropics,
       filename = file.path(fig.dir, "map_globaltropics.pdf"), height = 6, width = 14, bg = "transparent")
ggsave(map_globaltropics,
       filename = file.path(fig.dir, "map_globaltropics.jpg"), height = 1.5, width = 3.5, bg = "transparent")

map_neotropics <- ggplot() +
  geom_polygon(aes(y = lat, x = long, group = group, fill = factor(Neotropics)),
               data = tdwg_shp_small_fort2) + 
  scale_fill_manual(values = c("grey60", "grey20")) +
  map_theme + theme(legend.position = "none",
                    plot.background = element_rect(fill = "transparent", color = NA))
ggsave(map_neotropics, filename = file.path(fig.dir, "map_neotropics.pdf"), height = 6, width = 14)
ggsave(map_neotropics, filename = file.path(fig.dir, "map_neotropics.jpg"), height = 1.5, width = 3.5)

map_afrotropics <- ggplot() +
  geom_polygon(aes(y = lat, x = long, group = group, fill = factor(Afrotropics)),
                                    data = tdwg_shp_small_fort2) + 
  scale_fill_manual(values = c("grey60", "grey20")) +
  map_theme + theme(legend.position = "none",
                    plot.background = element_rect(fill = "transparent", color = NA))
ggsave(map_afrotropics, filename = file.path(fig.dir, "map_afrotropics.pdf"), height = 6, width = 14)
ggsave(map_afrotropics, filename = file.path(fig.dir, "map_afrotropics.jpg"), height = 1.5, width = 3.5)

map_indotropics <- ggplot() +
  geom_polygon(aes(y = lat, x = long, group = group, fill = factor(Indotropics)), data = tdwg_shp_small_fort2) + 
  scale_fill_manual(values = c("grey60", "grey20")) +
  map_theme + theme(legend.position = "none",
                    plot.background = element_rect(fill = "transparent", color = NA))
ggsave(map_indotropics, filename = file.path(fig.dir, "map_indotropics.pdf"), height = 6, width = 14)
ggsave(map_indotropics, filename = file.path(fig.dir, "map_indotropics.jpg"), height = 1.5, width = 3.5)

# Fruit size mapped  ========================
target.col <- c("LAT", "LONG", "LEVEL_3_CO", "max95FruitLengthFilled")
tdwg_maxFS <- melt(tdwg_final[,target.col], id.vars = c("LAT", "LONG", "LEVEL_3_CO"), measure.vars = "max95FruitLengthFilled")
q_probs = seq(0.0, 1.0, 0.125)
q_values = as.numeric(quantile(tdwg_maxFS$value, probs = q_probs, na.rm = T, type = 8))
max_q_mdpt <- vector()
for(i in 1:(length(q_values)-1)){ max_q_mdpt[i] <- (q_values[i] + q_values[i+1])/2 }
tdwg_maxFS$maxFS_q <- cut(tdwg_maxFS$value, q_values, include.lowest = T)
levels(tdwg_maxFS$maxFS_q) <- gsub(levels(tdwg_maxFS$maxFS_q), pattern = ",", replacement = " - ")
levels(tdwg_maxFS$maxFS_q) <- gsub(levels(tdwg_maxFS$maxFS_q), pattern = "\\(|\\[|\\]", replacement = "")

maxFS_p <- ggplot() +
  geom_polygon(aes(y = lat, x = long, group = group), fill = "grey80", size = 0.1, colour = "white",
               data = tdwg_shp2) +
  geom_point(aes(x = LONG, y = LAT,
                 fill = maxFS_q, size = maxFS_q),
             data =tdwg_maxFS, alpha = 0.9, colour = "black", pch = 21) +
  coord_fixed() +
  map_theme + theme(plot.title = element_blank()) +
  guides(size=FALSE, fill = guide_legend(title = "Maximum\n(95th percentile)\nfruit length (cm)",
                                           override.aes = list(size=sqrt(max_q_mdpt)))) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_size_manual(values = sqrt(max_q_mdpt)) +
  scale_fill_viridis(discrete = TRUE) +
  theme(legend.position = "right")

ggsave(file.path(fig.dir, "maxFS.pdf"), maxFS_p, width = 9, height = 5)

# Body size mapped  ========================
target.col <- c("LAT", "LONG", "LEVEL_3_CO", "curr_medianBodySize", "presNat_medianBodySize","curr_max95BodySize","presNat_max95BodySize")

tdwg_maxBS_melt <- melt(tdwg_final[,target.col],
                        id.vars = c("LAT", "LONG", "LEVEL_3_CO"),
                        measure.vars = c("curr_max95BodySize",
                                         "presNat_max95BodySize"))
tdwg_maxBS_melt$variable <- factor(tdwg_maxBS_melt$variable, levels = c("curr_max95BodySize", "presNat_max95BodySize"), labels = c("Current", "Present-natural"))

maxBS_arthm_min <- min(tdwg_maxBS_melt$value)/1000
maxBS_arthm_max <- max(tdwg_maxBS_melt$value)/1000

maxBS_p <- ggplot() + 
  geom_polygon(aes(y = lat, x = long, group = group), fill = "grey80", size = 0.1, colour = "white",
               data = tdwg_shp2) +
  geom_point(aes(y = LAT, x= LONG, fill = value/1000, size = log10(value)), data = tdwg_maxBS_melt, pch = 21, colour = "black") +
  facet_wrap(~variable, nrow = 2)+
  coord_fixed() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_viridis(discrete = FALSE,
                       breaks = c(maxBS_arthm_min, 0.1, 1, 10, 100, 1000, maxBS_arthm_max),
                       trans = log10_trans(),
                       labels = c(round(maxBS_arthm_min, 2), 0.1, 1, 10, 100, "1,000", "4,946")) +  
  scale_size_continuous(breaks = rev(c(3+log10(maxBS_arthm_min), 2, 3, 4, 5, 6, 3+log10(maxBS_arthm_max))),
                        labels = rev(c(0.03, 0.1, 1, 10, 100, "1,000", "4,946")),
                        range = c(0.5,3)) +
  guides(fill = guide_colorbar(title = "Maximum\n(95th percentile)\nbody mass (kg)",
                                 #label.theme = element_text(angle = 45, vjust = 0.5),
                                 order = 1, direction = "vertical"),
         size = guide_legend(title = "", order = 2)) +
  map_theme +
  theme(legend.key.width = unit(0.5, "cm"),
        plot.title = element_blank(),
        legend.spacing.x = unit(0.1, 'cm'),
        legend.position = "right",
        legend.text.align =  0.5,
        #legend.title = element_blank(),
        legend.box = "vertical")

ggsave(file.path(fig.dir, "maxBS.pdf"), maxBS_p, width = 9, height = 10)

## Model-averaged coefficients ===============
maxBS_modavg_res <- read.csv(file.path(res.dir, "maxBS_ols_modavg.csv"), stringsAsFactors = TRUE)
maxBS_modavg_res <- subset(maxBS_modavg_res, !coefficient == "(Intercept)" & 
                             Scenario %in% c("Current", "Present-Natural") &
                             Geographic.Scale %in% c("Global", "Afrotropics", "Neotropics", "Indotropics"))
maxBS_modavg_res$Variable <- factor(maxBS_modavg_res$coefficient,
                                       levels = c("curr_logMax95BS_scl", "pnat_logMax95BS_scl",
                                                  "globalPC1_scl", "globalPC2_scl", "globalPC3_scl",
                                                  "regionalPC1_scl", "regionalPC2_scl", "regionalPC3_scl",
                                                  "lgm_ens_Pano_scl", "lgm_ens_Tano_scl"),
                                       labels = c("Maximum body size", "Maximum body size",
                                                  "Climate PC1", "Climate PC2", "Climate PC3",
                                                  "Climate PC1", "Climate PC2", "Climate PC3",
                                                  "LGM Prec. anomaly", "LGM Temp. anomaly"))

maxBS_modavg_res$GeographicScale <- factor(maxBS_modavg_res$Geographic.Scale,
                                           levels = c("Global", "Afrotropics", "Neotropics", "Indotropics"),
                                           labels= c("Global", "Afrotropics", "Neotropics", "Indo-Australia"))

maxBS_modavg_plot <- ggplot(data = maxBS_modavg_res) + 
  geom_hline(aes(yintercept = 0), linetype = "dashed", colour = "grey50") +
  geom_point(aes(y = fullAvgCoef, x = Variable, colour = Scenario), position = position_dodge(0.6)) + 
  geom_errorbar(aes(ymin = fulllower2.5CI, ymax = fullupper97.5, x = Variable, colour = Scenario), position = position_dodge(0.6), width = 0.1 ) +
  scale_size_continuous(limits = c(0, 1), name = "Variable\nimportance") +
  labs(y = "Standardized\nmodel averaged coefficients") +
  facet_wrap(~ GeographicScale, drop = TRUE, nrow = 4, strip.position = "right") + 
  scale_color_manual(values = c("#FDB515", "#003262")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "grey30", fill = "transparent"))
ggsave(file.path(fig.dir, "maxBS_modavg.pdf"), maxBS_modavg_plot, height= 9, width = 8)

## Variance explained in full BS models ===============
varexp_theme <- theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
                      axis.title.y = element_text(size = 16),
                      panel.background = element_blank(),
                      strip.background = element_blank(),
                      strip.text = element_text(color = "grey20", size = 18))

maxBS_ols_modavg_res <- read.csv(file.path(res.dir, "maxBS_ols_modavg.csv"), stringsAsFactors = TRUE)
maxBS_ols_modavg_res <- subset(maxBS_ols_modavg_res, !coefficient == "(Intercept)" & 
                                 Geographic.Scale %in% c("Global", "Afrotropics", "Neotropics", "Indotropics"))
maxBS_levels <- c("curr_logMax95BS_scl", "pnat_logMax95BS_scl", "globalPC1_scl", "globalPC2_scl", "globalPC3_scl", "regionalPC1_scl", "regionalPC2_scl", "regionalPC3_scl","lgm_ens_Pano_scl", "lgm_ens_Tano_scl")
maxBS_labels <- c("Maximum body size", "Maximum body size","Climate PC1", "Climate PC2", "Climate PC3",
                  "Climate PC1", "Climate PC2", "Climate PC3","LGM Prec. Anom.", "LGM Temp. Anom.")

maxBS_ols_modavg_res$coefficient <- factor(maxBS_ols_modavg_res$coefficient,
                                           levels = maxBS_levels,
                                           labels = maxBS_labels)

maxBS_ols_modavg_res$Geographic.Scale <- factor(maxBS_ols_modavg_res$Geographic.Scale, levels = c("Global", "Afrotropics", "Neotropics", "Indotropics"), labels = c("Global", "Afrotropics", "Neotropics", "Indo-Australia"))
maxBS_ols_modavg_res$totalR2[!maxBS_ols_modavg_res$coefficient == "Maximum body size"] <- NA
maxBS_ols_modavg_res$totalR2plot <- sapply(X = maxBS_ols_modavg_res$totalR2, FUN = roundCorrectly)
maxBS_ols_modavg_res$colType <- 0
maxBS_ols_modavg_res$colType[maxBS_ols_modavg_res$coefficient == "Maximum body size" & maxBS_ols_modavg_res$Scenario == "Current"] <- 1
maxBS_ols_modavg_res$colType[maxBS_ols_modavg_res$coefficient == "Maximum body size" & maxBS_ols_modavg_res$Scenario == "Present-Natural"] <- 2

maxBS_col <- c("grey50",
               "#FDB515",
               "#003262")
maxBS_ols_modavg_plot <- ggplot(data = subset(maxBS_ols_modavg_res, !is.na(coefficient) & Scenario %in% c("Current", "Present-Natural"))) + 
  geom_bar(aes(y = varexp,  x = coefficient, fill = factor(colType)), stat = "identity") +
  facet_grid(Geographic.Scale ~ Scenario) +
  labs(y = "Variance explained", x = "") +
  geom_text(aes(y = Inf , x = Inf, label = ifelse(is.na(totalR2),"",totalR2plot)), hjust = 1, vjust = 1, parse = T) +
  guides(fill = FALSE) +
  scale_y_continuous(breaks = c(0,0.1, 0.2, 0.3)) +
  scale_fill_manual(values = maxBS_col) +
  varexp_theme
ggsave(file.path(fig.dir, "maxBS_ols_modavg_relaimpo.pdf"), maxBS_ols_modavg_plot, height= 12, width = 8)

## Partial residuals =========
maxBS_partialresid <- readRDS(file.path(res.dir, "maxBS_partialresid.rds"))
maxBS_partialresid <- maxBS_partialresid[c(1:16)]

for(i in 1:length(maxBS_partialresid)){
  if(i %in% c(1,5,9,13)){
    maxBS_partialresid[[i]]$points$Scenario <- "Current"
  }
  if(i %in% c(2,6,10,14)){
    maxBS_partialresid[[i]]$points$Scenario <- "Present-natural"
  }
  if(i %in% c(3,7,11,15)){
    maxBS_partialresid[[i]]$points$Scenario <- "Present-natural (Conservative)"
  }
  if(i %in% c(4,8,12,16)){
    maxBS_partialresid[[i]]$points$Scenario <- "Present-natural (Liberal)"
  }
  if(i %in% c(1,2,3,4)){
    maxBS_partialresid[[i]]$points$Geographic.scale <- "Global"
  }
  if(i %in% c(5,6,7,8)){
    maxBS_partialresid[[i]]$points$Geographic.scale <- "Neotropics"
  }
  if(i %in% 9:12){
    maxBS_partialresid[[i]]$points$Geographic.scale <- "Afrotropics"
  }
  if(i %in% 13:16){
    maxBS_partialresid[[i]]$points$Geographic.scale <- "Indo-Australia"
  }
}

geog_scale <- c("Global", "Afrotropics", "Neotropics", "Indo-Australia")
maxBS_partialresid_pts <- do.call("rbind", lapply(maxBS_partialresid[c(1,2,5,6,9,10,13,14)], FUN = function(x){ x$points }))
maxBS_partialresid_pts$Geographic.scale <- factor(maxBS_partialresid_pts$Geographic.scale,
                                                  levels = geog_scale)
write.csv(maxBS_partialresid_pts, file.path(res.dir, "maxBS_partialresid_pts.csv"), row.names = F)

maxBS_partialresid_slopes <- unlist(lapply(maxBS_partialresid[c(1,2,5,6,9,10,13,14)], FUN = function(x) { x$slope } ))
maxBS_partialresid_intercepts <- unlist(lapply(maxBS_partialresid[c(1,2,5,6,9,10,13,14)], FUN = function(x) { x$intercept } ))
maxBS_abline <- data.frame(slope = maxBS_partialresid_slopes,
                           intercept = maxBS_partialresid_intercepts,
                           Geographic.scale = factor(rep(c("Global", "Neotropics", "Afrotropics", "Indo-Australia"), each = 2), levels = geog_scale),
                           Scenario = rep(c("Current", "Present-natural"), 2))
write.csv(maxBS_abline, file.path(res.dir, "maxBS_abline.csv"), row.names = F)

maxBS_abline$slope[maxBS_abline$Geographic.scale == "Afrotropics"] <- NA
maxBS_abline$intercept[maxBS_abline$Geographic.scale == "Afrotropics"] <- NA
maxBS_presid_plot <- ggplot(aes(y = presid, x = x, group = Scenario, color = Scenario),
                            data = subset(maxBS_partialresid_pts, Scenario %in% c("Current", "Present-natural"))) + 
  geom_point(shape = 21) + 
  geom_abline(aes(slope = slope, intercept = intercept, colour = Scenario),
              data = subset(maxBS_abline, Scenario %in% c("Current", "Present-natural")), show.legend = F) + 
  facet_wrap(~Geographic.scale, nrow = 4, strip.position = "right") + 
  labs(x = "Log maximum body size\n(standard deviations)", y = "Log maximum fruit size \n(partial residuals)") +
  scale_colour_manual(values = c("#FDB515", "#003262")) +
  guides(color = guide_legend(override.aes=list(size = 5))) +
  theme(legend.position = "bottom",
        legend.key = element_blank(),
        legend.title = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(color = "grey30", fill = "transparent"))
ggsave(maxBS_presid_plot, filename = file.path(fig.dir, "maxBS_presid_plot.pdf"), height = 8, width = 3)

## Histograms of fruit size change =========
fruitsizechange <- read.csv(file.path(res.dir, "tdwgFruitSizeChange.csv")) # units are in cm since that is the original fruit length units

max95FSchange_hist <- ggplot(data = fruitsizechange) + geom_histogram(aes(changeInMax95FruitSize), bins = 15) +
  labs(x = "Projected difference in\nmaximum fruit length (cm)", y = "Count") +
  geom_vline(xintercept = mean(fruitsizechange$changeInMaxFruitSize, na.rm = T),
             size = 0.5, linetype = "dashed", color = "#ED4E33") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme(panel.background = element_blank())

max95FSchange_cons_hist <- ggplot(data = fruitsizechange) + geom_histogram(aes(changeInMax95FruitSize_cons), bins = 15) +
  labs(x = "Projected difference in\nmaximum fruit length (cm)", y = "Count") +
  geom_vline(xintercept = mean(fruitsizechange$changeInMaxFruitSize_cons, na.rm = T),
             size = 0.5, linetype = "dashed", color = "#ED4E33") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme(panel.background = element_blank())

maxFSchange_hist <- ggplot(data = fruitsizechange) + geom_histogram(aes(changeInMaxFruitSize), bins = 15) +
  labs(x = "Projected difference in\nmaximum fruit length (cm)", y = "Count") +
  geom_vline(xintercept = mean(fruitsizechange$changeInMaxFruitSize, na.rm = T),
             size = 0.5, linetype = "dashed", color = "#ED4E33") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme(panel.background = element_blank())

maxFSchange_cons_hist <- ggplot(data = fruitsizechange) + geom_histogram(aes(changeInMaxFruitSize_cons), bins = 15) +
  labs(x = "Projected difference in\nmaximum fruit length (cm)", y = "Count") +
  geom_vline(xintercept = mean(fruitsizechange$changeInMaxFruitSize_cons, na.rm = T),
             size = 0.5, linetype = "dashed", color = "#ED4E33") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme(panel.background = element_blank())

## Change in fruit size mapped =======
# Maximum fruit size
maxFSchange_zCuts <- quantile(fruitsizechange$changeInMaxFruitSize, probs = seq(0,1,length.out = 6), na.rm = TRUE)
fruitsizechange$maxFSchange_categ <-  cut(fruitsizechange$changeInMaxFruitSize, breaks = maxFSchange_zCuts, include.lowest = T, ordered_result = T)
maxFSchange_zmidpt <- vector()
for(i in 1:(length(maxFSchange_zCuts)-1)){ maxFSchange_zmidpt[i] <- (maxFSchange_zCuts[i] + maxFSchange_zCuts[i+1])/2 }
levels(fruitsizechange$maxFSchange_categ) <- cleanCuts(levels(fruitsizechange$maxFSchange_categ))
maxFSchange_lib_map <- ggplot() +
  geom_polygon(aes(y = lat, x = long, group = group), data = tdwg_shp2) +
  geom_point(aes(y = LAT, x= LONG, size = maxFSchange_categ, fill = maxFSchange_categ), data = fruitsizechange, pch = 21, colour = "grey80") + 
  guides(size = FALSE,
         fill = guide_legend(title.position = "top",
                               title.hjust = 0.5,
                               title = "Change in maximum fruit length (cm)",
                               override.aes = list(size = seq(1,3, length.out = 5)))) +#list(size = rescale(maxFSchange_zmidpt, to = c(1,5))))) +
  coord_fixed() +
  map_theme +
  ggtitle("High defaunation scenario") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.box = 'vertical',
        legend.key = element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10)) +
  scale_size_manual(values = seq(1,3, length.out = 5)) +#rescale(maxFSchange_zmidpt, to = c(1,5))) +
  scale_fill_manual(values = brewer.pal(9,"RdBu")[c(9,7,6,3,1)])
ggsave(file.path(fig.dir, "maxFSchange_lib_map.pdf"), maxFSchange_lib_map, width = 9, height = 6, device = "pdf")

maxFSchangecons_zCuts <- quantile(fruitsizechange$changeInMaxFruitSize_cons, probs = seq(0,1,length.out = 6), na.rm = TRUE)
fruitsizechange$maxFSchange_cons_categ <-  cut(fruitsizechange$changeInMaxFruitSize_cons, breaks = maxFSchangecons_zCuts, include.lowest = T, ordered_result = T)
maxFSchangecons_zmidpt <- vector()
for(i in 1:(length(maxFSchange_zCuts)-1)){ maxFSchangecons_zmidpt[i] <- (maxFSchangecons_zCuts[i] + maxFSchangecons_zCuts[i+1])/2 }
levels(fruitsizechange$maxFSchange_cons_categ) <- cleanCuts(levels(fruitsizechange$maxFSchange_cons_categ))
maxFSchange_cons_map <- ggplot() +
  geom_polygon(aes(y = lat, x = long, group = group), data = tdwg_shp2) +
  geom_point(aes(y = LAT, x= LONG, size = maxFSchange_cons_categ, fill = maxFSchange_cons_categ), data = fruitsizechange, pch = 21, colour = "grey80") + 
  guides(size = FALSE,
         fill = guide_legend(title = "Change in maximum fruit length (cm)",
                               title.position = "top",
                               title.hjust = 0.5,
                               override.aes = list(size = seq(1,3, length.out = 5)))) + #list(size = rescale(maxFSchangecons_zmidpt, to = c(1,5))))) +
  coord_fixed() +
  map_theme +
  ggtitle("Low defaunation scenario") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.box = 'vertical',
        legend.key = element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10)) +
  scale_size_manual(values = seq(1,3, length.out = 5)) +#rescale(maxFSchangecons_zmidpt, to = c(1,5))) +
  scale_fill_manual(values = brewer.pal(9,"RdBu")[c(9,7,6,3,1)])
ggsave(file.path(fig.dir, "maxFSchange_cons_map.pdf"), maxFSchange_cons_map, width = 9, height = 6)

maxFSchange_plot <- plot_grid(maxFSchange_lib_map,
                        maxFSchange_cons_map,
                        nrow = 2)
ggsave(file.path(fig.dir, "maxFSchange_map.pdf"), maxFSchange_plot, width = 9, height = 10)

# 95th percentile maps
max95FSchange_zCuts <- quantile(fruitsizechange$changeInMax95FruitSize, probs = seq(0,1,length.out = 6), na.rm = TRUE)
fruitsizechange$max95FSchange_categ <-  cut(fruitsizechange$changeInMax95FruitSize,
                                            breaks = max95FSchange_zCuts, include.lowest = T, ordered_result = T)
max95FSchange_zmidpt <- vector()
for(i in 1:(length(max95FSchange_zCuts)-1)){ max95FSchange_zmidpt[i] <- (max95FSchange_zCuts[i] + max95FSchange_zCuts[i+1])/2 }
levels(fruitsizechange$max95FSchange_categ) <- cleanCuts(levels(fruitsizechange$max95FSchange_categ))
max95FSchange_lib_map <- ggplot() +
  geom_polygon(aes(y = lat, x = long, group = group), data = tdwg_shp2) +
  geom_point(aes(y = LAT, x= LONG, size = max95FSchange_categ, fill = max95FSchange_categ), data = fruitsizechange, pch = 21, colour = "grey80") + 
  guides(size = FALSE,
         fill = guide_legend(title.position = "top",
                               title.hjust = 0.5,
                               title = "Change in maximum (95th percentile) fruit length (cm)",
                               override.aes = list(size = seq(1,3, length.out = 5)))) + #list(size = rescale(max95FSchange_zmidpt, to = c(1,5))))) +
  coord_fixed() +
  map_theme +
  ggtitle("High defaunation scenario") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.box = 'vertical',
        legend.key = element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10)) +
  scale_size_manual(values = seq(1,3, length.out = 5)) +#rescale(max95FSchange_zmidpt, to = c(1,5))) +
  scale_fill_manual(values = brewer.pal(9,"RdBu")[c(9,7,6,3,1)])
ggsave(file.path(fig.dir, "max95FSchange_lib_map.pdf"), max95FSchange_lib_map, width = 9, height = 6)

max95FSchangecons_zCuts <- quantile(fruitsizechange$changeInMax95FruitSize_cons, probs = seq(0,1,length.out = 6), na.rm = TRUE)
fruitsizechange$max95FSchange_cons_categ <-  cut(fruitsizechange$changeInMax95FruitSize_cons,
                                                 breaks = max95FSchangecons_zCuts, include.lowest = T, ordered_result = T)
max95FSchangecons_zmidpt <- vector()
for(i in 1:(length(max95FSchangecons_zCuts)-1)){ max95FSchangecons_zmidpt[i] <- (max95FSchangecons_zCuts[i] + max95FSchangecons_zCuts[i+1])/2 }
levels(fruitsizechange$max95FSchange_cons_categ) <- cleanCuts(levels(fruitsizechange$max95FSchange_cons_categ))
max95FSchange_cons_map <- ggplot() +
  geom_polygon(aes(y = lat, x = long, group = group), data = tdwg_shp2) +
  geom_point(aes(y = LAT, x= LONG, size = max95FSchange_cons_categ, fill = max95FSchange_cons_categ), data = fruitsizechange, pch = 21, colour = "grey80") + 
  guides(size = FALSE,
         fill = guide_legend(title = "Change in maximum (95th percentile) fruit length (cm)",
                               title.position = "top",
                               title.hjust = 0.5,
                               override.aes = list(size = seq(1,3, length.out = 5)))) +#list(size = rescale(max95FSchangecons_zmidpt, to = c(1,5))))) +
  coord_fixed() +
  map_theme +
  ggtitle("Low defaunation scenario") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.box = 'vertical',
        legend.key = element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10)) +
  scale_size_manual(values = seq(1,3, length.out = 5)) + #rescale(max95FSchangecons_zmidpt, to = c(1,5))) +
  scale_fill_manual(values = brewer.pal(9,"RdBu")[c(9,7,6,3,1)])
ggsave(file.path(fig.dir, "max95FSchange_cons_map.pdf"), max95FSchange_cons_map, width = 9, height = 6)

max95FSchange_plot <- plot_grid(max95FSchange_lib_map,
                                max95FSchange_cons_map,
                          nrow = 2)
ggsave(file.path(fig.dir, "max95FSchange_map.pdf"), max95FSchange_plot, width = 9, height = 10)

## Extinction and body size ==========
tdwg_final <- read.csv(file.path(data.dir,"tdwg_final.csv"))
tdwg_env <- read.csv(file.path(data.dir, "TDWG_Environment_AllData_2019Feb.csv"))
presnat_mamm <- read.csv(file.path(res.dir,"mammal_presNat_occ_trait.csv"))
presnat_mamm <- subset(presnat_mamm, LEVEL_3_CO %in% tdwg_final$LEVEL_3_CO)
presnat_mamm <- merge(x = presnat_mamm,  y = tdwg_env[c("LEVEL_3_CO", "REALM_LONG")], by = "LEVEL_3_CO")
presnat_mamm$REALM_LONG <- factor(presnat_mamm$REALM_LONG, levels = c("Neotropics", "Afrotropics", "IndoMalay", "Australasia"), labels = c("Neotropics", "Afrotropics", "Indo-Australia", "Indo-Australia"))

target_col <- c("SpecName", "Mass.g", "IUCN.Status.1.2", "REALM_LONG", "Default", "Liberal", "Conservative")

mamm_df  <- presnat_mamm[target_col]
mamm_df <- mamm_df[!duplicated(mamm_df),]

mamm_df$ConsStat <- ifelse(mamm_df$Default == "Y" & mamm_df$Liberal == "Y" & mamm_df$Conservative == "Y",
                           "Conservative",
                           ifelse(mamm_df$Default == "Y" & mamm_df$Liberal == "Y",
                                  "Default", "Liberal"))

mamm_df$CurrVsPnat <- ifelse(mamm_df$IUCN.Status.1.2 %in% c("DD", "LC", "NT", "EN", "VU", "CR"),
                             "Current",
                             "Late Pleistocene +\nHolocene extinctions")
mamm_df$FutrVsCurr <- ifelse(mamm_df$IUCN.Status.1.2 %in% c("NT", "DD", "LC"),
                             "Low or unknown\nextinction risk",
                             "High\nextinction risk") # make sure to subset out all extinct species (i.e., status: EW, EX, EP)
mamm_df$FutrVsCurr <- factor(mamm_df$FutrVsCurr, levels = c("Low or unknown\nextinction risk", "High\nextinction risk"))

mamm_df2 <- mamm_df
mamm_df2$REALM_LONG <- "Global"
mamm_df3 <- rbind(mamm_df, mamm_df2)
mamm_df3$REALM_LONG <- factor(mamm_df3$REALM_LONG, levels = c("Global", "Afrotropics", "Neotropics", "Indo-Australia"))
mamm_df3$ConsStat <- factor(mamm_df3$ConsStat, levels = c("Default", "Liberal", "Conservative"),
                            labels = c('Default',
                                       'Liberal',
                                       "Conservative"))
mamm_df3 <- mamm_df3[!duplicated(mamm_df3),]

CurrVsPnatBSHist <- ggplot(aes(fill = CurrVsPnat, x = log10(Mass.g/1000)), data = mamm_df3) +
  geom_histogram(binwidth = 0.35) +
  facet_wrap(~REALM_LONG, nrow = 4, scales = "free_y") +
  scale_fill_manual(values = wes_palette("Royal1", n = 2)) +
  scale_y_continuous(name = "Number of mammalian frugivores", expand = c(0,0)) +
  labs(x=expression(paste(Log[10], " body mass (kg)", sep = " ")))+
  scale_x_continuous(expand = c(0,0)) +
  theme(axis.title = element_text(size = 12),
        strip.background = element_blank(),
        strip.text = element_text(size = 12, colour = "grey20"),
        panel.background = element_blank(), 
        legend.position = "bottom", 
        legend.title = element_blank() )
ggsave(CurrVsPnatBSHist, filename = file.path(fig.dir, "CurrVsPnatBSHist.pdf"), width = 5, height = 9)

CurrVsFutrBSHist <- ggplot(aes(fill = FutrVsCurr, x = log10(Mass.g/1000)), data = subset(mamm_df3, !IUCN.Status.1.2 %in% c("EX", "EP", "EW"))) +
  geom_histogram(binwidth = 0.5) +
  facet_wrap(~REALM_LONG, nrow = 4, scales = "free") +
  scale_fill_manual(values = wes_palette("Royal1", n = 2)) +
  scale_y_continuous(name = "Number of\nmammalian frugivores", expand = c(0,0)) +
  labs(x=expression(paste(Log[10], " body mass (kg)", sep = " ")))+
  scale_x_continuous(expand = c(0,0)) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12, colour = "grey20"),
        panel.background = element_blank(), 
        legend.position = "bottom", 
        legend.title = element_blank())
ggsave(CurrVsFutrBSHist, filename = file.path(fig.dir, "CurrVsFutrBSHist.pdf"), width = 5, height = 9)

write.csv(mamm_df3, file.path(res.dir, "mamm_df.csv"), row.names = T)

## FIG 1: Maximum fruit size, body size and changes in body size ==========
fig1 <- plot_grid(maxFS_p + theme(legend.title = element_text(size = 14),
                                  legend.text = element_text(size = 11)),
                  plot_grid(maxBS_p + theme(legend.title =  element_text(size = 14),
                                            legend.key.height = unit(0.8, "cm"),
                                            legend.text = element_text(size = 11)),
                                            #legend.title.align = 0.5),
                            CurrVsPnatBSHist,
                            nrow = 1, rel_widths = c(0.75, 0.25), labels = c("b","c"), label_size = 20),
                  nrow = 2, rel_heights = c(0.8, 1.2), labels = "a", label_size = 20)
ggsave(fig1, filename = file.path(fig.dir,"fig1_prelim.pdf"), height = 11, width = 12)

## FIG 2: Partial residual + model coeff ==========
presid_leg <- get_legend(maxBS_presid_plot + 
                         guides(color = guide_legend(override.aes = list(shape = "square", size = 7))) +
                         theme(legend.position = "bottom",
                               legend.text = element_text(size = 10),
                                 legend.key = element_blank(),
                                 legend.box.background = element_blank()))
fig2 <- plot_grid(plot_grid(maxBS_modavg_plot +
                            theme(legend.position = "none",
                            legend.background = element_rect(fill = "transparent"),
                            legend.key = element_blank(),
                            legend.key.width = unit(0.5,"cm"),
                            legend.title = element_blank(),
                            panel.grid = element_blank(),
                            plot.background = element_rect(fill = "transparent", color = NA)),
                    maxBS_presid_plot + 
                      theme(legend.position = "none",
                            plot.background = element_rect(fill = "transparent", color = NA),
                            panel.grid = element_blank()),
                            axis = "b", align = "h",
                            rel_widths = c(1.2,0.6),
                            labels = "auto", nrow = 1),
                  presid_leg, nrow = 2, rel_heights = c(1, 0.1))

ggsave(file.path(fig.dir, "fig2_presid_modcoeff.pdf"), fig2, width = 7, height = 7)

## FIG 3: Variance explained ===============
fig3 <- maxBS_ols_modavg_plot + theme(axis.title.y = element_text(size = 26),
                                axis.text.y = element_text(size = 16),
                                axis.text.x = element_text(size = 16, angle = 50),
                                plot.background = element_rect(fill = "transparent", color = NA),
                                panel.grid = element_blank())
ggsave(file.path(fig.dir, "fig3_varexp.pdf"), fig3, width = 10, height = 10)

## FIG 4: Change in fruit size ===============
max95FSchange_comb <- ggdraw() +
  draw_plot(max95FSchange_lib_map +
              #guides(colour = guide_legend(title.position = "top", title = "Change in maximum (95th percentile) fruit length (cm)")) +
              theme(legend.position = "bottom",
                                   legend.title = element_text(size = 8),
                                   legend.text = element_text(size = 6)), 0, 0, 1,1) + 
  draw_plot(max95FSchange_hist + theme(text = element_text(size = 6)), 0.1, 0.2, 0.18, 0.25)

max95FSchangecons_comb <- ggdraw() +
  draw_plot(max95FSchange_cons_map + 
              theme(legend.position = "bottom",
                                        legend.title = element_text(size = 8),
                                        legend.text = element_text(size = 6)), 0,0, 1,1) +
  draw_plot(max95FSchange_cons_hist + theme(text = element_text(size = 6)), 0.1, 0.2, 0.18, 0.25)


fig4 <- plot_grid(CurrVsFutrBSHist + theme(title = element_text(size = 8),
                                           axis.title = element_text(size = 8),
                                           axis.text = element_text(size = 6),
                                           legend.text = element_text(size = 8),
                                           legend.position = "bottom"), 
                  plot_grid(max95FSchange_comb,
                            max95FSchangecons_comb,
                            rel_heights = c(1,1), nrow = 2,
                            labels = c("b", "c")),
                  nrow = 1, rel_widths = c(0.4, 1.2), scale = c(1,1),
                  labels = c("a"))

ggsave(file.path(fig.dir, "fig4_max95FSchange.pdf"), fig4, width = 10, height = 8)

## SUPP FIG 1: Plot sensitivity of gap filling =========
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

## SUPP FIG 2: Effect of removing Australia ====================
maxBS_ols_modavg_res <- read.csv(file.path(res.dir, "maxBS_ols_modavg.csv"), stringsAsFactors = TRUE)
medBS_ols_modavg_res <- read.csv(file.path(res.dir, "medBS_ols_modavg.csv"), stringsAsFactors = TRUE)

maxBS_ols_modavg_res$Model <- "Maximum Body Size"
medBS_ols_modavg_res$Model <- "Median Body Size"

ols_modavg_res <- rbind(maxBS_ols_modavg_res, medBS_ols_modavg_res)
ols_modavg_res_subset <- subset(ols_modavg_res, Geographic.Scale %in% c("Indotropics", "Indotropics (No Australia + New Guinea)") &
                                  !coefficient == "(Intercept)" &
                                  Scenario %in% c("Current", "Present-Natural"))
coef_levels <- c("curr_logMedBS_scl", "pnat_logMedBS_scl","curr_logMax95BS_scl", "pnat_logMax95BS_scl",
                 "regionalPC1_noAus_scl", "regionalPC2_noAus_scl", "regionalPC3_noAus_scl", "regionalPC1_scl", "regionalPC2_scl", "regionalPC3_scl","lgm_ens_Pano_scl", "lgm_ens_Tano_scl")
coef_labels <- c("Median body size", "Median body size", "Maximum body size", "Maximum body size","Climate PC1", "Climate PC2", "Climate PC3",
                 "Climate PC1", "Climate PC2", "Climate PC3","LGM Prec. Anom.", "LGM Temp. Anom.")

ols_modavg_res_subset$Geographic.Scale <- factor(ols_modavg_res_subset$Geographic.Scale, 
                                                 levels = c("Indotropics", "Indotropics (No Australia + New Guinea)"),
                                                 labels = c("Indo-Australia", "Indotropics\n(Indo-Australia but Australia +\nNew Guinea excluded)"))

ols_modavg_res_subset$Variable <- factor(ols_modavg_res_subset$coefficient,
                                         levels = coef_levels,
                                         labels = coef_labels)

noAus_plot <- 
  ggplot(data = ols_modavg_res_subset) + 
  geom_hline(aes(yintercept = 0), linetype = "dashed", colour = "grey50") +
  geom_point(aes(y = fullAvgCoef, x = Variable, colour = Scenario, shape = Geographic.Scale), position = position_dodge(0.6)) + 
  geom_errorbar(aes(ymin = fulllower2.5CI, ymax = fullupper97.5, x = Variable, colour = Scenario, shape = Geographic.Scale), position = position_dodge(0.6), width = 0.1 ) +
  scale_size_continuous(limits = c(0, 1), name = "Variable\nimportance") +
  labs(y = "Standardized\nmodel averaged coefficients") +
  facet_wrap(~Model, drop = TRUE, nrow = 4, strip.position = "right", scales = "free_x") + 
  scale_color_manual(values = c("#FDB515", "#003262")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        strip.text = element_blank(),
        strip.background = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "grey30", fill = "transparent"))
ggsave(noAus_plot, filename = file.path(fig.dir, "indotropics_sensitivity.pdf"), height = 5, width = 10)

## SUPP FIG 3: Correlation plots for variables =========
tdwg_global_corr <- read.csv(file.path(res.dir, "tdwg_global_corr.csv"), row.names = 1, check.names = F)
tdwg_nw_corr <- read.csv(file.path(res.dir, "tdwg_nw_corr.csv"), row.names = 1, check.names = F)
tdwg_oww_corr <- read.csv(file.path(res.dir, "tdwg_oww_corr.csv"), row.names = 1, check.names = F)
tdwg_owe_corr <- read.csv(file.path(res.dir, "tdwg_owe_corr.csv"), row.names = 1, check.names = F)

createCorrHeatmap <- function(cor){
  # Takes a correlation matrix and plots a prettyheatmap  
  cor_df <- as.data.frame(cor)
  cor_df$Variable <- rownames(cor)
  cor_df_melt <- melt(cor_df)
  cor_df_melt$value_bin <- cut(cor_df_melt$value, breaks = seq(-1,1,0.25))
  qlabs  <- c("-1 \u2013 -0.75",
              "-0.75 \u2013 -0.5",
              " -0.5 \u2013 -0.25",
              "-0.25 \u2013 0   ",
              "0 \u2013 0.25",
              "0.25 \u2013 0.5 ",
              "0.5 \u2013 0.75",
              "0.75 \u2013 1")
  levels(cor_df_melt$value_bin) <- qlabs
  cor_df_melt$value_abs <- abs(cor_df_melt$value)
  cor_df_melt$Variable <- factor(cor_df_melt$Variable, levels  = rev(rownames(cor)))
  cor_df_melt$variable <- factor(cor_df_melt$variable, levels  = rownames(cor))
  ggplot(cor_df_melt) + 
    geom_point(aes(y= Variable, x= variable, colour = value_bin, size =value_abs), pch = 16) + 
    geom_text(aes(y= Variable, x= variable, label = round(value,2)), size =2) +
    # guide scales by area and not by radius
    guides(colour = guide_legend(title = "Correlation",override.aes = list(size = 2+c(8.75, 6.75, 3.75, 1.25,
                                                                                      1.25, 3.75, 6.75, 8.75)*8/10)),
           size = F) +
    scale_radius(limits = c(0, 1), range = c(2,10)) +
    scale_colour_manual(values = brewer.pal("RdYlBu", n= 8),
                        limits = qlabs) + 
    coord_fixed() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title = element_blank(),
          panel.background = element_blank(),
          legend.position = "right")
}


corr_comb_plot <- plot_grid( plotlist = list(
  createCorrHeatmap(tdwg_global_corr) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + ggtitle("Global"),
  createCorrHeatmap(tdwg_nw_corr) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + ggtitle("Neotropics"),
  createCorrHeatmap(tdwg_oww_corr) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + ggtitle("Afrotropics"),
  createCorrHeatmap(tdwg_owe_corr) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + ggtitle("Indo-Australia")
), nrow = 2)
corr_comb_plot_wleg <- plot_grid(corr_comb_plot, 
                                 get_legend(createCorrHeatmap(tdwg_global_corr) + 
                                              theme(legend.key = element_blank(),
                                                    legend.text.align = 0.5)), nrow = 1, rel_widths = c(0.8, 0.1))

ggsave(file.path(fig.dir, "corr_comb_plot_wleg.pdf"), corr_comb_plot_wleg, width = 13, height = 9, device = cairo_pdf)

## SUPP FIG 4: Body size distribution of extinct mammal frugivores  =========
FrugClassPlot <-
  ggplot(aes(fill = ConsStat, x = log10(Mass.g/1000)), data = subset(mamm_df3, IUCN.Status.1.2 == "EP")) +
  geom_histogram(binwidth = 0.5) +
  facet_wrap(~REALM_LONG, nrow = 4, scales = "free_y") +
  scale_fill_manual(values = brewer.pal(name = "YlGnBu", n = 9)[c(3,5,7)]) + #wes_palette("Royal1", n = 3)) +
  scale_y_continuous(name = "Number of mammalian frugivores", expand = c(0,0)) +
  labs(x=expression(paste(Log[10], " body mass (kg)", sep = " "))) +
  scale_x_continuous(expand = c(0,0)) +
  theme(axis.title = element_text(size = 12),
        strip.background = element_blank(),
        strip.text = element_text(size = 12, colour = "grey20"),
        panel.background = element_blank(),
        legend.position = "bottom", 
        legend.title = element_blank())
ggsave(FrugClassPlot, filename = file.path(fig.dir, "FrugClassPlot.pdf"), width = 5, height = 9)

## SUPP FIG 5: Parametric bootstrap for CTMC models  =========
# Plot empirical rate versus rates derived from simulations using the same rates
dimarco_parboot <- readRDS(file.path(res.dir, "dimarco_parboot.rds"))
hoffmann_parboot <- readRDS(file.path(res.dir, "hoffmann_parboot.rds"))

dimarco_par_sim <- dimarco_parboot[[1]]
dimarco_par_obs <- dimarco_parboot[[2]]

hoffmann_par_sim <- hoffmann_parboot[[1]]
hoffmann_par_obs <- hoffmann_parboot[[2]]

dimarco_parboot_plot <- ggplot() +
  geom_histogram(aes(x = value), bins = 30, data = dimarco_par_sim) +
  geom_vline(aes(xintercept = value), data = dimarco_par_obs) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  labs(y = "No. of simulations", x = "Instantaneous rates", title = "Di Marco et al. 2014") +
  facet_wrap(~par_label) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

hoffmann_parboot_plot <- ggplot() +
  geom_histogram(aes(x = value), bins = 30, data = hoffmann_par_sim) +
  geom_vline(aes(xintercept = value), data = hoffmann_par_obs) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  labs(y = "No. of simulations", x = "Instantaneous rates", title = "Hoffmann et al. 2010") +
  facet_wrap(~par_label) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

parbootplot <- cowplot::plot_grid(hoffmann_parboot_plot, dimarco_parboot_plot, labels = "auto")

ggsave(parbootplot, filename = file.path(fig.dir, "parbootplot.pdf"), width = 12, height = 8)

## SUPP FIG 6: Plot true maximum fruit size change =========
maxFSchange_comb <- ggdraw() +
  draw_plot(maxFSchange_lib_map +
              theme(legend.position = "bottom",
                    legend.title = element_text(size = 8),
                    legend.text = element_text(size = 6)), 0, 0, 1,1) + 
  draw_plot(maxFSchange_hist + theme(text = element_text(size = 6)), 0.1, 0.2, 0.18, 0.25)

maxFSchangecons_comb <- ggdraw() +
  draw_plot(maxFSchange_cons_map + 
              theme(legend.position = "bottom",
                    legend.title = element_text(size = 8),
                    legend.text = element_text(size = 6)), 0,0, 1,1) +
  draw_plot(maxFSchange_cons_hist + theme(text = element_text(size = 6)), 0.1, 0.2, 0.18, 0.25)

maxFS_combined_plot <- plot_grid(maxFSchange_comb,
                                 maxFSchangecons_comb,
                                 rel_heights = c(1,1), nrow = 2)
ggsave(maxFS_combined_plot, filename = file.path(fig.dir, "maxFS_combined_plot.pdf"), width = 9, height = 10)

## Neighbourhood definition  ===============
nb_def <- readRDS(file.path(res.dir, "nb_def.rds"))
names(nb_def[[1]]) <- 1:length(nb_def[[1]])

nb_df <- ldply(.data = nb_def[[1]],.id = "from",
      .fun = function(x){ data.frame("to" = x)})

nb_coords <- data.frame(nb_def[[2]])
names(nb_coords) <- c("LONG", "LAT")
nb_coords$id <- 1:nrow(nb_coords)

nb_df2 <- merge(x = nb_df, y = nb_coords, by.x = "from", by.y = "id", all = T)
nb_df3 <- merge(x = nb_df2, y = nb_coords, by.x = "to", by.y = "id", suffixes = c("_from", "_to"), all = T)

nb_plot <- ggplot() + 
  geom_polygon(aes(x= long, y = lat, group = group), data = tdwg_shp2, fill = "grey80", colour = "white", size = 0.1) +
  geom_point(aes(x = LONG, y = LAT), data = nb_coords, pch = 21, size = 3 ) +
  geom_segment(aes(x = LONG_from, y = LAT_from, xend = LONG_to, yend = LAT_to ), data = nb_df3, size = 0.3) +
  theme_map()
ggsave(nb_plot, filename = file.path(fig.dir, "nb_def.pdf"), width=14, height = 6)

## OTHER PLOTS  ===============
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


test <- ggplot(data = fruitsizechange) + 
  geom_point(aes(y = changeInMaxFruitSize_cons, x = (curr_maxBodySize - futr_cons_maxBodySize)/1000)) +
  geom_abline(aes(slope = 0, intercept = 0), linetype = "dashed", colour = "red") +
  labs(y = "Change in fruit size\n(lower extinction prob. estimates)", x = "Decrease in true maximum body size (kg)\n(lower extinction prob. estimates)")
test2 <- ggplot(data = fruitsizechange) + 
  geom_point(aes(y = changeInMaxFruitSize_cons, x = (curr_max95BodySize - futr_cons_max95BodySize)/1000)) +
  geom_abline(aes(slope = 0, intercept = 0), linetype = "dashed", colour = "red") +
  labs(y = "Change in fruit size\n(lower extinction prob. estimates)", x = "Decrease in 95th percentile body size (kg)\n(lower extinction prob. estimates)")

ggsave(plot_grid(test, test2, nrow = 1), filename = "~/Desktop/fruitsizechange.pdf")

ggplot(aes(y = (curr_maxBodySize - futr_cons_maxBodySize)/1000, 
           x = (curr_max95BodySize - futr_cons_max95BodySize)/1000 ), data = fruitsizechange) + geom_point() +
geom_smooth(method = "lm")


mean(tdwg_final$presNat_max95BodySize) / mean(tdwg_final$curr_max95BodySize)
mean(tdwg_final$presNat_max95BodySize_cons) / mean(tdwg_final$curr_max95BodySize)
mean(tdwg_final$presNat_max95BodySize_cons)
list.files(res.dir)
curr_occ_trait <- read.csv(file.path(res.dir, "mammal_curr_occ_trait.csv"))
test <- subset(curr_occ_trait, LEVEL_3_CO %in% subset(tdwg_final, REALM_LONG == "Neotropics")$LEVEL_3_CO)

test <- read.csv(file.path(res.dir, "mamm_df.csv"))
test_nw <- subset(test, REALM_LONG == "Neotropics" & CurrVsPnat == "Late Pleistocene +\nHolocene extinctions")
nrow(subset(test_nw, Mass.g > 1000000 & ConsStat == "Default"))

