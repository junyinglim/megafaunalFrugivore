## Analyze Data
# Look at extinction risk among palms

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
library(rredlist)

## Import palm dataset ========================
palm_occ <- read.csv(file.path(data.dir, "palms_in_tdwg3.csv"))
palm_trait <- read.csv(file.path(data.dir, "PalmTraits_10.csv"))
#palm_trait$SpecName <- gsub(palm_trait$SpecName, pattern= " ", replacement = "_")

palm_redlist <- ddply(.data = palm_trait,
                      .variables = .(SpecName),
                      .fun = function(x) { rl_history(name = x$SpecName,
                                                      key = "548ba4d434dfc552f9945f666e2e0f0d65cdbe18cc5d6bb46c8176eff8c42995")}, .progress= "text")
 


# Madagascar palms ========================
palm_trait <- read.csv(file.path(data.dir,"PalmTraits_10.csv"), na.strings = "")
palm_trait$RedList2012 <- factor(palm_trait$RedList2012, levels = c("LC", "NT", "VU", "EN", "CR", "DD"))
palm_trait$RedList2012_categ <- factor(palm_trait$RedList2012,
                                       levels = c("LC", "NT", "VU", "EN", "CR", "DD"),
                                       labels = c(1,1,2,2,2,0))

madagascar_FS_IUCN <- ggplot(data = subset(palm_trait, !is.na(RedList2012))) +
  geom_point(aes(y = log(AverageFruitLength_cm), x = RedList2012))

madagascar_EOO <- ggplot(aes(y = log10(AverageFruitLength_cm),
                             x = log10(EOO)),
                         data = subset(palm_trait, !is.na(EOO))) +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_text_repel(aes(label = paste0(substr(accGenus,1,1), ". ", accSpecies)),
                  size = 2) +
  labs(x = "Log10 extent of occurrence (km^2)")


madagascar_AOO <- ggplot(aes(y = log10(AverageFruitLength_cm),
                             x = log10(AOO)),
                         data = subset(palm_trait, !is.na(AOO))) +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_text_repel(aes(label = paste0(substr(accGenus,1,1), ". ", accSpecies)),
                  size = 2) +
  labs(x = "Log10 area of occupancy (km^2)")

madagascar_AOO_IUCN <- ggplot(aes(x = RedList2012,
                             y = log10(AOO)),
                         data = subset(palm_trait, !is.na(AOO))) +
  geom_point()

madagascar_EOO_IUCN <- ggplot(aes(x = RedList2012,
                                  y = log10(EOO)),
                              data = subset(palm_trait, !is.na(AOO))) +
  geom_point()
  
ggsave(file.path(fig.dir, "madagascar_IUCN.pdf"), madagascar_FS_IUCN)
ggsave(file.path(fig.dir, "madagascar_AOO.pdf"), madagascar_AOO)
ggsave(file.path(fig.dir, "madagascar_EOO.pdf"), madagascar_EOO)
ggsave(file.path(fig.dir, "madagascar_AOO_IUCN.pdf"), madagascar_AOO_IUCN)
ggsave(file.path(fig.dir, "madagascar_EOO_IUCN.pdf"), madagascar_EOO_IUCN)
