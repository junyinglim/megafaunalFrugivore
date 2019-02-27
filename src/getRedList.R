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
 


