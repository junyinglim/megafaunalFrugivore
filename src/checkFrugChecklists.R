# Check frugivore checklist
library(plyr)
options(stringsAsFactors = FALSE, row.names = FALSE)

main.dir <- "/Users/junyinglim/Dropbox/Projects/2019/palms/projects/megafaunalFrugivore"
projdata.dir <- file.path(main.dir, "data")
tdwg_data <- read.csv(file.path(projdata.dir, "TDWG_Environment_AllData_2014Dec.csv"))

phylacine_metadata <- read.csv("~/Dropbox/Projects/2019/palms/data/frugivores/PHYLACINE/Spatial_metadata.csv")

presNat_overlap <- readRDS(file.path(projdata.dir, "mammal_presnat_occ_raw.rds"))

nochange_splist <- phylacine_metadata$Binomial.1.2[phylacine_metadata$Change.In.Cells == 0] # taxa that have no range modifications made

presNat_overlap_subset <- subset(presNat_overlap, !SpecName %in% nochange_splist)

# Flagging aberrant occurrences due to the low resolution of raster maps
presNat_overlap_res <- merge(presNat_overlap_subset, tdwg_data[c("LEVEL_3_CO", "ISISLAND", "CONTINENT", "GeologicalOrigin")])

presNat_overlap_res$ISISLAND[presNat_overlap_res$LEVEL_3_CO %in% c("SUM", "JAW", "BOR", "SRL", "NWG")] <- 0 # sumatra, jawa, borneo, sri lanka, and new guinea (SIC = sicilia)
presNat_overlap_res$ISISLAND[presNat_overlap_res$GeologicalOrigin == "continental"] <- 0


# for each species, for each grid cell, are there 1) multiple tdwg countries assigned, and 2) is at least one from a continent, and at least another one a mainland. If true, then return that grid cell.

presNat_overlap_res$CONTINENT_COARSE <- as.vector(presNat_overlap_res$CONTINENT)
presNat_overlap_res$CONTINENT_COARSE[presNat_overlap_res$CONTINENT_COARSE == "AFRICA"] <- "AFRICA"
presNat_overlap_res$CONTINENT_COARSE[presNat_overlap_res$CONTINENT_COARSE %in% c("NORTHERN AMERICA", "SOUTHERN AMERICA", "ASIA-TROPICAL")] <- "AMERICA"
presNat_overlap_res$CONTINENT_COARSE[presNat_overlap_res$CONTINENT_COARSE %in% c("EUROPE", "ASIA-TEMPERATE", "ASIA-TROPICAL")] <- "EURASIA"
presNat_overlap_res$CONTINENT_COARSE[presNat_overlap_res$CONTINENT_COARSE %in% c("AUSTRALASIA")] <- "AUSTRALASIA"

# For this purpose there are only four "continents"  (Eurasia, Africa, Australia and America) 

table(presNat_overlap_res$CONTINENT_COARSE)

checkGridCell <- function(x){
  x$island_flag <- ifelse(sum(c(0,1) %in% x$ISISLAND) == 2, 1, 0)
  x$diffcont_flag <- ifelse(length(unique(x$CONTINENT_COARSE)) == 2, 1, 0)
  return(x)
}
test <- ddply(.data = presNat_overlap_res,
              .variables = .(SpecName),
              .fun = function(x){
                flaglist <- ddply(.data = x,
                                  .variables = .(gridcell),
                                  .fun = checkGridCell)
                return(flaglist)},
              .progress = "text")
test2 <- subset(test, island_flag == 1 | diffcont_flag == 1)
dim(test)
dim(test2)
write.csv(test2, "~/Desktop/flaggedPresentNaturalRanges.csv", row.names = FALSE)
subset(test, SpecName == "Loxodonta_africana" & LEVEL_3_CO == "CNY")

library(rgdal)
mammalMaps <- readOGR("/Users/junyinglim/Desktop/IUCN\ 2016-3/MAMMALS.shp")
mammalMaps
