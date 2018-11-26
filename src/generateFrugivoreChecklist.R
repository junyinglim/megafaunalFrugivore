# generateFrugivoreChecklist.R
## Packages ========================
library(raster)
library(rgdal)
library(rgeos)
library(plyr)

## Directories ========================
main.dir <- "/Users/junyinglim/Dropbox/Projects/2019/palms/projects/megafaunalFrugivore"
projdata.dir <- file.path(main.dir, "data")

rawdata.dir <- "~/Dropbox/Projects/2019/palms/data"
frug.dir <- file.path(rawdata.dir, "frugivores")

## Import diet data ========================
#palmOcc <- read.csv(file.path("~/Dropbox/Projects/2019/palms/data/", "WorldChecklistData", "palms_in_tdwg3.csv"))
#palmAreas <- as.vector(unique(palmOcc$Area_code_L3)) # tdwg units with palms

# IUCN taxonomy
source(file.path(main.dir, "src","resolveTaxonomy.R"))
iucnTaxonRef <- readRDS(file.path(projdata.dir, "iucnSynTab.rds"))

# Import MammalDiet v1.0 dataset
mammalDiet <- read.delim(file.path(frug.dir, "mammalDiet_v1.0.txt"), stringsAsFactors = FALSE)
mammalDiet$SpecName <- paste(mammalDiet$Genus, mammalDiet$Species, sep = "_")

# Create species list of obligate and facultative frugivores
mammaldiet_allFrug_SpList <- subset(mammalDiet, Fruit >= 1)$SpecName
mammaldiet_oblgFrug_SpList <- subset(mammalDiet, Fruit == 1)$SpecName

temp <- resolveTaxonomy(x = gsub("_", " ", mammaldiet_allFrug_SpList),
                        ref=iucnTaxonRef, syn.col = "name", acc.col = "accName")
sum(is.na(temp)) # 26 names in mammaldiet that are not in iucn taxonomy
temp[is.na(temp)] <- mammaldiet_allFrug_SpList[is.na(temp)]
mammaldiet_allFrug_SpList <- gsub(" ", "_", temp)

# Import phylacine dataset
phylacine.dir <- file.path(frug.dir, "PHYLACINE")
phylacine_trait <- read.csv(file.path(phylacine.dir, "Trait_data.csv"), stringsAsFactors = FALSE)
spatial_metadata <- read.csv(file.path(phylacine.dir, "Spatial_metadata.csv"), stringsAsFactors = FALSE)

phylacine_presentNat_oblgHerb_SpList <- subset(phylacine_trait, Diet.Plant == 100 & IUCN.Status.1.2 == "EP")$Binomial.1.2
length(phylacine_presentNat_oblgHerb_SpList) # 178 extinct herbivores

## Create TDWG level checklists 
# Import tdwg polygons
tdwg.dir <- file.path(rawdata.dir, "TDWG", "level3")
tdwg_shape <- readOGR(file.path(tdwg.dir, "level3.shp"))
tdwg_shape@data$id <- rownames(tdwg_shape@data)
tdwg_shape_reproj <- spTransform(tdwg_shape, crs("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

# Define function for finding overlap of mammal ranges with tdwg boundaries
overlapTDWG <- function(x, tdwg, dir){
  # Finds overlap between raster files and tdwg polygon IDs
  # Args:
  #   x = filename
  #   dir = directory with raster files
  # Returns:
  #   list of dataframes with each species and all tdwg countries it occurs, in long format
  sp_raster <- raster(file.path(dir,x))
  sp_polys <- rasterToPolygons(sp_raster, fun = function(x){x>0})
  overlapList <- over(sp_polys, tdwg, returnList = TRUE)
  names(overlapList) <- 1:length(overlapList)
  overlapDF <- ldply(overlapList, .id = "gridcell") # note the grid cells are unique within species but not across species
  overlapDF$SpecName <- gsub(x, pattern = "\\.tif", replacement = "")
  return(overlapDF)
}

# Identify species with no range modifications
nochange_splist <- subset(phylacine_metadata, Change.In.Cells == 0)$Binomial.1.2

# Original IUCN 
nochange_splist

# Create present-natural ranges
presNat.dir <- file.path(phylacine.dir, "Present_natural")
presNatRanges_fname <- list.files(presNat.dir)

rangeNoData_spList <- spatial_metadata$Binomial.1.2[spatial_metadata$Number.Cells.Present.Natural.Range == 0 &  spatial_metadata$Number.Cells.Current.Range == 0] # species with empty range
excl_fname <- paste0(rangeNoData_spList, ".tif")
presNatRanges_fname <- presNatRanges_fname[!presNatRanges_fname %in%  excl_fname] # exclude files with zero range

presNat_fname <- paste0(c(mammaldiet_allFrug_SpList, phylacine_presentNat_oblgHerb_SpList), ".tif")
presNat_fname_subset <- presNat_fname[presNat_fname %in% presNatRanges_fname] 

presNat_overlap <- ldply(.data = presNat_fname_subset,
                      .fun = overlapTDWG,
                      dir = presNat.dir,
                      tdwg = tdwg_shape_reproj, buffer = FALSE,
                      .progress = "text")

saveRDS(presNat_overlap, file.path(projdata.dir, "mammal_presnat_occ_raw.rds"))

write.table(presNat_list, file = file.path(projdata.dir, "mammal_presNat_occ.txt"), row.names = FALSE, sep = "\t")

# Summary statistics
length(rangeNoData_spList) # 15 with empty ranges for both present natural 
length(mammaldiet_allFrug_SpList) # 1950 all frugivores
length(phylacine_presentNat_oblgHerb_SpList) # 178 extinct herbivores
sum(paste0(mammaldiet_allFrug_SpList, ".tif") %in% presNatRanges_fname) # 1820 with data; 1894 if synonyms resolved
sum(paste0(phylacine_presentNat_oblgHerb_SpList, ".tif") %in% presNatRanges_fname) # 178 with data
length(presNat_fname_subset) # should be the sum of the previous two lines




# Create current ranges
curr.dir <- file.path(phylacine.dir, "Current")
currRanges_fname <- list.files(curr.dir)
rangeNoData_spList <- spatial_metadata$Binomial.1.2[spatial_metadata$Number.Cells.Current.Range == 0] # species with empty range
curr_fname <- paste0(mammaldiet_allFrug_SpList, ".tif")
excl_fname <- paste0(rangeNoData_spList, ".tif")
currRanges_fname <- currRanges_fname[!currRanges_fname %in% excl_fname]

curr_fname_subset <- curr_fname[curr_fname %in% currRanges_fname] # only include those with a non-empty range 
length(rangeNoData_spList)
length(curr_fname)
sum(curr_fname %in% currRanges_fname)
length(curr_fname_subset)
curr_list <- ldply(.data = curr_fname_subset,
                   .fun = overlapTDWG,
                   dir = curr.dir,
                   tdwg = tdwg_shape_reproj, buffer = FALSE,
                   .progress = "text")
write.table(curr_list, file = file.path(projdata.dir, "mammal_current_occ.txt"), row.names = FALSE, sep = "\t")

# Identify which grid cells overlap with both continental and 
main.dir <- "/Users/junyinglim/Dropbox/Projects/2019/palms/projects/megafaunalFrugivore"
projdata.dir <- file.path(main.dir, "data")
tdwg_data <- read.csv(file.path(projdata.dir, "TDWG_Environment_AllData_2014Dec.csv"))

presNat_overlap <- readRDS(file.path(projdata.dir, "mammal_presnat_occ_raw.rds"))

# Flagging aberrant occurrences due to the low resolution of raster maps

presNat_overlap_res <- merge(presNat_overlap, tdwg_data[c("LEVEL_3_CO", "ISISLAND", "CONTINENT")])
presNat_overlap_res$ISISLAND[presNat_overlap_res$LEVEL_3_CO %in% c("SUM", "JAW", "BOR", "SRL", "NWG")] <- 0 # sumatra, jawa, borneo, sri lanka, and new guinea (SIC = sicilia)

# for each species, for each grid cell, are there 1) multiple tdwg countries assigned, and 2) is at least one from a continent, and at least another one a mainland. If true, then return that grid cell.

presNat_overlap_res$CONTINENT_COARSE <- as.vector(presNat_overlap_res$CONTINENT)
presNat_overlap_res$CONTINENT_COARSE[presNat_overlap_res$CONTINENT_COARSE == "AFRICA"] <- "AFRICA"
presNat_overlap_res$CONTINENT_COARSE[presNat_overlap_res$CONTINENT_COARSE %in% c("EUROPE", "ASIA-TEMPERATE", "ASIA-TROPICAL")] <- "EURASIA"
presNat_overlap_res$CONTINENT_COARSE[presNat_overlap_res$CONTINENT_COARSE %in% c("NORTHERN AMERICA", "SOUTHERN AMERICA", "ASIA-TROPICAL")] <- "AMERICA"
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
