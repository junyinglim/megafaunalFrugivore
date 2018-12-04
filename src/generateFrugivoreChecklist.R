# generateFrugivoreChecklist.R
# Using present-natural ranges from PHYLACINE (Faurby et al 2018) database and IUCN range map, generate frugivore checklists for both current and present-natural conditions

## Packages ========================
library(raster); library(rgdal); library(rgeos)
library(plyr); library(sf)
options(stringsAsFactors = FALSE)

## Directories ========================
main.dir <- "/Users/junyinglim/Dropbox/Projects/2019/palms/projects/megafaunalFrugivore"
projdata.dir <- file.path(main.dir, "data")

rawdata.dir <- "~/Dropbox/Projects/2019/palms/data"
frug.dir <- file.path(rawdata.dir, "frugivores")

## Import diet data ========================
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
phylacine_presentNat_facHerb_SpList <- subset(phylacine_trait, Diet.Plant < 100 & Diet.Plant >= 50 & IUCN.Status.1.2 == "EP")$Binomial.1.2
length(phylacine_presentNat_oblgHerb_SpList) # 178 extinct obligate herbivores
length(phylacine_presentNat_facHerb_SpList) # 24 extinct occassional herbivores

## Create TDWG level checklists ======================== 
# Import tdwg polygons
tdwg.dir <- file.path(rawdata.dir, "TDWG", "level3")
tdwg_shape <- readOGR(file.path(tdwg.dir, "level3.shp"))
tdwg_shape@data$id <- rownames(tdwg_shape@data)
tdwg_shape_reproj <- spTransform(tdwg_shape, crs("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

# Define function for finding overlap of mammal ranges with tdwg boundaries
overlapTDWG <- function(x, tdwg, dir = NULL, type = "raster", polygon = NULL){
  # Finds overlap between raster files and tdwg polygon IDs
  # Args:
  #   x = filename
  #   dir = directory with raster files
  # Returns:
  #   list of dataframes with each species and all tdwg countries it occurs, in long format
  if(type == "raster"){
    sp_raster <- raster(file.path(dir,x))
    sp_polys <- rasterToPolygons(sp_raster, fun = function(x){x>0})
    overlapList <- over(sp_polys, tdwg, returnList = TRUE)
    names(overlapList) <- 1:length(overlapList)
    overlapDF <- ldply(overlapList, .id = "gridcell") # note the grid cells are unique within species but not across species
    overlapDF$SpecName <- gsub(x, pattern = "\\.tif", replacement = "")  
  } else {
    sp_polys <- as(subset(polygon, binomial == x), "Spatial")
    overlapList <- over(sp_polys, tdwg, returnList = TRUE)
    names(overlapList) <- 1:length(overlapList)
    overlapDF <- ldply(overlapList, .id = "polygon")
    overlapDF$SpecName <- gsub(x, pattern = " ", replacement = "_")  
  }
  return(overlapDF)
}

# Identify species with no range modifications
phylacine_metadata <- read.csv(file.path(phylacine.dir, "Spatial_metadata.csv"))
phylacine_synonym <- read.csv(file.path(phylacine.dir, "Synonymy_table_valid_species_only.csv"))
phylacine_synonym$IUCN.Binomial <- paste(phylacine_synonym$IUCN.2016.3.Genus, phylacine_synonym$IUCN.2016.3.Species, sep = " ")

zerorange_splist <- subset(phylacine_metadata, Number.Cells.Current.Range == 0 & Number.Cells.Present.Natural.Range == 0)$Binomial.1.2
nochange_splist <- subset(phylacine_metadata, Change.In.Cells == 0 & (!Binomial.1.2 %in% zerorange_splist))$Binomial.1.2
nochange_iucn_splist <- resolveTaxonomy(nochange_splist, ref = phylacine_synonym, syn.col = "Binomial.1.2", acc.col = "IUCN.Binomial")

# Import original IUCN range map polygons
iucn_range_faurby_polys <- read_sf(file.path(frug.dir, "IUCN\ 2016-3", "MAMMALS.shp"))
iucn_range_polys <- read_sf("~/Desktop/TERRESTRIAL_MAMMALS/TERRESTRIAL_MAMMALS.shp")

sum(nochange_iucn_splist %in% iucn_range_polys$binomial)
length(nochange_iucn_splist)

## Create present-natural ranges  ======================== 
presNat_dir <- file.path(phylacine.dir, "Present_natural")
presNatRanges_fname <- list.files(presNat_dir)
excl_fname <- paste0(c(zerorange_splist, nochange_splist), ".tif")
presNatRanges_fname <- presNatRanges_fname[!presNatRanges_fname %in% excl_fname] # exclude files with zero range

presNat_fname <- paste0(c(mammaldiet_allFrug_SpList, phylacine_presentNat_oblgHerb_SpList), ".tif")
presNat_fname_subset <- presNat_fname[presNat_fname %in% presNatRanges_fname]

presNat_overlap <- ldply(.data = presNat_fname_subset, .fun = overlapTDWG,
                      dir = presNat_dir, tdwg = tdwg_shape_reproj, 
                      .progress = "text")
saveRDS(presNat_overlap, file.path(projdata.dir, "mammal_presnat_occ_raw.rds"))

## Create present-natural for facultative herbivores extinct taxa  ======================== 
presNat_fname2 <- paste0(phylacine_presentNat_facHerb_SpList, ".tif")
presNat_fname2_subset <- presNat_fname2[presNat_fname2 %in% presNatRanges_fname]
presNat_facherb_overlap <- ldply(.data = presNat_fname2_subset, .fun = overlapTDWG,
                                             dir = presNat_dir, tdwg = tdwg_shape_reproj,
                                             .progress = "text")
saveRDS(presNat_facherb_overlap, file.path(projdata.dir, "mammal_presnat_facherb_occ_raw.rds"))

## Create "no change" ranges  ======================== 
nochange_frug_splist <- mammaldiet_allFrug_SpList[mammaldiet_allFrug_SpList %in% nochange_splist] 

#setdiff(gsub(x=nochange_frug_splist, "_", " "), iucn_range_polys$binomial) # 9 species not represented in IUCN shape file but are in Faurby's 2016 version
res <- list()
for(i in 1:length(nochange_frug_splist)){
  print(i)
  print(nochange_frug_splist[i])
  x = gsub("_", " ", nochange_frug_splist[i])
  if(x %in% c("Canis aureus", "Graomys centralis", "Nectomys magdalenae", "Tolypeutes matacus", "Tolypeutes tricinctus", "Mimon koepckeae", "Pteropus pelewensis", "Rattus arfakienis", "Petinomys sagitta")){
    res[[i]] <- overlapTDWG(x, type = "polygon", tdwg = tdwg_shape, polygon = iucn_range_faurby_polys)
  } else if(x %in% c("Peromyscus caniceps", "Peromyscus dickeyi", "Peromyscus pseudocrinitus", "Peromyscus sejugis") ){
    res[[i]] <- data.frame("polygon" = 1, "LEVEL_3_CO" = "MXN", "ID" = 209, "LEVEL_NAME" = "Mexico Northwest", "REGION_NAM" = "Mexico", "CONTINENT" = "NORTHERN AMERICA", "id" = 209, "SpecName" = x)
  } else if(x %in% c("Octodon pacificus") ){
    res[[i]] <- data.frame("polygon" = 1, "LEVEL_3_CO" = "CLC", "ID" = 67, "LEVEL_NAME" = "Chile Central", "REGION_NAM" = "Southern South America", "CONTINENT" = "SOUTHERN AMERICA", "id" = 67, "SpecName" = x)
  } else if(x %in% c("Pteropus howensis") ){
    res[[i]] <- data.frame("polygon" = 1, "LEVEL_3_CO" = "SOL", "ID" = 299, "LEVEL_NAME" = "Solomon Is.", "REGION_NAM" = "Papuasia", "CONTINENT" = "ASIA-TROPICAL", "id" = 299, "SpecName" = x)
  } else {
    res[[i]] <- overlapTDWG(x, type = "polygon", tdwg = tdwg_shape, polygon = iucn_range_polys)
  }
}
nochange_overlap <- do.call("rbind",res)
nochange_overlap$SpecName <- gsub(x = nochange_overlap$SpecName, " ", "_")
saveRDS(nochange_overlap, file.path(projdata.dir, "mammal_nochange_occ_raw.rds"))

# Summary statistics
length(rangeNoData_spList) # 15 with empty ranges for both present natural 
length(mammaldiet_allFrug_SpList) # 1950 all frugivores
length(phylacine_presentNat_oblgHerb_SpList) # 178 extinct herbivores
sum(paste0(mammaldiet_allFrug_SpList, ".tif") %in% presNatRanges_fname) # 1820 with data; 1894 if synonyms resolved
sum(paste0(phylacine_presentNat_oblgHerb_SpList, ".tif") %in% presNatRanges_fname) # 178 with data
length(presNat_fname_subset) # should be the sum of the previous two lines

## Create current ranges  ======================== 
curr.dir <- file.path(phylacine.dir, "Current")
currRanges_fname <- list.files(curr.dir)
curr_fname <- paste0(mammaldiet_allFrug_SpList, ".tif")
rangeNoData_spList <- spatial_metadata$Binomial.1.2[spatial_metadata$Number.Cells.Current.Range == 0]
excl_fname <- paste0(c(zerorange_splist, nochange_splist, rangeNoData_spList), ".tif")

currRanges_fname <- currRanges_fname[!currRanges_fname %in% excl_fname]
curr_fname_subset <- curr_fname[curr_fname %in% currRanges_fname]

curr_overlap <- ldply(.data = curr_fname_subset, .fun = overlapTDWG,
                      dir = curr.dir, type = "raster",
                      tdwg = tdwg_shape_reproj, .progress = "text")
saveRDS(curr_overlap, file.path(projdata.dir, "mammal_current_occ_raw.rds"))
