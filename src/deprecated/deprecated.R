## DEPRECATED CODE 


## EXPLORE MAMMALDIET V2 DATASET ================
library(readxl)
mammalDiet2 <- read_excel(file.path(frug.dir, "Combined MammalDIET and MammalDIET 2 dataset.xlsx"))
mammalDiet2 <- as.data.frame(mammalDiet2)
subset(mammalDiet2)

sum(mammalDiet2$binomial %in% mammalDiet$SpecName)

nrow(mammalDiet2)
nrow(mammalDiet)
head(mammalDiet2)
7216+4927
nrow(iucnMaps)

## IUCN TAXONOMIC REFERENCE DATA ================
# Import IUCN taxonomic information
iucnTaxonRef <- read.csv(file.path(frug.dir, "iucnMammalTaxonList.csv"))
iucnTaxonRef$SpecName <- paste(iucnTaxonRef$Genus, iucnTaxonRef$Species)
sum(iucnSpList %in% iucnTaxonRef$SpecName)
length(iucnSpList) # 10 species that are not in the taxon dataset, not even at the species level
iucnSpList[!iucnSpList %in% iucnTaxonRef$SpecName] # 2 species with spatial data but not taxonomic data

# Create synonymy table
iucnSynTab <- readRDS(file.path(main.dir, "iucnSynTab.rds"))

# Determine IUCN range map overlap with mammaldiet dataset
mammalDiet$accName <- iucnSynTab$SpecName[match(mammalDiet$SpecName, iucnSynTab$name)]
mammalDiet$accName[is.na(mammalDiet$accName)] <- mammalDiet$SpecName[is.na(mammalDiet$accName)]

sum(mammalDiet$SpecName %in% iucnSpList) # before updating names
sum(mammalDiet$accName %in% iucnSpList) # after updating names

iucnSpList[! iucnSpList %in% iucnSynTab$SpecName] # 2 speces where maps exist but not in the taxonomic reference

# 
faurbySpList2 <- faurbySpList
faurbySpList2[which(!is.na(match(faurbySpList, iucnSynTab$name)))] <- as.vector(iucnSynTab$SpecName[na.omit(match(faurbySpList, iucnSynTab$name))])

sum(faurbySpList %in% iucnSpList)
sum(faurbySpList2 %in% iucnSpList) # 649 matched to IUCN range maps

#
sum(faurbySpList %in% iucnSynTab$SpecName)
sum(faurbySpList2 %in% iucnSynTab$SpecName) 

length(iucnSpList)
sum(iucnSpList %in% iucnSynTab$SpecName)
faurbySpList2[!faurbySpList2 %in% iucnSynTab$SpecName]

## UNDERSTANDING THE FAURBY NUMBERS




## ESTIMATE THE TDWG COVERAGE FOR FAURBY DATASET ================
faurbyFilesSubset <- faurbyFiles[which(faurbyList %in% frugivoreList)]
faurbySpList <- gsub(faurbyFiles, pattern = "\\.tif", replacement = "")
library(raster)
library(rgdal)

# Import tdwg boundaries
tdwg_shape <- readOGR(file.path(main.dir, "TDWG", "level3", "level3.shp"))

tdwg_shape_reproj <- spTransform(tdwg_shape, crs("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

overlapTDWG <- function(x, tdwg){
  temp <- raster(file.path(faurby.dir,x))
  tempPoly <- rasterToPolygons(temp, fun = function(x){x>0})
  overlap <- over(tempPoly, tdwg)
  LEVEL_3_CO <- unique(overlap$LEVEL_3_CO) # even a small overlap is counted here
  return(data.frame("SpecName" = gsub(x, pattern = ".tif", replacement = ""), "LEVEL_3_CO" = LEVEL_3_CO))
}

faurbyOccList <- lapply(faurbyFilesSubset, FUN = overlapTDWG, tdwg = tdwg_shape_reproj)
faurbyOcc <- do.call("rbind",faurbyOccList)
write.csv(faurbyOcc, file.path(frug.dir, "faubryOcc.csv"))

# Cross-reference Faurby list with the TDWG mammal occ
library(plyr)
test <- subset(mammalOccSubset, SpecName %in% faurbyList)
mammalOccCountryCount <- ddply(.data = test, .variables = .(SpecName), .fun = function(x){ data.frame(ntaxa = nrow(x)) })
faurbyOccCountryCount <- ddply(.data = faurbyOcc, .variables = .(SpecName), .fun = function(x){ data.frame(ntaxa = nrow(x)) })

test2 <- merge(x = mammalOccCountryCount, y = faurbyOccCountryCount, by = "SpecName")

plot(ntaxa.y~ntaxa.x, data = test2)
abline(0,1)

# Which frugivore species will have SMALLER present-natural ranges than current ranges
subset(test2, ntaxa.x > ntaxa.y)


## IUCN RANGE MAPS ========================
#iucnMaps <- readOGR(file.path("~/Desktop/TERRESTRIAL_MAMMALS/TERRESTRIAL_MAMMALS.shp"))
#saveRDS(iucnMaps, "~/Desktop/iucnMammals.rds")
iucnMaps <- readRDS(file.path(frug.dir, "iucnMammals.rds"))
#write.csv(iucnMaps@data, "~/Desktop/iucnMetadata.csv")

# Collapse subspecies by species
iucnMaps@data$binomial <- as.vector(iucnMaps@data$binomial)
subspInd <- grep(iucnMaps@data$binomial, pattern = "ssp")
iucnMaps@data$binomial[subspInd] <- str_split(iucnMaps@data$binomial[subspInd], n = 2, pattern = " ssp\\.", simplify = TRUE)[,1]

length(iucnMaps@data$binomial) # 12655 polygons
sum(duplicated(iucnMaps@data$binomial)) # 7224 duplicates (i.e., ranges in different polygons or subspecies ranges)
length(unique(iucnMaps@data$binomial)) # 5431 species
iucnSpList <- as.vector(unique(iucnMaps@data$binomial))

sum(mammalDiet$SpecName %in% iucnSpList) # 4927 species in the dataset
mammalDiet$SpecName[!mammalDiet$SpecName %in% iucnSpList]

## MERGE DATASETS  ========================
iucnSpList_acc <- resolveTaxonomy(iucnSpList, ref = iucnTaxonRef, syn.col = "name", acc.col = "accName")
sum(iucnSpList_acc  %in% iucnTaxonRef$accName) # length(iucnSpList)
length(unique(iucnTaxonRef$accName))
iucnSpList[is.na(iucnSpList_acc)]


sum(allFrugivoreSpList %in% iucnSpList) / length(allFrugivoreSpList)


sum((allFrugivoreSpList) %in% iucnSpList)


# Clean up taxonomy
test <- resolveTaxonomy(allFrugivoreDietName, ref = iucnTaxonRef, syn.col = "name", acc.col = "accName")
sum(is.na(test)) # 559 with no name in the iucn
allFrugivoreDietName[is.na(test)]



# Import Faurby tif dataset
faurby.dir <- file.path(frug.dir, "FaurbyMaster", "SPECIES_FILES")
faurbyFiles <- list.files(faurby.dir, pattern = "\\.tif") # 1000 taxa 
faurbyList <- gsub(faurbyFiles, pattern = ".tif", replacement = "")

sum(frugivoreList %in% faurbyList) # 229 current species in Faurby list

sum(faurbyList %in% mammalDiet$SpecName) # 711 in the faurby list in mammalDiet
sum(faurbyList %in% frugivoreList) # 229 extant frugivores in Faurby dataset

setdiff(faurbyList, frugivoreList) # 772 in fauby not in frugivore list
faurbyList == frugivoreList
setdiff(faurbyList, mammalDiet$SpecName)


faurbyMeta <- read.csv(file.path(frug.dir, "faurbyMetadataJYL.csv"), stringsAsFactors = FALSE)
faurbyMeta$SpecName <- gsub(faurbyMeta$Bininomial.name, pattern = "_", replacement = " ")
test <- resolveTaxonomy(x = faurbyMeta$SpecName, ref = iucnSynTab, syn.col = "name", acc.col ="accName")
test[is.na(test)] <- faurbyMeta$SpecName[is.na(test)]

test2 <- resolveTaxonomy(x = faurbyList, ref = iucnSynTab, syn.col = "name", acc.col ="accName")
test2[is.na(test2)] <- faurbyList[is.na(test2)]

table(faurbyMeta$Taxonomy)


extinctTaxonNames <- subset(faurbyMeta, Taxonomy == "Added")$SpecName
sum(extinctTaxonNames %in% faurbyList)
extinctTaxonNames[!extinctTaxonNames %in% faurbyList]

