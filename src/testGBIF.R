library(plyr)
main.dir <- "~/Dropbox/projects/2019/palms/projects/megafaunalFrugivore/"
data.dir <- file.path(main.dir, "data")
fig.dir <- file.path(main.dir, "figs")


## On the adequacy of GBIF data ====================
gbif_occ <- read.delim(file.path(data.dir, "gbif_0002199-191105090559680.csv")) # search with Arecaceae as family

gbif_occ_subset <- subset(gbif_occ, taxonRank == "SPECIES" & decimalLatitude != "" & decimalLongitude != "")

species_Nocc <- ddply(.data = gbif_occ_subset, .variables = .(species),
                      .fun = function(x){ nOcc <- nrow(x); data.frame( nOcc  ) })
species_Nocc$species <- as.vector(species_Nocc$species) 
palm_trait <- read.csv(file.path(data.dir, "PalmTraits_10.csv"))

sum(species_Nocc$species %in% palm_trait$SpecName) # 1687 species are in palm traits
length(species_Nocc$species) # 1849 species

library(Taxonstand)
# For taxa that do not match the trait database, check for synonyms
nomatch <- species_Nocc$species[!species_Nocc$species %in% palm_trait$SpecName]
synonym_tab <- TPL(splist = as.vector(nomatch))

synonym_tab_subset <- subset(synonym_tab, Plant.Name.Index = T)
synonym_tab_subset$AcceptedSp <- paste(synonym_tab_subset$New.Genus, synonym_tab_subset$New.Species, sep = " ")

index <- match(species_Nocc$species, synonym_tab_subset$Taxon)
species_Nocc$species[which(!is.na(index))] <- synonym_tab_subset$AcceptedSp[index[!is.na(index)]]

# Merge the duplicates after synonyms have been corrected
species_Nocc_final <- ddply(.data = species_Nocc, .variables = .(species),
                          .fun = summarize,
                          nOcc = sum(nOcc))

sum(species_Nocc_final$species %in% palm_trait$SpecName) # 1720
length(species_Nocc_final$species)

species_Nocc_final$PalmTraits <- species_Nocc_final$species %in% palm_trait$SpecName
species_Nocc_final$OccGroup <- ifelse(species_Nocc_final$nOcc > 1, "More than 1 record", "Singleton")

# Plot
library(ggplot2)
gbif_stat <- ggplot() +
  geom_bar( aes(x = factor(ifelse(nOcc > 1, "More than 1 record", "Singleton")),
                fill = PalmTraits), data = species_Nocc_final, stat = "count") +
  labs(y = "No of species") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 12),
        legend.position = "bottom")
  
ggsave(gbif_stat, filename = file.path(fig.dir, "gbif_summary.pdf"), width = 4, height = 4)

species_Nocc_final_multi <- subset(species_Nocc_final, OccGroup == "More than 1 record")
species_Nocc_final_multi$species %in% gbif_occ_subset$species

library(raster)
r<- raster(ext = extent(-180, 180, -90, 90), res=c(2.5,2.5))
p <- rasterToPolygons(r)


# Import tdwg units
library(sp)
library(rgdal)
palm_occ <- read.csv(file.path(data.dir, "palms_in_tdwg3.csv"))
palm_occ$SpecName <- gsub(palm_occ$SpecName, pattern = "_", replacement = " ")
palm_occ_subset <- subset(palm_occ, SpecName %in% species_Nocc_final$species )

tdwg_shp_raw<- readOGR("~/Dropbox/Projects/2019/palms/data/TDWG/level3/level3.shp")
tdwg_shp_raw2 <- spTransform(tdwg_shp_raw, CRS(proj4string(p)))

getGrid <- function(x){
  listofcountries <- x$Area_code_L3
  id <- rownames(over(subset(tdwg_shp_raw2, LEVEL_3_CO %in% listofcountries), y= p, returnList = T)[[1]])
  data.frame(id)
}

palm_occ_grid <- ddply(.data = palm_occ_subset, .variables = .(SpecName), .fun = getGrid, .progress = "text" )

getGridOcc 

subset(listofcountries)
test <- subset(tdwg_shp_raw2, LEVEL_3_CO == "THA")
head(tdwg_shp_raw@data)
proj4string(p)

## On frugivore classifications ========================
rawdata.dir <- "~/Dropbox/Projects/2019/palms/data"
frug.dir <- file.path(rawdata.dir, "frugivores")

mammalDiet <- read.delim(file.path(frug.dir, "mammalDiet_v1.0.txt"), stringsAsFactors = FALSE)
mammalDiet$SpecName <- paste(mammalDiet$Genus, mammalDiet$Species, sep = "_")


mammalDietFamilyList <- unique(mammalDiet$Family)

mammalDietFam_summary <- ddply(.data = mammalDiet, .variables = .(Family),
                               summarize,
                               famPropFrug =  sum(Fruit %in% 1:2, na.rm = T)/length(Fruit),
                               famPropOblgFrug = sum(Fruit == 1, na.rm = T)/ length(Fruit))

mammalDietOrder_summary <- ddply(.data = mammalDiet, .variables = .(Order),
                                 summarize,
                                 orderPropFrug =  sum(Fruit %in% 1:2, na.rm = T)/length(Fruit),
                                 orderPropOblgFrug = sum(Fruit == 1, na.rm = T)/ length(Fruit))


# Create species list of obligate and facultative frugivores
mammaldiet_allFrug_SpList <- subset(mammalDiet, Fruit >= 1)$SpecName
mammaldiet_oblgFrug_SpList <- subset(mammalDiet, Fruit == 1)$SpecName
mammaldiet_facFrug_SpList <- subset(mammalDiet, Fruit %in% c(3))$SpecName
#ggplot(data = mammalDiet) + geom_bar(aes(x = Order, fill = factor(Fruit)), stat= "count") +
#  theme(axis.text = element_text(angle = 90))


phylacine.dir <- file.path(frug.dir, "PHYLACINE")
phylacine_trait <- read.csv(file.path(phylacine.dir, "Trait_data.csv"), stringsAsFactors = F)
extinctTaxa <- subset(phylacine_trait, IUCN.Status.1.2 == "EP")
extinctTaxa$Family <- toupper(extinctTaxa$Family.1.2)
extinctTaxa$Order <- toupper(extinctTaxa$Order.1.2)

extinctTaxa_Fam <- merge(x= extinctTaxa, y = mammalDietFam_summary, by = "Family", all.x = TRUE)
extinctTaxa_All <- merge(x = extinctTaxa_Fam, y = mammalDietOrder_summary, by = "Order", all.x = TRUE)

library(ggrepel)
ggplot(data = extinctTaxa_All) +
  geom_point(aes(y = famPropFrug, x = Diet.Plant, color = Order)) +
  scale_colour_brewer(palette = "ra")
  #geom_text_repel(aes(y = famPropFrug, x = Diet.Plant, label = Binomial.1.2), size = 1)
?scale_colour_brewer

#ggplot(data = extinctTaxa_All) + geom_boxplot(aes(y = orderPropFrug, x = factor(ifelse(Diet.Plant >= 50, 1, 0))))

# Which taxa are 
write.csv(extinctTaxa_All, file = "~/Desktop/PleistoceneHerbivoires.csv", row.names = F)

unique(extinctTaxa$Order)[! unique(extinctTaxa$Order) %in% mammalDietOrder_summary$Order]
subset(extinctTaxa, Order %in% c("LITOPTERNA", "NOTOUNGULATA"))
unique(extinctTaxa$Family)[! unique(extinctTaxa$Family) %in% mammalDietFam_summary$Family]
subset(mammalDietFam_summary, Family == "DASYPODIDAE")
subset(mammalDiet, Family =="MEGALONYCHIDAE")


