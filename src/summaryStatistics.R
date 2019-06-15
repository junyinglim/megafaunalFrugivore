## GENERATE SUMMARY STATISTICS

## Packages ========================
rm(list = ls())

## Directories ========================
main.dir <- "/Users/junyinglim/Dropbox/Projects/2019/palms/projects/megafaunalFrugivore"
data.dir <- file.path(main.dir, "data")
res.dir <- file.path(main.dir, "results")
src.dir <- file.path(main.dir, "src")
fig.dir <- file.path(main.dir, "figs")
rawdata.dir <- "~/Dropbox/Projects/2019/palms/data"
frug.dir <- file.path(rawdata.dir, "frugivores")
options(stringsAsFactors =FALSE)

## Import data ========================
# Import results
tdwg_final <- read.csv(file.path(res.dir, "tdwgFruitSizeChange.csv"))
curr_occ_trait <- read.csv(file.path(res.dir, "mammal_curr_occ_trait.csv"))
pnat_occ_trait <- read.csv(file.path(res.dir, "mammal_presnat_occ_trait.csv"))
palm_occ_trait <- read.csv(file.path(res.dir, "tdwg_palm_occ_trait.csv"))
dim(tdwg_final)

# Define realms
globList <- tdwg_final$LEVEL_3_CO
afrolist <- subset(tdwg_final, THREEREALM == "OWWest")$LEVEL_3_CO
indolist <- subset(tdwg_final, THREEREALM == "OWEast")$LEVEL_3_CO
neolist <- subset(tdwg_final, THREEREALM == "NewWorld")$LEVEL_3_CO

# Realm subsets
curr_occ_trait_glob <- subset(curr_occ_trait, LEVEL_3_CO %in% globList)
pnat_occ_trait_glob <- subset(pnat_occ_trait, LEVEL_3_CO %in% globList)

## Fruit and body size variation ========================
mean(tdwg_final$medianFruitLengthFilled)
sd(tdwg_final$medianFruitLengthFilled)
range(tdwg_final$medianFruitLengthFilled)

mean(tdwg_final$max95FruitLengthFilled)
sd(tdwg_final$max95FruitLengthFilled)
range(tdwg_final$max95FruitLengthFilled)

mean(tdwg_final$curr_medianBodySize)
sd(tdwg_final$curr_medianBodySize)
range(tdwg_final$curr_medianBodySize)

mean(tdwg_final$presNat_medianBodySize)
sd(tdwg_final$presNat_medianBodySize)
range(tdwg_final$presNat_medianBodySize)

mean(tdwg_final$curr_max95BodySize)
sd(tdwg_final$curr_max95BodySize)
range(tdwg_final$curr_max95BodySize)

mean(tdwg_final$presNat_max95BodySize)
sd(tdwg_final$presNat_max95BodySize)
range(tdwg_final$presNat_max95BodySize)

## Largest frugivores in each realm ========================
afro_curr_trait <- subset(curr_occ_trait, LEVEL_3_CO %in% afrolist)[c("SpecName","Mass.g")]
afro_curr_trait <- afro_curr_trait[!duplicated(afro_curr_trait),]
head(afro_curr_trait[order(afro_curr_trait$Mass.g, decreasing = T),], n = 10)

indo_curr_trait <- subset(curr_occ_trait, LEVEL_3_CO %in% indolist)[c("SpecName","Mass.g")]
indo_curr_trait <- indo_curr_trait[!duplicated(indo_curr_trait),]
head(indo_curr_trait[order(indo_curr_trait$Mass.g, decreasing = T),], n = 10)

neo_curr_trait <- subset(curr_occ_trait, LEVEL_3_CO %in% neolist)[c("SpecName","Mass.g")]
neo_curr_trait <- neo_curr_trait[!duplicated(neo_curr_trait),]
head(neo_curr_trait[order(neo_curr_trait$Mass.g, decreasing = T),], n = 10)

afro_pnat_trait <- subset(pnat_occ_trait, LEVEL_3_CO %in% afrolist)[c("SpecName","Mass.g")]
afro_pnat_trait <- afro_pnat_trait[!duplicated(afro_pnat_trait),]
head(afro_pnat_trait[order(afro_pnat_trait$Mass.g, decreasing = T),], n = 10)

indo_pnat_trait <- subset(pnat_occ_trait, LEVEL_3_CO %in% indolist)[c("SpecName","Mass.g")]
indo_pnat_trait <- indo_pnat_trait[!duplicated(indo_pnat_trait),]
head(indo_pnat_trait[order(indo_pnat_trait$Mass.g, decreasing = T),], n = 10)

neo_pnat_trait <- subset(pnat_occ_trait, LEVEL_3_CO %in% neolist)[c("SpecName","Mass.g")]
neo_pnat_trait <- neo_pnat_trait[!duplicated(neo_pnat_trait),]
head(neo_pnat_trait[order(neo_pnat_trait$Mass.g, decreasing = T),], n = 10)
dim(neo_pnat_trait[neo_pnat_trait$Mass.g > 1000000,])

tapply(tdwg_final$curr_max95BodySize, INDEX = tdwg_final$THREEREALM, FUN = mean)
tapply(tdwg_final$presNat_max95BodySize, INDEX = tdwg_final$THREEREALM, FUN = mean)

## Changds in fruit size  ========================
mean(tdwg_final$changeInMedFruitSize, na.rm = T)
sd(tdwg_final$changeInMedFruitSize, na.rm = T)

mean(tdwg_final$changeInMaxFruitSize, na.rm = T)
sd(tdwg_final$changeInMaxFruitSize, na.rm = T)
tdwg_final$changeInMaxFruitSize



# Mean changes; largest fruits are going to be disproportionately affected so degree of change is likely to be much larger than suggested by the mean values
subset(tdwg_final, THREEREALM == "OWEast" & max95FruitLengthFilled > 10)
palmTrait <- read.csv(file.path(res.dir, "tdwg_palm_occ_trait.csv"))
quantile(subset(palmTrait, Area_code_L3 == "LSI")$AverageFruitLength_cm_filled, probs = 0.95)

length(tdwg_final$medianFruitLengthFilled)
mean(tdwg_final$curr_medianBodySize)
mean(tdwg_final$presNat_medianBodySize)
mean(tdwg_final$curr_max95BodySize)
mean(tdwg_final$presNat_max95BodySize)
summary(glob_pnat_medBS_sar_mod)







## SCRATCH SHEET
levels(fruitsizechange$medFSchange_categ)
subset(fruitsizechange,  changeInMaxFruitSize > 3)[c("curr_max95BodySize", "futr_maxBodySize")]

mammal_curr_occ_trait <- read.csv(file.path(res.dir, "mammal_curr_occ_trait.csv"))
mammal_pnat_occ_trait <- read.csv(file.path(res.dir, "mammal_presnat_occ_trait.csv"))
subset(mammal_curr_occ_trait, LEVEL_3_CO == "CBD")[c("SpecName", "IUCN.Status.1.2", "Mass.g")]

mean(fruitsizechange$changeInMedFruitSize)

names(fruitsizechange)
subset(tdwg_final, CONTINENT == "AUSTRALASIA")$presNat_medianBodySize
subset(tdwg_final, LEVEL_3_CO == "BZN")$presNat_max95BodySize
subset(tdwg_final, LEVEL_3_CO == "BZN")$curr_max95BodySize
subset(tdwg_final, LEVEL_3_CO == "BZN")$futr_maxBodySize
subset(fruitsizechange, LEVEL_3_CO == "BZN")$changeInMaxFruitSize

hist(subset(mammal_curr_occ_trait, LEVEL_3_CO == "BZN")$Mass.g)
hist(log(subset(mammal_curr_occ_trait, LEVEL_3_CO == "BZN")$Mass.g))
quantile(subset(mammal_curr_occ_trait, LEVEL_3_CO == "BZN")$Mass.g, probs = 0.95)
quantile(log(subset(mammal_curr_occ_trait, LEVEL_3_CO == "BZN")$Mass.g), probs = 0.95)

subset(mammal_curr_occ_trait, LEVEL_3_CO == "BZN" & IUCN.Status.1.2 %in% c("VU", "EN", "CR"))$Mass.g

subset(mammal_curr_occ_trait, LEVEL_3_CO == "BZN" & IUCN.Status.1.2 == "VU")
quantile(subset(mammal_curr_occ_trait, LEVEL_3_CO == "BZN")$Mass.g, probs = 0.95)
quantile(subset(mammal_curr_occ_trait, LEVEL_3_CO == "BZN" & (!IUCN.Status.1.2 %in% c(c("VU","CR","EW","EN","EX"))))$Mass.g, probs = 0.95)
table(mammal_curr_occ_trait$IUCN.Status.1.2)

IUCN.Status.1.2 %in% c("CR","EW","EN","EX")




ggplot() +
  geom_polygon(aes(y = lat, x = long, group = group), data = tdwg_shp2) +
  coord_cartesian(ylim = c(-14, 24), xlim = c(80, 160))
