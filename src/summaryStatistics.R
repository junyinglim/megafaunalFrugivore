## GENERATE SUMMARY STATISTICS

## Packages ========================
rm(list = ls())

## Directories ========================
main.dir <- "~/Dropbox/projects/2019/palms/projects/megafaunalFrugivore"
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



## How many extinct species included 
tdwg_final_glob <- read.csv(file.path(res.dir, "tdwgFruitSizeChange.csv"))
countryList <- tdwg_final_glob$LEVEL_3_CO

mammal_presnat_occ_trait <- read.csv(file.path(res.dir, "mammal_presnat_occ_trait.csv"))
mammal_curr_occ_trait <- read.csv(file.path(res.dir, "mammal_curr_occ_trait.csv"))
length(unique(subset(mammal_presnat_occ_trait, IUCN.Status.1.2 %in% c("EW", "EX") & LEVEL_3_CO %in% countryList)$SpecName)) # 23 EX or EW

default_EP <- unique(subset(mammal_presnat_occ_trait, IUCN.Status.1.2 %in% c("EP") & LEVEL_3_CO %in% countryList & Liberal == "Y")$SpecName)
cons_EP <- unique(subset(mammal_presnat_occ_trait, IUCN.Status.1.2 %in% c("EP") & LEVEL_3_CO %in% countryList & Cons == "Y")$SpecName)
supercons_EP <- unique(subset(mammal_presnat_occ_trait, IUCN.Status.1.2 %in% c("EP") & LEVEL_3_CO %in% countryList & SuperCons == "Y")$SpecName) 

length(default_EP) # 177 EP
length(cons_EP) # 119 EP
length(supercons_EP) # 67 EP

length(unique(subset(mammal_curr_occ_trait, LEVEL_3_CO %in% countryList)$SpecName)) # 1604 taxa in the current 
length(unique(subset(mammal_presnat_occ_trait, LEVEL_3_CO %in% countryList)$SpecName)) # 1811 in the present-natural

length(unique(subset(mammal_curr_occ_trait, IUCN.Status.1.2 %in% c("CR", "DD", "EN", "LC", "NT","VU"))$SpecName))
length(unique(subset(mammal_presnat_occ_trait, IUCN.Status.1.2 %in% c("CR", "DD", "EN", "LC", "NT","VU"))$SpecName))

phylacine_trait <- read.csv("~/Dropbox/projects/2019/palms/data/frugivores/PHYLACINE/Trait_data.csv")
curr_status <- table(subset(phylacine_trait, Binomial.1.2 %in% unique(mammal_curr_occ_trait$SpecName))$IUCN.Status.1.2)
sum(curr_status[names(curr_status) %in% c("CR", "DD", "EN", "LC", "NT", "VU")])

pnat_status <- table(subset(phylacine_trait, Binomial.1.2 %in% unique(mammal_presnat_occ_trait$SpecName))$IUCN.Status.1.2)
sum(pnat_status[names(pnat_status) %in% c("CR", "DD", "EN", "LC", "NT", "VU")])

