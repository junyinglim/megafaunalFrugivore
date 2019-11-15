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
subset(mammalDietFam_summary, Family == "DASYPODIDAE")

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


## Determining extinction rates ========================
# Di Marco et al (2014) data was obtained by digitizing Table 1
dimarcodata <- matrix( c(262, 24, 15, 6, 3, 0,
                         3, 20, 18, 4, 1, 0,
                         1, 4, 52, 28, 3, 1,
                         1, 0, 2, 27, 10, 0,
                         0, 0, 1, 2, 9, 1,
                         0, 0, 0, 0, 0, 0), byrow = T, ncol = 6)

## Hoffmann et al (2010) data was obtained by cross-referencing Table S5 and S8
# to determine the distribution of Red List category changes across species at the start
# and end of each time point

# hoffmann_data<- read.csv("~/Desktop/Hoffmann_genuinechanges.csv")
# hoffmann_mam_data <- subset(hoffmann_data, Class == "Mammalia")
# hoffmann_mam_data$CategoryStart <- gsub(hoffmann_mam_data$Category.at.start.of.period, pattern = "CR\\(PE\\)|CR\\(PEW\\)|EW", replacement = "EX")
# hoffmann_mam_data$CategoryEnd <- gsub(hoffmann_mam_data$Category.at.end.of.period, pattern = "CR\\(PE\\)|CR\\(PEW\\)|EW", replacement = "EX")
# hoffmann_mam_data$StatusChange <- paste(hoffmann_mam_data$CategoryStart, hoffmann_mam_data$CategoryEnd, sep = "_")
# table(hoffmann_mam_data$StatusChange)
# subset(hoffmann_mam_data, StatusChange %in% c("EX_EN", "EX_CR"))

hoffmanndata <- matrix( c(3112, 39, 2, 2, 0, 0,
                          7, 281, 47, 4, 1, 0,
                          1, 5, 438, 39, 2, 0,
                          0, 1, 3, 400, 31, 0,
                          0, 0, 2, 3, 129, 4,
                          0, 0, 0, 1, 1, 98), byrow = T, ncol = 6)


redlisttrans <- function(LC_NT, NT_LC, NT_VU, 
                         VU_EN, VU_NT, EN_CR,
                         EN_VU, CR_EX, CR_EN,
                         EX_CR, t){
 # Generates transition probabilities given a set of instantaneous rates of change between red list categories
 # Args:
 #     LC_NT = transition rate between LC and NT category
 #     t = amount of time
 # Returns:
 #     matrix, transition probabilities
 LC_LC <- -1 * LC_NT
 NT_NT <- -1 * (NT_LC + NT_VU)
 VU_VU <- -1 * (VU_NT + VU_EN)
 EN_EN <- -1 * (EN_VU + EN_CR)
 CR_CR <- -1 * (CR_EN + CR_EX)
 EX_EX <- -1 * EX_CR
 Q <- matrix( c(LC_LC, LC_NT, 0,     0,     0,     0,
                NT_LC, NT_NT, NT_VU, 0,     0,     0,
                0,     VU_NT, VU_VU, VU_EN,     0,     0,
                0,     0,     EN_VU, EN_EN, EN_CR,     0,
                0,     0,     0,     CR_EN, CR_CR, CR_EX,
                0,     0,     0,     0,     EX_CR, EX_EX), byrow = T, ncol = 6)
  return(expm::expm(Q * t, order = 8))
}
  
redlistLik <- function(data, par, EX_CR = NULL, ...){
  # Calculate the likelihood given a matrix of numbers representing change (or not) in Red List status
  # 
  #data = redlistdata
  #par = rep(1E-6, 9)
  if(is.null(EX_CR)){
    P <- redlisttrans(par[1], par[2], par[3], par[4], par[5], par[6], par[7], par[8], par[9], par[10], ...)  
  } else {
    P <- redlisttrans(par[1], par[2], par[3], par[4], par[5], par[6], par[7], par[8], par[9], EX_CR = EX_CR, ...)  
  }
  
  P_s <- P / rowSums(P) # standardize rows
  #print(log(P_s[data > 0]))
  return(-1 * sum(log(P_s[data > 0]) * data[data>0]))
}

init_par <- rep(0.01, 10)
dimarcoML <- optim(par = init_par[1:9], fn = redlistLik, method = "L-BFGS-B",
           lower = 1E-6, upper = 1, data = dimarcodata, t = 33, EX_CR = 0)

hoffmannML <- optim(par = init_par, fn = redlistLik, method = "L-BFGS-B",
                    lower = 1E-6, upper = 1, data = hoffmanndata, t = 12)
# testdata <- hoffmanndata
# testdata[6,] <- 0
# hoffmannML2 <- optim(par = init_par[1:9], fn = redlistLik, method = "L-BFGS-B",
#                     lower = 1E-6, upper = 1, data = testdata, t = 12, EX_CR = 0)

# Under this model, what is the probability of going extinct after 100 years if you start as a LC
dimarcoML_pExt <- redlisttrans(dimarcoML$par[1], dimarcoML$par[2], dimarcoML$par[3],
                               dimarcoML$par[4], dimarcoML$par[5], dimarcoML$par[6],
                               dimarcoML$par[7], dimarcoML$par[8], dimarcoML$par[9],
                               EX_CR = 0, t= 100)

hoffmannML_pExt <- redlisttrans(hoffmannML$par[1], hoffmannML$par[2], hoffmannML$par[3],
                               hoffmannML$par[4], hoffmannML$par[5], hoffmannML$par[6],
                               hoffmannML$par[7], hoffmannML$par[8], hoffmannML$par[9],
                               hoffmannML$par[10], t= 100)

ctmc_pExt <- data.frame(IUCN.Status = c("LC", "NT", "VU", "EN", "CR"),
                        dimarco = dimarcoML_pExt[1:5,6],
                        hoffmann = hoffmannML_pExt[1:5,6])

write.csv(ctmc_pExt, file.path(data.dir, "ctmc_pExt.csv"),row.names = FALSE)

# Parametric bootstrap
dimarco_Pt <- redlisttrans(dimarcoML$par[1],
                           dimarcoML$par[2],
                           dimarcoML$par[3],
                           dimarcoML$par[4],
                           dimarcoML$par[5],
                    dimarcoML$par[6],
                    dimarcoML$par[7], 
                    dimarcoML$par[8],
                    dimarcoML$par[9], t= 33, EX_CR = 0)

hoffmann_Pt <- redlisttrans(hoffmannML$par[1],
                              hoffmannML$par[2],
                              hoffmannML$par[3],
                              hoffmannML$par[4],
                              hoffmannML$par[5],
                              hoffmannML$par[6],
                              hoffmannML$par[7],
                              hoffmannML$par[8],
                              hoffmannML$par[9],
                              hoffmannML$par[10], t = 12)



redlist_categ <- c("LC", "NT", "VU", "EN", "CR", "EX")
dimarcodata_init <- rowSums(dimarcodata)
hoffmanndata_init <- rowSums(hoffmanndata)

dimarco_sim <- list()
hoffmann_sim <- list()

for(i in 1:1000){
  dimarco_matrix <- matrix(nrow = 6, ncol = 6, dimnames = list(redlist_categ, redlist_categ))
  hoffmann_matrix <- matrix(nrow = 6, ncol = 6, dimnames = list(redlist_categ, redlist_categ))
  for(j in 1:6){
    dimarco_sim_categ <- sample( x = c("LC", "NT", "VU", "EN", "CR", "EX"), size = dimarcodata_init[j],
                         prob = dimarco_Pt[j, ], replace = T)
    hoffmann_sim_categ <- sample( x = c("LC", "NT", "VU", "EN", "CR", "EX"), size = hoffmanndata_init[j],
                                  prob = hoffmann_Pt[j, ], replace = T)
    
    dimarco_matrix[j,] <- table(factor(dimarco_sim_categ, levels = redlist_categ))
    hoffmann_matrix[j,] <- table(factor(hoffmann_sim_categ, levels = redlist_categ))
  }
  dimarco_sim[[i]] <- dimarco_matrix
  hoffmann_sim[[i]] <- hoffmann_matrix
}

dimarco_sim_MLpara <- list()
hoffmann_sim_MLpara <- list()
for(i in 1:1000){
  print(i)
  dimarco_sim_MLpara[[i]] <- 
    optim(par = init_par[1:9], fn = redlistLik, method = "L-BFGS-B",
          lower = 1E-6, upper = 1, data = dimarco_sim[[i]], EX_CR = 0, t = 33)$par
  
  hoffmann_sim_MLpara[[i]] <- 
    optim(par = init_par, fn = redlistLik, method = "L-BFGS-B",
          lower = 1E-6, upper = 1, data = hoffmann_sim[[i]], t = 12)$par
}

par_labels <- c("LC to NT", "NT to LC", "NT to VU", "VU to EN", "VU to NT", "EN to CR", "EN to VU", "CR to EX", "CR to EN", "EX to CR")
dimarco_par_sim <- data.frame( value = unlist(dimarco_sim_MLpara) ) 
dimarco_par_sim$par_label <- rep(par_labels[1:9], 1000)

dimarco_par_obs <- data.frame( value = dimarcoML$par)
dimarco_par_obs$par_label <- par_labels[1:9]

hoffmann_par_sim <- data.frame( value  = unlist(hoffmann_sim_MLpara) )
hoffmann_par_sim$par_label <- rep(par_labels[1:10], 1000)

hoffmann_par_obs <- data.frame(value = hoffmannML$par)
hoffmann_par_obs$par_label <- par_labels

dimarco_parboot <- ggplot() +
  geom_histogram(aes(x = value), bins = 30, data = dimarco_par_sim) +
  geom_vline(aes(xintercept = value), data = dimarco_par_obs) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  labs(y = "No. of simulations", x = "Instantaneous rates", title = "Di Marco et al. 2014") +
  facet_wrap(~par_label) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

hoffmann_parboot <- ggplot() +
  geom_histogram(aes(x = value), bins = 30, data = hoffmann_par_sim) +
  geom_vline(aes(xintercept = value), data = hoffmann_par_obs) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  labs(y = "No. of simulations", x = "Instantaneous rates", title = "Hoffmann et al. 2010") +
  facet_wrap(~par_label) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

library(cowplot)
parbootplot <- plot_grid(hoffmann_parboot, dimarco_parboot, labels = "auto")

ggsave(parbootplot, filename = file.path(fig.dir, "parbootplot.pdf"), width = 12, height = 8)

ctmc_pExt
0.1/0.013891699
0.6723/0.088429901
0.9990/0.204767517
