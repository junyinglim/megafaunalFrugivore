## Analyze Data
# Generate summary statistics at the TDWG-level scale

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

## Import palm dataset ========================
tdwg_final <- read.csv(file.path(data.dir,"tdwg_final.csv"))

## PRELIMINARY PLOTS  ========================
source(file.path(src.dir, "ggmodavg.R"))

# tdwg_shp <- readOGR("~/Dropbox/Projects/2019/palms/data/TDWG/level3/level3.shp")
# tdwg_shp@data$id <- rownames(tdwg_shp@data)
# tdwg_fort <- fortify(tdwg_shp, region = "id")
# tdwg_fort <- merge(tdwg_fort, tdwg_shp@data, by = "id")
# saveRDS(tdwg_fort, file.path(data.dir, "tdwg_shp_fort.rds") )

tdwg_shp <- readRDS(file.path(data.dir, "tdwg_shp_fort.rds"))

tdwg_final$lgm_ens_Tano <- (tdwg_final$lgm_ens_Tmean/10 - 273.15) - tdwg_final$bio1_mean/10
tdwg_final$lgm_ens_Pano <- (tdwg_final$lgm_ens_Pmean/10 ) - tdwg_final$bio12_mean

# # Plot LGM
# z <- merge(tdwg_shp, tdwg_final[c("LEVEL_3_CO", "lgm_ens_Tano","lgm_ens_Pano")], by = "LEVEL_3_CO", all.x = TRUE)
# z2 <- z[order(z$order),]
# tdwg_final$LEVEL_3_CO
# subset(tdwg_final, LEVEL_3_CO == "ZAI")
# lgmTano <- ggplot(data = z2) +
#   geom_polygon(aes(y = lat, x = long, group = group, fill = lgm_ens_Tano )) +
#   scale_fill_viridis(na.value = "grey50") +
#   theme(panel.background =element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
#         axis.title = element_blank())
# ggsave(file.path(fig.dir, "lgmTano.pdf"),lgmTano, width = 10, height = 5)
# 
# lgmPano <- ggplot(data = subset(z2, !LEVEL_3_CO == "ANT")) +
#   geom_polygon(aes(y = lat, x = long, group = group, fill = lgm_ens_Pano )) +
#   scale_fill_viridis(na.value = "grey50") +
#   theme(panel.background =element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
#         axis.title = element_blank())
# ggsave(file.path(fig.dir, "lgmPano.pdf"), lgmPano, width = 10, height = 5)


# tdwg_gg <- melt(tdwg_final, id.vars = c("LEVEL_3_CO", "LEVEL_NAME", "THREEREALM", "medianFruitLengthFilled"), measure.vars = c("presNat_maxBodySize", "curr_maxBodySize")) 
# y <- ggplot(aes(y = log(maxFruitLengthFilled), x = log(value)), data= tdwg_gg) + geom_point(aes(color = THREEREALM)) + facet_wrap(~variable, nrow =2) + geom_text_repel(aes(label = LEVEL_NAME)) + geom_smooth(aes(colour = THREEREALM),method = "lm", se = F)
# ggsave("~/Desktop/maxBodySize.pdf", y, width = 12, height = 12)
# 
# y <- ggplot(aes(y = log(medianFruitLengthFilled), x = deltaMedianBodySize), data= tdwg_final) + geom_point(aes(color = THREEREALM))+ geom_text_repel(aes(label = LEVEL_NAME)) + geom_smooth(aes(colour = THREEREALM),method = "lm", se = F)
# ggsave("~/Desktop/deltaBodySize.pdf", y, width = 12, height = 12)
# # areas with the largest fruit have had the largest losses
# 
# # Proportional loss of 
# y <- ggplot(aes(y = log(medianFruitLengthFilled), x = log(deltaMedianBodySize/presNat_maxBodySize + 1) ), data= tdwg_final) + geom_point(aes(color = THREEREALM))+ geom_text_repel(aes(label = LEVEL_NAME)) + geom_smooth(aes(colour = THREEREALM),method = "lm", se = F)
# ggsave("~/Desktop/propDeltaBodySize.pdf", y, width = 12, height = 12)
# 
# summary(glm(megapalm_nsp ~ THREEREALM*curr_megaHerb_nSp, data= tdwg_final))
# summary(glm(megapalm_nsp ~ THREEREALM*presNat_megaHerb_nSp, data = tdwg_final))
# ggplot(aes(y=megapalm_nsp, x =presNat_megaHerb_nSp-curr_megaHerb_nSp, color = THREEREALM), data = tdwg_final) + geom_point()
# ggplot(aes(y=megapalm_nsp, x =curr_megaHerb_nSp, color = THREEREALM), data = tdwg_final) + geom_point() + geom_smooth(method = "glm")



## DATA SUBSETS ============
tdwg_final_glob <- subset(tdwg_final, REALM_LONG %in% c("Afrotropics", "Neotropics", "IndoMalay", "Australasia"))
tdwg_final_nw <- subset(tdwg_final, REALM_LONG == "Neotropics")
tdwg_final_oww <- subset(tdwg_final, REALM_LONG == "Afrotropics")
tdwg_final_owe <- subset(tdwg_final, REALM_LONG %in% c("Australasia","IndoMalay"))

# Relationships with dependent variable 
z <- melt(tdwg_final_glob, id.var = "medianFruitLengthFilled", measure.vars = c("curr_medianBodySize","globalPC1","globalPC2","globalPC3", "lgm_ens_Tano", "lgm_ens_Pano"))
ggplot(aes(y = log(medianFruitLengthFilled), x = value), data = z) + geom_point() + facet_wrap(~ variable, scale = "free") + geom_smooth(method = "loess")


## MEDIAN BODY SIZE ============
glob_curr_mod <- lm(log(medianFruitLengthFilled) ~ scale(log(curr_medianBodySize)) + 
                      scale(globalPC1) + 
                      scale(globalPC2) + 
                      scale(globalPC3) +
                      scale(soilcount) +
                      scale(lgm_ens_Tano) + 
                      scale(lgm_ens_Pano), 
                    data =tdwg_final_glob, na.action = "na.fail")
glob_pnat_mod <- update(glob_curr_mod, ~.-scale(log(curr_medianBodySize))+
                          scale(log(presNat_medianBodySize)))
summary(glob_curr_mod)
vif(glob_curr_mod)
vif(glob_pnat_mod)
glob_curr_step_mod <- step(glob_curr_mod); summary(glob_curr_step_mod)
glob_pnat_step_mod <- step(glob_pnat_mod); summary(glob_pnat_step_mod)

nw_curr_mod <- lm(log(medianFruitLengthFilled) ~ scale(log(curr_medianBodySize)) + scale(globalPC1) + scale(globalPC2) + scale(globalPC3) + scale(soilcount) + scale(lgm_ens_Tano) + scale(lgm_ens_Pano), data = tdwg_final_nw, na.action = "na.fail")
nw_pnat_mod <- update(nw_curr_mod, ~.-scale(log(curr_medianBodySize)) + log(presNat_medianBodySize))
vif(nw_curr_mod) # check variance inflation factors for collinearities
vif(nw_pnat_mod)
nw_curr_step_mod <- step(nw_curr_mod); summary(nw_curr_step_mod)
nw_pnat_step_mod <- step(nw_pnat_mod); summary(nw_pnat_step_mod)

oww_curr_mod <- lm(log(medianFruitLengthFilled) ~ scale(log(curr_medianBodySize)) + scale(regionalPC1) + scale(regionalPC2) + scale(regionalPC3) + scale(soilcount) + scale(lgm_ens_Tano) + scale(lgm_ens_Pano), data = tdwg_final_oww, na.action = "na.fail")
oww_pnat_mod <- update(oww_curr_mod, ~.-scale(log(curr_medianBodySize)) + scale(log(presNat_medianBodySize)))
vif(oww_curr_mod)
vif(oww_pnat_mod)
oww_curr_step_mod <- step(oww_curr_mod)
oww_pnat_step_mod <- step(oww_pnat_mod)
summary(oww_curr_step_mod)
summary(oww_pnat_step_mod)

owe_curr_mod <- lm(log(medianFruitLengthFilled) ~ scale(log(curr_medianBodySize)) + scale(regionalPC1) + scale(regionalPC2) + scale(regionalPC3) + scale(soilcount) + scale(lgm_ens_Tano) + scale(lgm_ens_Pano), data = tdwg_final_owe, na.action = "na.fail")
owe_pnat_mod <- update(owe_curr_mod, ~.-scale(log(curr_medianBodySize)) + scale(log(presNat_medianBodySize)))
vif(owe_curr_mod)
vif(owe_pnat_mod)
owe_curr_step_mod <- step(owe_curr_mod); summary(owe_curr_step_mod)
owe_pnat_step_mod <- step(owe_pnat_mod); summary(owe_pnat_step_mod)

# Generate output table
medianModSummary <- rbind(
  summarizeLM(glob_curr_step_mod, scale = "Global", scenario = "Current"),
  summarizeLM(glob_pnat_step_mod, scale = "Global", scenario = "Present-natural"),
  summarizeLM(nw_curr_step_mod, scale = "Neotropics", scenario = "Current"),
  summarizeLM(nw_pnat_step_mod, scale = "Neotropics", scenario = "Present-natural"),
  summarizeLM(oww_curr_step_mod, scale = "Afrotropics", scenario = "Current"),
  summarizeLM(oww_pnat_step_mod, scale = "Afrotropics", scenario = "Present-natural"),
  summarizeLM(owe_curr_step_mod, scale = "Indotropics", scenario = "Current"),
  summarizeLM(owe_pnat_step_mod, scale = "Indotropics", scenario = "Present-natural")
)
write.csv(medianModSummary, file.path(res.dir, "medianModSummary.csv"), fileEncoding="UTF-8",
          row.names = F)

## MAXIMUM BODY SIZE ============
glob_curr_mod <- lm(log(maxFruitLengthFilled) ~ scale(log(curr_maxBodySize)) + scale(globalPC1) + scale(globalPC2) + scale(globalPC3) + scale(soilcount) + scale(lgm_ens_Tano) + scale(lgm_ens_Pano), data =tdwg_final_glob, na.action = "na.fail")
glob_pnat_mod <- update(glob_curr_mod, ~.- scale(log(curr_maxBodySize))+ scale(log(presNat_maxBodySize)))
vif(glob_curr_mod)
vif(glob_pnat_mod)
glob_curr_step_mod <- step(glob_curr_mod); summary(glob_curr_step_mod)
glob_pnat_step_mod <- step(glob_pnat_mod); summary(glob_pnat_step_mod)

nw_curr_mod <- lm(log(maxFruitLengthFilled) ~ scale(log(curr_maxBodySize)) + scale(regionalPC1) + scale(regionalPC2) + scale(regionalPC3) + scale(soilcount) + scale(lgm_ens_Tano) + scale(lgm_ens_Pano), data = tdwg_final_nw, na.action = "na.fail")
nw_pnat_mod <- update(nw_curr_mod, ~.- scale(log(curr_maxBodySize))+ scale(log(presNat_maxBodySize)))
vif(nw_curr_mod) # check variance inflation factors for collinearities
vif(nw_pnat_mod)
nw_curr_step_mod <- step(nw_curr_mod); summary(nw_curr_step_mod) 
nw_pnat_step_mod <- step(nw_pnat_mod); summary(nw_pnat_step_mod)


oww_curr_mod <- lm(log(maxFruitLengthFilled) ~ scale(log(curr_maxBodySize)) + scale(regionalPC1) + scale(regionalPC2) + scale(regionalPC3) + scale(soilcount) + scale(lgm_ens_Tano) + scale(lgm_ens_Pano), data = tdwg_final_oww, na.action = "na.fail")
oww_pnat_mod <- update(oww_curr_mod, ~.-scale(log(curr_maxBodySize)) + scale(log(presNat_maxBodySize)))
vif(oww_curr_mod)
vif(oww_pnat_mod)
oww_curr_step_mod <- step(oww_curr_mod); summary(oww_curr_step_mod)
oww_pnat_step_mod <- step(oww_pnat_mod); summary(oww_pnat_step_mod)


owe_curr_mod <- lm(log(maxFruitLengthFilled) ~ scale(log(curr_maxBodySize)) + scale(regionalPC1) + scale(regionalPC2) + scale(regionalPC3) + scale(soilcount) + scale(lgm_ens_Tano) + scale(lgm_ens_Pano), data = tdwg_final_owe, na.action = "na.fail")
owe_pnat_mod <- update(owe_curr_mod, ~.-scale(log(curr_maxBodySize)) + scale(log(presNat_maxBodySize)))
vif(owe_curr_mod)
vif(owe_pnat_mod)
owe_curr_step_mod<- step(owe_curr_mod); summary(owe_curr_step_mod)
owe_pnat_step_mod<- step(owe_pnat_mod); summary(owe_pnat_step_mod)

maxModSummary <- rbind(
  summarizeLM(glob_curr_step_mod, scale = "Global", scenario = "Current"),
  summarizeLM(glob_pnat_step_mod, scale = "Global", scenario = "Present-natural"),
  summarizeLM(nw_curr_step_mod, scale = "Neotropics", scenario = "Current"),
  summarizeLM(nw_pnat_step_mod, scale = "Neotropics", scenario = "Present-natural"),
  summarizeLM(oww_curr_step_mod, scale = "Afrotropics", scenario = "Current"),
  summarizeLM(oww_pnat_step_mod, scale = "Afrotropics", scenario = "Present-natural"),
  summarizeLM(owe_curr_step_mod, scale = "Indotropics", scenario = "Current"),
  summarizeLM(owe_pnat_step_mod, scale = "Indotropics", scenario = "Present-natural")
)
write.csv(maxModSummary, file.path(res.dir, "maxModSummary.csv"), fileEncoding="UTF-8",
          row.names = F)


## STANDARD DEVIATION IN BODY SIZE ============
tdwg_final_glob_range <- subset(tdwg_final_glob, palm_nSp > 1 & presNat_nSp > 1 & curr_nSp > 1)

ggplot(aes(y=fruit_realmrange_z, x= curr_realmrange_z, color = THREEREALM), data = tdwg_final_glob_range) + geom_point() + geom_text_repel(aes(label = LEVEL_NAME)) + geom_smooth(method = "lm")
ggplot(aes(y=fruit_realmrange_z, x= pnat_realmrange_z, color = THREEREALM), data = tdwg_final_glob_range) + geom_point() + geom_smooth(method = "lm")

glob_curr_mod <- lm(fruit_globalrange_z ~ scale(curr_globalrange_z) + scale(globalPC1) + scale(globalPC2) + scale(globalPC3) + scale(soilcount) + scale(lgm_ens_Tano) + scale(lgm_ens_Pano), data = tdwg_final_glob_range )
glob_pnat_mod <- update(glob_curr_mod, ~.- scale(curr_globalrange_z)+ scale(pnat_globalrange_z))
vif(glob_curr_mod)
vif(glob_pnat_mod)
glob_curr_step_mod <- step(glob_curr_mod); summary(glob_curr_step_mod)
glob_pnat_step_mod <- step(glob_pnat_mod); summary(glob_pnat_step_mod)

nw_curr_mod <- lm(fruit_realmrange_z ~ scale(curr_realmrange_z) + scale(regionalPC1) + scale(regionalPC2) + scale(regionalPC3) + scale(soilcount) + scale(lgm_ens_Tano) + scale(lgm_ens_Pano), data = subset(tdwg_final_glob_range, THREEREALM == "NewWorld"))
nw_pnat_mod <- update(nw_curr_mod, ~.- scale(curr_realmrange_z) + scale(pnat_realmrange_z))
vif(nw_curr_mod) # check variance inflation factors for collinearities
vif(nw_pnat_mod)
nw_curr_step_mod <- step(nw_curr_mod); summary(nw_curr_step_mod) 
nw_pnat_step_mod <- step(nw_pnat_mod); summary(nw_pnat_step_mod)

oww_curr_mod <- lm(fruit_realmrange_z ~ scale(curr_realmrange_z) + scale(regionalPC1) + scale(regionalPC2) + scale(regionalPC3) + scale(soilcount) + scale(lgm_ens_Tano) + scale(lgm_ens_Pano), data = subset(tdwg_final_glob_range, THREEREALM == "OWWest"))
oww_pnat_mod <- update(oww_curr_mod, ~.- scale(curr_realmrange_z) + scale(pnat_realmrange_z))
vif(oww_curr_mod)
vif(oww_pnat_mod)
oww_curr_step_mod <- step(oww_curr_mod); summary(oww_curr_step_mod) 
oww_pnat_step_mod <- step(oww_pnat_mod); summary(oww_pnat_step_mod)

owe_curr_mod <- lm(fruit_realmrange_z ~ scale(curr_realmrange_z) + scale(regionalPC1) + scale(regionalPC2) + scale(regionalPC3) + scale(soilcount) + scale(lgm_ens_Tano) + scale(lgm_ens_Pano), data = subset(tdwg_final_glob_range, THREEREALM == "OWEast"))
owe_pnat_mod <- update(owe_curr_mod, ~.- scale(curr_realmrange_z) + scale(pnat_realmrange_z))
vif(oww_curr_mod)
vif(oww_pnat_mod)
owe_curr_step_mod <- step(owe_curr_mod); summary(owe_curr_step_mod) 
owe_pnat_step_mod <- step(owe_pnat_mod); summary(owe_pnat_step_mod)

rangeModSummary <- rbind(
  summarizeLM(glob_curr_step_mod, scale = "Global", scenario = "Current"),
  summarizeLM(glob_pnat_step_mod, scale = "Global", scenario = "Present-natural"),
  summarizeLM(nw_curr_step_mod, scale = "Neotropics", scenario = "Current"),
  summarizeLM(nw_pnat_step_mod, scale = "Neotropics", scenario = "Present-natural"),
  summarizeLM(oww_curr_step_mod, scale = "Afrotropics", scenario = "Current"),
  summarizeLM(oww_pnat_step_mod, scale = "Afrotropics", scenario = "Present-natural"),
  summarizeLM(owe_curr_step_mod, scale = "Indotropics", scenario = "Current"),
  summarizeLM(owe_pnat_step_mod, scale = "Indotropics", scenario = "Present-natural")
)
write.csv(rangeModSummary, file.path(res.dir, "rangeModSummary.csv"), fileEncoding="UTF-8",
          row.names = F)



# NOTES: Might be worth it to include bootstrapped coefficient values as well since the distributionof independent variable is quite skewed.


## SPECIES DIVERSITY
library(pscl)
head(tdwg_final_glob)
sum(tdwg_final_glob$megapalm_nsp == 0, na.rm = TRUE)

glob_sp_mod <- zeroinfl(megapalm_nsp ~ presNat_megaHerb_nSp, data = tdwg_final_glob, dist = "poisson")


glob_sp_mod2 <- glm(megapalm_nsp ~ presNat_megaHerb_nSp+ scale(globalPC1) + scale(globalPC2) + scale(globalPC3) + scale(lgm_ens_Tano) + scale(lgm_ens_Pano)+ THREEREALM, data = tdwg_final_glob, family = "poisson")
glob_sp_mod3 <- glm(megapalm_nsp ~ curr_megaHerb_nSp+ scale(globalPC1) + scale(globalPC2) + scale(globalPC3) + scale(lgm_ens_Tano) + scale(lgm_ens_Pano) + THREEREALM, data = tdwg_final_glob, family = "poisson")
AIC(glob_sp_mod2)

nw_sp_mod1 <- glm(megapalm_nsp ~ presNat_megaHerb_nSp + scale(regionalPC1) + scale(regionalPC2) + scale(regionalPC3) + scale(lgm_ens_Tano) + scale(lgm_ens_Pano), data = tdwg_final_nw, family = "poisson")
nw_sp_mod2 <- glm(megapalm_nsp ~ curr_megaHerb_nSp+ scale(regionalPC1) + scale(regionalPC2) + scale(regionalPC3) + scale(lgm_ens_Tano) + scale(lgm_ens_Pano), data = tdwg_final_nw, family = "poisson")
AIC(nw_sp_mod1)
AIC(nw_sp_mod2)


nw_prop_mod <- glm(propMegaPalm ~ propMegaMam_presnat+ scale(regionalPC1) + scale(regionalPC2) + scale(regionalPC3) + scale(lgm_ens_Tano) + scale(lgm_ens_Pano), data = tdwg_final_nw, family = "binomial")
nw_prop_mod2 <- glm(propMegaPalm ~ propMegaMam_curr+ scale(regionalPC1) + scale(regionalPC2) + scale(regionalPC3) + scale(lgm_ens_Tano) + scale(lgm_ens_Pano), data = tdwg_final_nw, family = "binomial")
summary(step(nw_prop_mod2))



ggplot(aes(y = megapalm_nsp, x = presNat_megaHerb_nSp), data = tdwg_final_glob) +
  geom_point(aes(color = THREEREALM))
ggplot(aes(y = megapalm_nsp, x = curr_megaHerb_nSp), data = tdwg_final_glob) +
  geom_point(aes(color = THREEREALM))


# AIC for present natural is better than that in curr mega herbivore diversity
glob_sp_mod <- glm(megapalm_nsp ~ curr_megaHerb_nSp + presNat_megaHerb_nSp + THREEREALM, data = tdwg_final_glob, family = "poisson")
vif(glob_sp_mod)
residuals


# Model averaging --------
oww_mod_dr <- dredge(oww_mod)
oww_mod_avg <- model.avg(oww_mod_dr)
oww_mod_ri <- summarizeRelImportance(oww_mod_avg)
oww_plot <- plotRelImportance(oww_mod_ri)+ geom_hline(yintercept = 0,linetype = "dotted" )
ggsave(file.path(fig.dir, "oww_mod_>5.pdf"), oww_plot, height = 3, width = 6)





## ASKFJASF +=======
# Number of megafaunal fruit species and the number of megafaunal species
global_mod <- glm(propMegaPalm ~ scale(log(curr_medianBodySize)) + scale(log(presNat_medianBodySize)) + scale(curr_megaHerb_nSp) + scale(presNat_megaHerb_nSp) + scale(globalPC1) + scale(globalPC2) + scale(globalPC3) + scale(lgm_ens_Tano) + scale(lgm_ens_Pano), data = tdwg_final2, na.action = "na.fail", family = "binomial")
global_mod_dr <- dredge(global_mod)
global_mod_avg <- model.avg(global_mod_dr)
global_mod_ri <- summarizeRelImportance(global_mod_avg)
global_plot <- plotRelImportance(global_mod_ri) + geom_hline(yintercept = 0,linetype = "dotted" )

summary(step(global_mod))
ggplot(aes(y = megapalm_nsp, x = log(AREA_KM2)), data = tdwg_final2) + geom_point() + geom_smooth(method = "glm")
curr_megaHerb_nSp
global_mod <- glm(megapalm_nsp ~ presNat_megaHerb_nSp + curr_megaHerb_nSp + log(AREA_KM2) + log(presNat_medianBodySize) + log(curr_medianBodySize), data = tdwg_final2, na.action = "na.fail")
summary(x)


# Current + Global -------
global_curr_data <- tdwg_final2
global_curr_mod <- lm(log(medianFruitLengthFilled) ~ scale(log(curr_medianBodySize)) + globalPC1 + globalPC2 + globalPC3 + lgm_ens_Tano + lgm_ens_Pano,
                      data = global_curr_data,
                      na.action = "na.fail")
global_curr_mod_dr <- dredge(global_curr_mod)
global_curr_mod_avg <- model.avg(global_curr_mod_dr)
global_curr_mod_avg_sum <- summarizeRelImportance(global_curr_mod_avg)
global_curr_mod_plot <- plotRelImportance(global_curr_mod_avg_sum)
ggsave(file.path(fig.dir, "global_curr_relimpt.pdf"),
       global_curr_mod_plot, width = 8, height = 4)

presid <- resid(global_curr_mod, type = "partial")
z <- data.frame(presid = presid[,1], var = scale(log(tdwg_final2$curr_medianBodySize))) 
global_curr_presid_plot <- ggplot(aes(y = presid, x =var),data = z)+ geom_point(shape = 1) + geom_smooth(method = "lm") + labs(y = "Log median fruit size\n(Partial residuals)", x = "Scaled log median body size (Current)") + theme(panel.background = element_blank())
ggsave(file.path(fig.dir, "global_curr_presid.pdf"), global_curr_presid_plot, width = 8, height = 4)

# Present natural global-------
target.col <- c("medianFruitLengthFilled", "presNat_medianBodySize", "NPP_mean", "globalPC1", "globalPC2", "globalPC3", "lgm_ens_Tano", "lgm_ens_Pano", "THREEREALM")
global_pnat_data <- tdwg_final2
global_pnat_mod <- lm(log(medianFruitLengthFilled) ~ scale(log(presNat_medianBodySize)) + globalPC1 + globalPC2 + globalPC3 + lgm_ens_Tano + lgm_ens_Pano,
                      data = global_pnat_data,
                      na.action = "na.fail")
global_pnat_mod_dr <- dredge(global_pnat_mod)
global_pnat_mod_avg <- model.avg(global_pnat_mod_dr)
global_pnat_mod_avg_sum <- summarizeRelImportance(global_pnat_mod_avg)
global_pnat_mod_plot <- plotRelImportance(global_pnat_mod_avg_sum)
ggsave(file.path(fig.dir, "global_pnat_relimpt.pdf"),
       global_pnat_mod_plot, width = 8, height = 4)

presid <- resid(global_pnat_mod, type = "partial")
z <- data.frame(presid = presid[,1], var = scale(log(global_pnat_data$presNat_medianBodySize)))
global_pnat_presid_plot <- ggplot(aes(y = presid, x =var),data = z)+ geom_point(shape = 1) + geom_smooth(method = "lm") + labs(y = "Mean fruit size\n(Partial residuals)", x = "Log median body size (Current)") + theme(panel.background = element_blank())
ggsave(file.path(fig.dir, "global_pnat_presid.pdf"), 
       global_pnat_presid_plot, width = 8, height = 4)


# Current + New World -------
tdwg_final_nw <- subset(tdwg_final2, THREEREALM == "NewWorld")
nw_curr_data <- tdwg_final_nw
nw_curr_mod <- lm(log(meanFruitLengthFilled) ~ scale(log(curr_medianBodySize)) + regionalPC1 + regionalPC2 + regionalPC3 + lgm_ens_Tano + lgm_ens_Pano,
                      data = nw_curr_data,
                      na.action = "na.fail")
nw_curr_mod_dr <- dredge(nw_curr_mod)
nw_curr_mod_avg <- model.avg(nw_curr_mod_dr)
nw_curr_mod_avg_sum <- summarizeRelImportance(nw_curr_mod_avg)
nw_curr_mod_plot <- plotRelImportance(nw_curr_mod_avg_sum)
ggsave(file.path(fig.dir, "nw_curr_relimpt.pdf"),
       nw_curr_mod_plot, width = 8, height = 4)

presid <- resid(nw_curr_mod, type = "partial")
z <- data.frame(presid = presid[,1], var = scale(log(nw_curr_data$curr_medianBodySize))) 
nw_curr_presid_plot <- ggplot(aes(y = presid, x =var),data = z)+ geom_point(shape = 1) + geom_smooth(method = "lm") + labs(y = "Mean fruit size\n(Partial residuals)", x = "Log median body size (Current)") + theme(panel.background = element_blank())
ggsave(file.path(fig.dir, "nw_curr_presid.pdf"), 
       nw_curr_presid_plot, width = 8, height = 4)

# Present-natural + New World -------
nw_pnat_data <- tdwg_final_nw
nw_pnat_mod <- lm(log(meanFruitLengthFilled) ~ scale(log(presNat_medianBodySize)) + regionalPC1 + regionalPC2 + regionalPC3 + lgm_ens_Tano + lgm_ens_Pano,
                  data = nw_pnat_data,
                  na.action = "na.fail")
nw_pnat_mod_dr <- dredge(nw_pnat_mod)
nw_pnat_mod_avg <- model.avg(nw_pnat_mod_dr)
nw_pnat_mod_avg_sum <- summarizeRelImportance(nw_pnat_mod_avg)
nw_pnat_mod_plot <- plotRelImportance(nw_pnat_mod_avg_sum)
ggsave(file.path(fig.dir, "nw_pnat_relimpt.pdf"),
       nw_pnat_mod_plot, width = 8, height = 4)

presid <- resid(nw_pnat_mod, type = "partial")
z <- data.frame(presid = presid[,1], var = scale(log(nw_pnat_data$presNat_medianBodySize))) 
nw_pnat_presid_plot <- ggplot(aes(y = presid, x =var),data = z)+ geom_point(shape = 1) + geom_smooth(method = "lm") + labs(y = "Mean fruit size\n(Partial residuals)", x = "Log median body size (Current)") + theme(panel.background = element_blank())
ggsave(file.path(fig.dir, "nw_pnat_presid.pdf"), 
       nw_pnat_presid_plot, width = 8, height = 4)

# Current Old World East  -------
tdwg_final_owe <- subset(tdwg_final2, THREEREALM == "OWEast")
owe_curr_data <- tdwg_final_owe
owe_curr_mod <- lm(log(meanFruitLengthFilled) ~ scale(log(curr_medianBodySize)) + NPP_mean + regionalPC1 + regionalPC2 + regionalPC3 + lgm_ens_Tano + lgm_ens_Pano,
                  data = owe_curr_data,
                  na.action = "na.fail")
owe_curr_mod_dr <- dredge(owe_curr_mod)
owe_curr_mod_avg <- model.avg(owe_curr_mod_dr)
owe_curr_mod_avg_sum <- summarizeRelImportance(owe_curr_mod_avg)
owe_curr_mod_plot <- plotRelImportance(owe_curr_mod_avg_sum)
ggsave(file.path(fig.dir, "owe_curr_relimpt.pdf"),
       owe_curr_mod_plot, width = 8, height = 4)

presid <- resid(owe_curr_mod, type = "partial")
z <- data.frame(presid = presid[,1], var = scale(log(owe_curr_data$curr_medianBodySize))) 
owe_curr_presid_plot <- ggplot(aes(y = presid, x =var),data = z)+ geom_point(shape = 1) + geom_smooth(method = "lm") + labs(y = "Mean fruit size\n(Partial residuals)", x = "Log median body size (Current)") + theme(panel.background = element_blank())
ggsave(file.path(fig.dir, "owe_curr_presid.pdf"), 
       owe_curr_presid_plot, width = 8, height = 4)

# Present Natural Old World East  -------
owe_pnat_data <- tdwg_final_owe
owe_pnat_mod <- lm(log(meanFruitLengthFilled) ~ scale(log(presNat_medianBodySize)) + NPP_mean + regionalPC1 + regionalPC2 + regionalPC3 + lgm_ens_Tano + lgm_ens_Pano,
                   data = owe_pnat_data,
                   na.action = "na.fail")
owe_pnat_mod_dr <- dredge(owe_pnat_mod)
owe_pnat_mod_avg <- model.avg(owe_pnat_mod_dr)
owe_pnat_mod_avg_sum <- summarizeRelImportance(owe_pnat_mod_avg)
owe_pnat_mod_plot <- plotRelImportance(owe_pnat_mod_avg_sum)
ggsave(file.path(fig.dir, "owe_pnat_relimpt.pdf"),
       owe_pnat_mod_plot, width = 8, height = 4)

presid <- resid(owe_pnat_mod, type = "partial")
z <- data.frame(presid = presid[,1], var = scale(log(owe_pnat_data$presNat_medianBodySize))) 
owe_pnat_presid_plot <- ggplot(aes(y = presid, x =var),data = z)+ geom_point(shape = 1) + geom_smooth(method = "lm") + labs(y = "Mean fruit size\n(Partial residuals)", x = "Log median body size (Current)") + theme(panel.background = element_blank())
ggsave(file.path(fig.dir, "owe_pnat_presid.pdf"), 
       owe_pnat_presid_plot, width = 8, height = 4)

# Current Old World West  -------
target.col <- c("meanFruitLengthFilled", "curr_medianBodySize", "NPP_mean", "PC1", "PC2", "PC3", "lgm_ens_Tano", "lgm_ens_Pano", "THREEREALM")
tdwg_final_oww <- subset(tdwg_final2, THREEREALM == "OWWest")

oww_curr_data <- tdwg_final_oww
oww_curr_mod <- lm(log(meanFruitLengthFilled) ~ scale(log(curr_medianBodySize)) + scale(regionalPC1) + scale(regionalPC2) + scale(regionalPC3) + scale(lgm_ens_Tano) + scale(lgm_ens_Pano),
                   data = oww_curr_data,
                   na.action = "na.fail")
oww_curr_mod_dr <- dredge(oww_curr_mod)
oww_curr_mod_avg <- model.avg(oww_curr_mod_dr)
oww_curr_mod_avg_sum <- summarizeRelImportance(oww_curr_mod_avg)
oww_curr_mod_plot <- plotRelImportance(oww_curr_mod_avg_sum)
ggsave(file.path(fig.dir, "oww_curr_relimpt.pdf"),
       oww_curr_mod_plot, width = 8, height = 4)

presid <- resid(oww_curr_mod, type = "partial")
z <- data.frame(presid = presid[,1], var = scale(log(oww_curr_data$curr_medianBodySize))) 
oww_curr_presid_plot <- ggplot(aes(y = presid, x =var),data = z)+ geom_point(shape = 1) + geom_smooth(method = "lm") + labs(y = "Mean fruit size\n(Partial residuals)", x = "Log median body size (Current)") + theme(panel.background = element_blank())
ggsave(file.path(fig.dir, "oww_curr_presid.pdf"), 
       oww_curr_presid_plot, width = 8, height = 4)

# Present Natural Old World West -------
target.col <- c("meanFruitLengthFilled", "presNat_medianBodySize", "NPP_mean", "PC1", "PC2", "PC3", "lgm_ens_Tano", "lgm_ens_Pano", "THREEREALM")
tdwg_final_oww <- subset(tdwg_final2, CONTINENT == "AFRICA") #THREEREALM == "OWWest")
oww_pnat_data <- tdwg_final_oww

oww_pnat_mod <- lm(log(medianFruitLengthFilled) ~ scale(log(presNat_medianBodySize)) + scale(NPP_mean) + scale(regionalPC1) + scale(regionalPC2) + scale(regionalPC3) + scale(lgm_ens_Tano) + scale(lgm_ens_Pano),
                   data = oww_pnat_data,
                   na.action = "na.fail")
oww_pnat_mod_dr <- dredge(oww_pnat_mod)
oww_pnat_mod_avg <- model.avg(oww_pnat_mod_dr)
oww_pnat_mod_avg_sum <- summarizeRelImportance(oww_pnat_mod_avg)
oww_pnat_mod_plot <- plotRelImportance(oww_pnat_mod_avg_sum)
ggsave(file.path(fig.dir, "oww_pnat_relimpt.pdf"),
       oww_pnat_mod_plot, width = 8, height = 4)

presid <- resid(oww_pnat_mod, type = "partial")
z <- data.frame(presid = presid[,1],
                var = scale(log(oww_pnat_data$presNat_medianBodySize)),
                label = oww_pnat_data$LEVEL_NAME) 
oww_pnat_presid_plot <-
  ggplot(aes(y = presid, x =var), data = z) +
  geom_point(shape = 1) +
  geom_text_repel(aes(label = label)) +
  geom_smooth(method = "loess") + labs(y = "Median fruit size\n(Partial residuals)", x = "Scaled log median body size\n(standard deviations)") +
  theme(panel.background = element_blank())
ggsave(file.path(fig.dir, "oww_pnat_presid.pdf"), 
       oww_pnat_presid_plot, width = 12, height = 6)

library(spdep)
fitted.sarlm()
fitted.

# SEM ===========
library(lavaan)
library(semPlot)
tdwg_final$FruitSize <- log10(tdwg_final$medianFruitLengthFilled)
tdwg_final$CurrBodySize <- log10(tdwg_final$curr_medianBodySize)
tdwg_final$PresNatBodySize <- log10(tdwg_final$presNat_medianBodySize)
tdwg_final$NPP <- tdwg_final$NPP_mean / 10000
tdwg_final$LGM_Tano <- tdwg_final$lgm_ens_Tano / 10
tdwg_final$LGM_Pano <- tdwg_final$lgm_ens_Pano / 100
tdwg_sem_nw <- subset(tdwg_final, THREEREALM == "NewWorld")

varlist <- c("FruitSize", "CurrBodySize", "PresNatBodySize", "soilcount", "NPP", "LGM_Tano", 
             "LGM_Pano", "regionalPC1", "regionalPC2", "regionalPC3")
semMod1 <- 'FruitSize ~ CurrBodySize + soilcount + NPP + LGM_Tano + LGM_Pano + regionalPC1 + regionalPC2 + regionalPC3
CurrBodySize ~ NPP + regionalPC1 + regionalPC2 + regionalPC3 + LGM_Tano + LGM_Pano
NPP ~ soilcount + regionalPC1 + regionalPC2 + regionalPC3 + LGM_Pano + LGM_Tano
'
sem.fit.1 <- sem(semMod1, data=tdwg_sem_nw[varlist])
fitMeasures(sem.fit.1)
summary(sem.fit.1, stand=T, rsq=T, fit.measures=T)

semMod2 <- 'FruitSize ~ PresNatBodySize + soilcount + NPP + LGM_Tano + LGM_Pano + regionalPC1 + regionalPC2 + regionalPC3
PresNatBodySize ~ NPP + regionalPC1 + regionalPC2 + regionalPC3 + LGM_Tano + LGM_Pano
NPP ~ soilcount + regionalPC1 + regionalPC2 + regionalPC3 + LGM_Pano + LGM_Tano'
sem.fit.2 <- sem(semMod2, data = tdwg_sem_nw[varlist])
summary(sem.fit.2, stand=T, rsq=T, fit.measures=T)
semPaths(sem.fit.2, what = "est", residuals = F, intercepts = F, exoCov = FALSE, rotation =3, panelGroups = TRUE, title = T, fade = F)

lavTestLRT(sem.fit.1, sem.fit.2)
modindices(sem.fit.1)
modindices(sem.fit.2)

# Criteria for fit measures:
# Model Chi-Square with its df and p-value: prefer p-value greater than 0.05
# Root Mean Square Error of Approximation (RMSEA): prefer lower 90%CI to be < 0.05
# Comparative Fit Index (CFI): prefer value greater than 0.90
# Mod. indexes should be << 3.84; if higher then may suggest a missing path

# JYL notes:
# Not sure what to look out for in resid() or modindices
#modindices(sem.fit.1)
#resid(sem.fit.1, type="standardized")
