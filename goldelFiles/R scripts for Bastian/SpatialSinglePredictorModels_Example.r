##############################################################################
##
##  Palms: Spatial single predictor models for islands and mainlands (full dataset)
##
##  W.D. Kissling 
##
##  Updates:
##  August 2010: including LGM climates + REALM as predictors
##  Sept. 2010: including ensemble LGM climates and number of soil types,
##              + spatial models
##
##
##############################################################################



#----------------------------- PRELIMINARIES ----------------------------------#
#########################
#Workspace
rm(list=ls())
setwd("C://Workspace/R")

#########################
#Libraries
library(foreign)
library(spdep)
library(ncf)
library(rgdal)
library(sp)
library(maptools)


#########################
#Load data
palms<-read.table("TDWG_Environment_AllData_Reduced.txt", header=T, sep="\t")       #366 units, Antarctica and Marcus Island not included
str(palms)

#------------------------------- DATA -----------------------------------------#

#Full dataset
palms_full<-subset(palms, palms$PALMSR>0)                                     #194 units with palms
str(palms_full)
palms_full$REALM_LONG<-factor(palms_full$REALM_LONG)
str(palms_full)


#----------------------------- COORDINATES ------------------------------------#

######################
#Coordinates full dataset
coords_full<-cbind(palms_full$LONG,palms_full$LAT)
coords_full<-as.matrix(coords_full)
str(coords_full)
print(coords_full)

#-------------------------- NEIGHBOURHOODS ------------------------------------#

######################
#1 neighbour
nbk1_full<-knn2nb(knearneigh(coords_full, k=1, longlat = T))
summary(nbk1_full)
plot(nbk1_full, coords_full, pch=20)
nbk1_full_sym <- make.sym.nb(nbk1_full)    #make symmetric
sw1_full <- nb2listw(nbk1_full_sym)
plot(sw1_full,coords_full, pch=20)

######################
#Spatial weights with gabriel connection
gabr_nb <- graph2nb(gabrielneigh(coords_full), sym=TRUE)
summary(gabr_nb)
plot(gabr_nb, coords_full, pch=20)
gabr_full <- nb2listw(gabr_nb)
plot(gabr_full,coords_full, pch=20)
summary(gabr_full)

######################
#weighted by inverse distance
nbk_all<-dnearneigh(coords_full, d1=0, d2=10000, longlat = T)
#plot(nbk_all, coords, pch=20)
nbk_all_sym <- make.sym.nb(nbk_all)    #make symmetric
sw_all <- nb2listw(nbk_all_sym)
#plot(sw_all,coords, pch=20)
dsts <- nbdists(nbk_all, coords_full)
idw <- lapply(dsts, function(x) 1/(x/1000))
sw_ivd <- nb2listw(nbk_all, glist = idw, style = "B")

######################
#identify max distance
dsts1 <- unlist(nbdists(nbk1_full, coords_full, longlat = T))     #get 
summary(dsts1)
max(dsts1)

######################
#all sites with at least 1 neighbour (3500 km)
nbk_3500<-dnearneigh(coords_full, d1=0, d2=3500, longlat = T)
summary(nbk_3500, coords_full, pch=20)
nbk_3500_sym <- make.sym.nb(nbk_3500)    #make symmetric
sw_3500 <- nb2listw(nbk_3500_sym)
summary(sw_3500,coords_full, pch=20)



#------------------------ SPATIAL MODELS -------------------------------------#

#Data
str(palms_full)

######################
#Precipitation
glm_prec_full<-lm(log(palms_full$PALMSR)~log(palms_full$PREC_Sum)) 
summary(glm_prec_full)
sem_prec<-errorsarlm(glm_prec_full, listw=gabr_full, na.action=na.omit)
slm_prec_fit <- print.sarlm.pred(predict.sarlm(sem_prec))
cor(x=slm_prec_fit$fit, y=log(palms_full$PALMSR), method="pearson")^2     #Full model fit
cor(x=slm_prec_fit$trend, y=log(palms_full$PALMSR), method="pearson")^2   #Environment only
summary(sem_prec)
moran.mc(residuals(sem_prec), listw=sw_3500, alternative = "greater", nsim=1000)
#moran.mc(residuals(sem_prec), listw=sw1_full, alternative = "greater", nsim=1000)
#moran.mc(residuals(sem_prec), listw=sw_ivd, alternative = "greater", nsim=1000)

plot(log(palms_full$PALMSR)~log(palms_full$PREC_Sum), xlab="log(precipitation)", ylab="log(richness)", cex=1.5, cex.axis=1.5, cex.lab=1.5, pch=20)
abline(glm_prec_full, col="red", lwd =4)

######################
#Seasonality in precipitation
glm_prec_cv_full<-lm(log(palms_full$PALMSR)~log(palms_full$PREC_CV)) 
sem_prec_cv<-errorsarlm(glm_prec_cv_full, listw=gabr_full, na.action=na.omit)
slm_prec_cv_fit <- print.sarlm.pred(predict.sarlm(sem_prec_cv))
cor(x=slm_prec_cv_fit$fit, y=log(palms_full$PALMSR), method="pearson")^2     #Full model fit
cor(x=slm_prec_cv_fit$trend, y=log(palms_full$PALMSR), method="pearson")^2   #Environment only
summary(sem_prec_cv)
moran.mc(residuals(sem_prec_cv), listw=sw_3500, alternative = "greater", nsim=1000)


######################
#Precipitation of driest quarter
glm_prec_dry_quart_full<-lm(log(palms_full$PALMSR)~log(palms_full$P_drie_quart)) 
sem_prec_dry_quart<-errorsarlm(glm_prec_dry_quart_full, listw=gabr_full, na.action=na.omit)
slm_prec_dry_quart_fit <- print.sarlm.pred(predict.sarlm(sem_prec_dry_quart))
cor(x=slm_prec_dry_quart_fit$fit, y=log(palms_full$PALMSR), method="pearson")^2     #Full model fit
cor(x=slm_prec_dry_quart_fit$trend, y=log(palms_full$PALMSR), method="pearson")^2   #Environment only
summary(sem_prec_dry_quart)
moran.mc(residuals(sem_prec_dry_quart), listw=sw_3500, alternative = "greater", nsim=1000)


######################
#Temperature
glm_temp_full<-lm(log(palms_full$PALMSR)~log(palms_full$Tmean_mean/10+273.15)) 
plot(log(palms_full$PALMSR)~log(palms_full$Tmean_mean/10+273.15))
sem_temp<-errorsarlm(glm_temp_full, listw=gabr_full, na.action=na.omit)
sem_temp_fit <- print.sarlm.pred(predict.sarlm(sem_temp))
cor(x=sem_temp_fit$fit, y=log(palms_full$PALMSR), method="pearson")^2     #Full model fit
cor(x=sem_temp_fit$trend, y=log(palms_full$PALMSR), method="pearson")^2   #Environment only
summary(sem_temp)
moran.mc(residuals(sem_temp), listw=sw_3500, alternative = "greater", nsim=1000)


######################
#Seasonality in temperature
glm_temp_sd_full<-lm(log(palms_full$PALMSR)~log(palms_full$Temp_SD/10)) 
sem_temp_sd<-errorsarlm(glm_temp_sd_full, listw=gabr_full, na.action=na.omit)
sem_temp_sd_fit <- print.sarlm.pred(predict.sarlm(sem_temp_sd))
cor(x=sem_temp_sd_fit$fit, y=log(palms_full$PALMSR), method="pearson")^2     #Full model fit
cor(x=sem_temp_sd_fit$trend, y=log(palms_full$PALMSR), method="pearson")^2   #Environment only
summary(sem_temp_sd)
moran.mc(residuals(sem_temp_sd), listw=sw_3500, alternative = "greater", nsim=1000)


######################
#Temperature coldest quarter
glm_temp_cold_full<-lm(log(palms_full$PALMSR)~log(palms_full$T_cold_quart/10+273.15)) 
plot(log(palms_full$PALMSR)~log(palms_full$Tmin_cold_month/10+273.15))
sem_temp_cold<-errorsarlm(glm_temp_cold_full, listw=gabr_full, na.action=na.omit)
sem_temp_cold_fit <- print.sarlm.pred(predict.sarlm(sem_temp_cold))
cor(x=sem_temp_cold_fit$fit, y=log(palms_full$PALMSR), method="pearson")^2     #Full model fit
cor(x=sem_temp_cold_fit$trend, y=log(palms_full$PALMSR), method="pearson")^2   #Environment only
summary(sem_temp_cold)
moran.mc(residuals(sem_temp_cold), listw=sw_3500, alternative = "greater", nsim=1000)


######################
#Soil types
glm_soil<-lm(log(palms_full$PALMSR)~log(palms_full$soilcount)) 
sem_soil<-errorsarlm(glm_soil, listw=gabr_full, na.action=na.omit)
sem_soil_fit <- print.sarlm.pred(predict.sarlm(sem_soil))
cor(x=sem_soil_fit$fit, y=log(palms_full$PALMSR), method="pearson")^2     #Full model fit
cor(x=sem_soil_fit$trend, y=log(palms_full$PALMSR), method="pearson")^2   #Environment only
summary(sem_soil)
moran.mc(residuals(sem_soil), listw=sw_3500, alternative = "greater", nsim=1000)


######################
#Altitudinal range
glm_alt_range_full<-lm(log(palms_full$PALMSR)~log(palms_full$alt_range)) 
sem_alt_range<-errorsarlm(glm_alt_range_full, listw=gabr_full, na.action=na.omit)
sem_alt_range_fit <- print.sarlm.pred(predict.sarlm(sem_alt_range))
cor(x=sem_alt_range_fit$fit, y=log(palms_full$PALMSR), method="pearson")^2     #Full model fit
cor(x=sem_alt_range_fit$trend, y=log(palms_full$PALMSR), method="pearson")^2   #Environment only
summary(sem_alt_range)
moran.mc(residuals(sem_alt_range), listw=sw_3500, alternative = "greater", nsim=1000)

   
######################
#AREA
glm_area_full<-lm(log(palms_full$PALMSR)~log(palms_full$AREA_KM2)) 
sem_area<-errorsarlm(glm_area_full, listw=gabr_full, na.action=na.omit)
sem_area_fit <- print.sarlm.pred(predict.sarlm(sem_area))
cor(x=sem_area_fit$fit, y=log(palms_full$PALMSR), method="pearson")^2     #Full model fit
cor(x=sem_area_fit$trend, y=log(palms_full$PALMSR), method="pearson")^2   #Environment only
summary(sem_area)
moran.mc(residuals(sem_area), listw=sw_3500, alternative = "greater", nsim=1000)


######################
#LGM precipitation (ensemble)
glm_LGM_prec_ens<-lm(log(palms_full$PALMSR)~log(palms_full$ensLGM_Pmean)) 
sem_LGM_prec<-errorsarlm(glm_LGM_prec_ens, listw=gabr_full, na.action=na.omit)
sem_LGM_prec_fit <- print.sarlm.pred(predict.sarlm(sem_LGM_prec))
cor(x=sem_LGM_prec_fit$fit, y=log(palms_full$PALMSR), method="pearson")^2     #Full model fit
cor(x=sem_LGM_prec_fit$trend, y=log(palms_full$PALMSR), method="pearson")^2   #Environment only
summary(sem_LGM_prec)
moran.mc(residuals(sem_LGM_prec), listw=sw_3500, alternative = "greater", nsim=1000)


######################
#LGM temperature (ensemble)
glm_LGM_temp_ens<-lm(log(palms_full$PALMSR)~log(palms_full$ensLGM_Tmean/10+273.15))
sem_LGM_temp<-errorsarlm(glm_LGM_temp_ens, listw=gabr_full, na.action=na.omit)
sem_LGM_temp_fit <- print.sarlm.pred(predict.sarlm(sem_LGM_temp))
cor(x=sem_LGM_temp_fit$fit, y=log(palms_full$PALMSR), method="pearson")^2     #Full model fit
cor(x=sem_LGM_temp_fit$trend, y=log(palms_full$PALMSR), method="pearson")^2   #Environment only
summary(sem_LGM_temp)
moran.mc(residuals(sem_LGM_temp), listw=sw_3500, alternative = "greater", nsim=1000)


######################
#PREC anomaly (ensemble)
glm_ensAnoP_full<-lm(log(palms_full$PALMSR)~log(abs(palms_full$ensLGM_Pano)))
sem_LGM_prec_ano<-errorsarlm(glm_ensAnoP_full, listw=gabr_full, na.action=na.omit)
sem_LGM_prec_ano_fit <- print.sarlm.pred(predict.sarlm(sem_LGM_prec_ano))
cor(x=sem_LGM_prec_ano_fit$fit, y=log(palms_full$PALMSR), method="pearson")^2     #Full model fit
cor(x=sem_LGM_prec_ano_fit$trend, y=log(palms_full$PALMSR), method="pearson")^2   #Environment only
summary(sem_LGM_prec_ano)
moran.mc(residuals(sem_LGM_prec_ano), listw=sw_3500, alternative = "greater", nsim=1000)


######################
#TEMP anomaly (ensemble)
glm_ensAnoT_full<-lm(log(palms_full$PALMSR)~log(palms_full$ensLGM_Tano)) 
sem_LGM_temp_ano<-errorsarlm(glm_ensAnoT_full, listw=gabr_full, na.action=na.omit)
sem_LGM_temp_ano_fit <- print.sarlm.pred(predict.sarlm(sem_LGM_temp_ano))
cor(x=sem_LGM_temp_ano_fit$fit, y=log(palms_full$PALMSR), method="pearson")^2     #Full model fit
cor(x=sem_LGM_temp_ano_fit$trend, y=log(palms_full$PALMSR), method="pearson")^2   #Environment only
summary(sem_LGM_temp_ano)
moran.mc(residuals(sem_LGM_temp_ano), listw=sw_3500, alternative = "greater", nsim=1000)



######################
#Realm: New World vs. Old World East vs. Old World West
#glm_realm_full<-lm(log(palms_full$PALMSR)~as.factor(palms_full$THREEREALM))
#sem_realm<-errorsarlm(glm_realm_full, listw=gabr_full, na.action=na.omit)
#sem_realm_fit <- print.sarlm.pred(predict.sarlm(sem_realm))
#cor(x=sem_realm_fit$fit, y=log(palms_full$PALMSR), method="pearson")^2     #Full model fit
#cor(x=sem_realm_fit$trend, y=log(palms_full$PALMSR), method="pearson")^2   #Environment only
#summary(sem_realm)
#moran.mc(residuals(sem_realm), listw=sw1_full, alternative = "greater", nsim=1000)


######################
#Realm_Long
glm_realm_long_full<-lm(log(palms_full$PALMSR)~as.factor(palms_full$REALM_LONG))
sem_realm_long<-errorsarlm(glm_realm_long_full, listw=gabr_full, na.action=na.omit)
sem_realm_fit_long <- print.sarlm.pred(predict.sarlm(sem_realm_long))
cor(x=sem_realm_fit_long$fit, y=log(palms_full$PALMSR), method="pearson")^2     #Full model fit
cor(x=sem_realm_fit_long$trend, y=log(palms_full$PALMSR), method="pearson")^2   #Environment only
summary(sem_realm_long)
moran.mc(residuals(sem_realm_long), listw=sw_3500, alternative = "less", nsim=1000)

#boxplot(log(palms_mainland_full$PALMSR)~as.factor(palms_mainland_full$REALM_LONG))


######################
#Is island
glm_isisland_full<-lm(log(palms_full$PALMSR)~as.factor(palms_full$ISISLAND))
sem_island<-errorsarlm(glm_isisland_full, listw=gabr_full, na.action=na.omit)
sem_island_fit <- print.sarlm.pred(predict.sarlm(sem_island))
cor(x=sem_island_fit$fit, y=log(palms_full$PALMSR), method="pearson")^2     #Full model fit
cor(x=sem_island_fit$trend, y=log(palms_full$PALMSR), method="pearson")^2   #Environment only
summary(sem_island)
moran.mc(residuals(sem_island), listw=sw_3500, alternative = "greater", nsim=1000)


