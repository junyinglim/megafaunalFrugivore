##############################################################################
##
##  Palms: Spatial multiple-predictor models for islands and mainlands (full dataset)
##
##  W.D. Kissling 
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
#palms<-read.table("TDWG_Environment_AllData_Reduced.txt", header=T, sep="\t")       #366 units, Antarctica and Marcus Island not included
#str(palms)

palms<-read.table("TDWG_Environment_AllData_2010_10_Reduced.txt", header=T, sep="\t")       #366 units, Antarctica and Marcus Island not included
str(palms)



#------------------------------- DATA -----------------------------------------#

#Full dataset
palms_full<-subset(palms, palms$PALMSR>0)                                     #194 units with palms
str(palms_full)
palms_full$REALM_LONG<-factor(palms_full$REALM_LONG)
str(palms_full)

#Continents
palms_continent<-palms_full[palms_full$ISISLAND==0,]                      #128 mainlands
str(palms_continent)

#Island full
palms_islands<-palms_full[palms_full$ISISLAND==1&palms_full$PALMSR>1,]                        #66 islands
str(palms_islands)


#----------------------------- COORDINATES ------------------------------------#

######################
#Coordinates full dataset
coords_full<-cbind(palms_full$LONG,palms_full$LAT)
coords_full<-as.matrix(coords_full)
str(coords_full)
print(coords_full)

######################
#Coordinates continents
coords_cont<-cbind(palms_continent$LONG,palms_continent$LAT)
coords_cont<-as.matrix(coords_cont)
str(coords_cont)
print(coords_cont)

######################
#Coordinates islands
coords_isl<-cbind(palms_islands$LONG,palms_islands$LAT)
coords_isl<-as.matrix(coords_isl)
str(coords_isl)
print(coords_isl)


#-------------------------- NEIGHBOURHOODS ------------------------------------#

######################
#1 neighbour full
nbk1_full<-knn2nb(knearneigh(coords_full, k=1, longlat = T))
summary(nbk1_full)
plot(nbk1_full, coords_full, pch=20)
nbk1_full_sym <- make.sym.nb(nbk1_full)    #make symmetric
sw1_full <- nb2listw(nbk1_full_sym)
plot(sw1_full,coords_full, pch=20)

######################
#identify max distance
dsts1 <- unlist(nbdists(nbk1_full, coords_full, longlat = T))     #get 
summary(dsts1)
max(dsts1)

######################
#all sites with at least 1 neighbour (3000 km)
nbk_3500<-dnearneigh(coords_full, d1=0, d2=3500, longlat = T)
summary(nbk_3500, coords_full, pch=20)
nbk_3500_sym <- make.sym.nb(nbk_3500)    #make symmetric
sw_3500 <- nb2listw(nbk_3500_sym)
summary(sw_3500,coords_full, pch=20)


######################
#1 neighbour cont
nbk1_cont<-knn2nb(knearneigh(coords_cont, k=1, longlat = T))
summary(nbk1_cont)
plot(nbk1_cont, coords_cont, pch=20)
nbk1_cont_sym <- make.sym.nb(nbk1_cont)    #make symmetric
sw1_cont <- nb2listw(nbk1_cont_sym)
plot(sw1_cont,coords_cont, pch=20)

######################
#1 neighbour isl
nbk1_isl<-knn2nb(knearneigh(coords_isl, k=1, longlat = T))
summary(nbk1_isl)
plot(nbk1_isl, coords_isl, pch=20)
nbk1_isl_sym <- make.sym.nb(nbk1_isl)    #make symmetric
sw1_isl <- nb2listw(nbk1_isl_sym)
plot(sw1_isl,coords_isl, pch=20)

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
sw_ivd <- nb2listw(nbk_all, glist = idw)


#-------- ALL SITES USING SELECTED VARIABLES WITH SPEARMAN < 0.60 -------------#

#Data
str(palms_full)

#Global model
global_glm<-lm(log(palms_full$PALMSR)~
log(palms_full$PREC_Sum)+
log(palms_full$PREC_CV)+
log(palms_full$Tmean_mean/10+273.15)+
log(palms_full$alt_range)+
log(palms_full$AREA_KM2)+
log(abs(palms_full$ensLGM_Pano))+log(palms_full$ensLGM_Tano)+
as.factor(palms_full$REALM_LONG)+as.factor(palms_full$ISISLAND))

summary(global_glm)
step(global_glm)

#Selected
selected_glm<-lm(log(palms_full$PALMSR) ~ log(palms_full$PREC_Sum)+
log(palms_full$Tmean_mean/10 + 273.15) + log(palms_full$alt_range)+
log(palms_full$AREA_KM2) + log(palms_full$ensLGM_Tano)+
as.factor(palms_full$REALM_LONG))

summary(selected_glm)

#Spatial model
sem_all<-errorsarlm(selected_glm, listw=sw1_full, na.action=na.omit)
summary(sem_all)

#Standardized coefficients
sem_all$coef[2]*sd(log(palms_full$PREC_Sum))/sd(log(palms_full$PALMSR))
sem_all$coef[3]*sd(log(palms_full$Tmean_mean/10 + 273.15))/sd(log(palms_full$PALMSR))
sem_all$coef[4]*sd(log(palms_full$alt_range))/sd(log(palms_full$PALMSR))
sem_all$coef[5]*sd(log(palms_full$AREA_KM2))/sd(log(palms_full$PALMSR))
sem_all$coef[6]*sd(log(palms_full$ensLGM_Tano))/sd(log(palms_full$PALMSR))

#Model fit
sem_all_fit <- print.sarlm.pred(predict.sarlm(sem_all))
cor(x=sem_all_fit$fit, y=log(palms_full$PALMSR), method="pearson")^2     #Full model fit
cor(x=sem_all_fit$trend, y=log(palms_full$PALMSR), method="pearson")^2   #Environment only
summary(sem_all)
moran.mc(residuals(sem_all), listw=sw_3500, alternative = "less", nsim=1000)


