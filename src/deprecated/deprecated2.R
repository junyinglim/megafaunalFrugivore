
z$type <- "Current"
z2$type<- "Present-natural"
z3<- rbind(z,z2)
z4 <- ggplot(aes(y = presid, x =var,color=type),data = z3)+ geom_point(shape = 1) + geom_smooth(method = "lm") + labs(y = "Mean fruit size\n(Partial residuals)", x = "Scaled log median body size") + theme(panel.background = element_blank(), legend.justification = "center", legend.title = element_blank(), legend.position = c(0.8, 0.1)) + scale_colour_manual(values = wes_palette("Cavalcanti1",2)) 
ggsave(z4, filename = file.path(fig.dir, "presid.pdf"), width = 5, height =4.5)

nw_curr_mod_avg_stat$treatment <- "Current"
nw_pnat_mod_avg_stat$treatment <- "Present-natural"
nw_all_mod_avg_stat <- rbind(nw_curr_mod_avg_stat,nw_pnat_mod_avg_stat)

newlabel <- data.frame(old = c("scale(ensLGM_Pano)","scale(ensLGM_Tano)", "scale(log(curr_medianBodySize))", "scale(NPP_mean)", "scale(PC1)", "scale(PC2)", "scale(PC3)", "scale(log(presNat_medianBodySize))"), new = c("LGM Prec anomaly", "LGM Temp anomaly", "Median Mammal Body Size", "NPP", "PC1", "PC2", "PC3", "Median Mammal Body Size"))
nw_all_mod_avg_stat$coefficient_label <- newlabel$new[match(nw_all_mod_avg_stat$coefficient, newlabel$old)]
nw_all_mod_avg_stat$coefficient_label <- factor(nw_all_mod_avg_stat$coefficient_label, level = c("Median Mammal Body Size", "NPP", "PC1", "PC2", "PC3", "LGM Prec anomaly", "LGM Temp anomaly"))

nw_relimpt_plot <- ggplot(data = nw_all_mod_avg_stat) + geom_point(aes(y = fullAvgCoef, x= coefficient_label, size = importance, color = treatment),position=position_dodge(width=0.5)) + geom_errorbar(aes(ymin = lower2.5, ymax = upper97.5, x = coefficient_label, color = treatment),position=position_dodge(width=0.5), width = 0.2) + labs(y = "Coefficients") + theme(axis.text.x = element_text(angle = 45, hjust = 0.9), axis.title.x = element_blank()) + scale_size_continuous(breaks = seq(0,1,0.25), limits = c(0,1), name = "Variable Importance") + scale_colour_manual(name = NULL, values = wes_palette("Cavalcanti1",2))

ggsave(nw_relimpt_plot, filename = file.path(fig.dir,"nw_relimpt_plot.pdf"), width = 7, height = 4)




# Present natural, global PCA
target.col <- c("PC1", "PC2","PC3", "meanFruitLengthFilled", "curr_medianBodySize", "NPP_mean", "PC1", "PC2", "PC3", "ensLGM_Tano", "ensLGM_Pano")
tdwg_final3 <- tdwg_final2[complete.cases(tdwg_final2[target.col]),]
fullmod <- lm(log(meanFruitLengthFilled) ~ log(curr_medianBodySize) + NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tano + ensLGM_Pano, data = tdwg_final3, na.action = "na.fail")
fullmod_dr <- dredge(fullmod)
fullmod_avg <- model.avg(fullmod_dr)
summary(fullmod_avg)











fullmod_step <- step(fullmod)
crPlot(fullmod, variable = "log(curr_medianBodySize)", smooth = FALSE, main = "Global ", ylab = "Partial Residuals (log mean fruit size)", xlab = "Log Current Median Body Size" )
summary(fullmod_step)

# Present natural global, global PCA
fullmod2 <- lm(log(medianFruitLengthFilled) ~ log(presNat_medianBodySize) + NPP_mean + PC1 + PC2 + PC3 + log(AREA_KM2) + ensLGM_Tano + ensLGM_Pano, data = tdwg_final2)
crPlot(fullmod2, variable = "log(presNat_medianBodySize)", smooth = FALSE, main = "Global ", ylab = "Partial Residuals (log median fruit size)", xlab = "LogCURRENT Median Body Size")
summary(fullmod2)

# New world only PCA
tdwg_final2_NW <- subset(tdwg_final2, THREEREALM == "NewWorld")
nwPCA <- prcomp(x = tdwg_final2_NW[env.var], scale. = T, center = T)
tdwg_final2_NW$PC1 <- nwPCA$x[,1]
tdwg_final2_NW$PC2 <- nwPCA$x[,2]
tdwg_final2_NW$PC3 <- nwPCA$x[,3]

# Current NW, NW PCA 
nwfullmod1 <- lm(log(meanFruitLengthFilled) ~ log(curr_medianBodySize) + NPP_mean + PC1 + PC2 + PC3+ ensLGM_Tano + ensLGM_Pano, data = tdwg_final2_NW)
step(nwfullmod1) # not significant after step-wise deletion
crPlot(nwfullmod1, variable = "log(curr_meanBodySize)", smooth = FALSE, main = "New World only ", ylab = "Partial Residuals (log mean fruit size)", xlab = "Log Current Mean Body Size")
nwmod1 <- lm(log(meanFruitLengthFilled) ~ log(curr_meanBodySize) + NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tano+ ensLGM_Pano, data = tdwg_final2_NW, na.action = 'na.fail')
nw_curr_modavg <- model.avg(dredge(nwmod1))
summary(nw_curr_modavg)

# Present natural NW, NW PCA
nwmod2 <- lm(log(meanFruitLengthFilled) ~ log(presNat_meanBodySize) + NPP_mean + PC1 + PC2 + PC3 + ensLGM_Tano+ ensLGM_Pano, data = tdwg_final2_NW, na.action = 'na.fail')
crPlot(nwmod2, variable = "log(presNat_meanBodySize)", smooth = FALSE, main = "New World only ", ylab = "Partial Residuals (log mean fruit size)", xlab = "Log Present-natural Mean Body Size")
summary(nwmod2)
summary(model.avg(dredge(nwmod2)))

# Current NW, Global PCA
tdwg_final2_NWsubset <- subset(tdwg_final2, THREEREALM == "NewWorld")
#scale(log(medianFruitLengthFilled)) 
nwsubset_mod1 <- lm(propMegaPalm~ scale(log(presNat_medianBodySize)) + scale(NPP_mean) + scale(PC1) + scale(PC2) + scale(PC3)+ scale(ensLGM_Tano) + scale(ensLGM_Pano), data = subset(tdwg_final2, propMegaPalm >0 & THREEREALM == "NewWorld"))
summary(nwsubset_mod1)
crPlot(nwsubset_mod1, variable = "scale(log(presNat_medianBodySize))", smooth = FALSE, main = "New World only ", ylab = "Partial Residuals (log mean fruit size)", xlab = "Log Current Mean Body Size")

nwsubset_mod2 <- lm(scale(log(meanFruitLengthFilled)) ~ scale(log(presNat_medianBodySize)) + scale(NPP_mean) + scale(PC1) + scale(PC2) + scale(PC3)+ scale(ensLGM_Tano) + scale(ensLGM_Pano), data = tdwg_final2_NWsubset, na.action = "na.fail")
summary(nwsubset_mod2)
crPlot(nwsubset_mod1, variable = "log(curr_medianBodySize)", smooth = FALSE, main = "New World only ", ylab = "Partial Residuals (log mean fruit size)", xlab = "Log Current Mean Body Size")

summary(model.avg(dredge(nwsubset_mod2)))
confint(model.avg(dredge(nwsubset_mod2)))

nwsubset_mod2 <- lm(log(meanFruitLengthFilled) ~ log(presNat_meanBodySize) + NPP_mean + PC1 + PC2 + PC3+ ensLGM_Tano + ensLGM_Pano, data = tdwg_final2_NWsubset)
summary(nwsubset_mod2)
crPlot(nwsubset_mod2, variable = "log(presNat_meanBodySize)", smooth = FALSE, main = "New World only ", ylab = "Partial Residuals (log mean fruit size)", xlab = "Log Present-natural Mean Body Size")

# Change in median size
deltaSize_mod <- lm(log(meanFruitLengthFilled) ~ deltaMedianBodySize + NPP_mean + PC1 + PC2 + PC3+ ensLGM_Tano + ensLGM_Pano, data = tdwg_final2)
crPlot(deltaSize_mod, variable = "deltaMedianBodySize", smooth = FALSE, main = "New World only ", ylab = "Partial Residuals (log mean fruit size)", xlab = "Log Present-natural Mean Body Size")
summary(deltaSize_mod)

plot(deltaMedianBodySize ~ log(curr_medianBodySize), data = tdwg_final2)
plot(meanFruitLengthFilled ~ curr_meanBodySize, data = tdwg_final2)




test <- melt(tdwg_final, id.var = c("LEVEL_3_CO", "THREEREALM", "Geology_3categ","meanFruitLengthFilled", "megapalm_nsp","propMegaPalm", "AREA_KM2"), measure.vars = c("presNat_meanBodySize", "curr_meanBodySize", "presNat_megaHerb_nSp", "curr_megaHerb_nSp" , "propMegaMam_curr", "propMegaMam_presnat"))

ggplot(aes(y = meanFruitLengthFilled, x= presNat_meanBod)) + geom_point()




fruitSize_mammalBodySizePlot  <- ggplot(aes(y = propMegaPalm, x = value), data = subset(test, variable %in% c("propMegaMam_curr", "propMegaMam_presnat") & !LEVEL_3_CO %in% c("VNA", "MRQ") & Geology_3categ %in% c("continental", "mainland") & propMegaPalm > 0)) +
  geom_point() +
  #facet_wrap(~variable) +
  facet_wrap(variable~THREEREALM) +
  geom_smooth(method = "glm")
ggsave(fruitSize_mammalBodySizePlot, file = file.path(fig.dir, "fruitVsMammalPlot.pdf"))

summary(glm(megapalm_nsp ~ log(AREA_KM2), data = test, family = "poisson"))

tdwg_final$mamSdiff <- with(presNat_nSp-curr_nSp, data = tdwg_final)
tdwg_final$mamMdiff <- with(presNat_meanBodySize-curr_meanBodySize, data = tdwg_final)

summary(lm(meanFruitLengthFilled ~ curr_nSp, data = subset(tdwg_final, THREEREALM == "NewWorld")))

library(MuMIn)
tdwg_final_subset <- na.omit(tdwg_final[c("propMegaPalm", "megapalm_nsp", "meanFruitLengthFilled", "curr_meanBodySize", "presNat_meanBodySize", "NPP_mean", "PREC_Sum", "Tmean_mean", "soilcount", "Geology_3categ", "alt_range", "AREA_KM2", "mamMdiff", "mamSdiff")])
#tdwg_final_subset2 <- subset(tdwg_final_subset, Geology_3categ %in% c("mainland", "continental")) # 154 datapoints
fullmod <- lm(log(meanFruitLengthFilled) ~ curr_meanBodySize + presNat_meanBodySize + NPP_mean + PREC_Sum + Tmean_mean + soilcount + alt_range + mamMdiff + mamSdiff, data = tdwg_final_subset2, na.action = "na.fail")
moddredge <- dredge(fullmod, beta = "sd")
modavg <- model.avg(moddredge)
confint(modavg)
summary(modavg)


fullmod <- lm(propMegaPalm ~ curr_meanBodySize + presNat_meanBodySize + NPP_mean + PREC_Sum + Tmean_mean + soilcount + mamMdiff + mamSdiff, data = tdwg_final_subset, na.action = "na.fail")
moddredge <- dredge(fullmod, beta = "sd")
modavg <- model.avg(moddredge)
confint(modavg)
summary(modavg)

fullmod <- lm(log(megapalm_nsp) ~ curr_meanBodySize + presNat_meanBodySize + NPP_mean + PREC_Sum + Tmean_mean + soilcount + log(AREA_KM2) + mamMdiff + mamSdiff, data = subset(tdwg_final_subset2, megapalm_nsp >0), na.action = "na.fail", family = "poisson")
moddredge <- dredge(fullmod, beta = "sd")
modavg <- model.avg(moddredge)
confint(modavg)
summary(modavg)
ggplot() + geom_point(aes(y = log(meanFruitLengthFilled), x = mamMdiff, color = THREEREALM), data = tdwg_final)


hist(resid(fullmod))
plot(resid(fullmod)~ fitted(fullmod))
abline(0,0)
library(lme4); library(lmerTest)
summary(lm(FruitSizeStd~mamMdiff, data = tdwg_final))
plot(propMegaPalm~mamMdiff, data = subset(tdwg_final_subset2, megapalm_nsp >0))
plot(log(meanFruitLengthFilled)~mamSdiff, data = subset(tdwg_final_subset2))
boxplot(FruitSizeStd~THREEREALM, tdwg_final)

meanFruitByRealm <- tapply(tdwg_final$meanFruitLengthFilled, INDEX = tdwg_final$THREEREALM, FUN = mean)
meanFruitByRealm <- data.frame(meanFruitByRealm, THREEREALM = names(meanFruitByRealm))
tdwg_final <- merge(tdwg_final, meanFruitByRealm)
tdwg_final$FruitSizeStd <- (tdwg_final$meanFruitLengthFilled - tdwg_final$meanFruitByRealm) / sd(tdwg_final$meanFruitByRealm)


fruitSize_mammalBodySizePlot_noisland  <- ggplot(aes(y = log(meanFruitLengthFilled), x = log(value)), data = subset(test, Geology_3categ %in% c("continental", "mainland"))) +
  geom_point() +
  facet_wrap(variable~THREEREALM) +
  geom_smooth(method = "lm") 

ggsave(fruitSize_mammalBodySizePlot_noisland, file = file.path(fig.dir, "fruitVsMammalPlot_noisland.pdf"))



summary(lm(log(meanFruitLengthFilled)~ mamMdiff, data = subset(tdwg_final, Geology_3categ %in% c("continental", "mainland"))))

summary(lm(log(meanFruitLengthFilled)~ log(presNat_meanBodySize), data = subset(tdwg_final, THREEREALM == "NewWorld" & ! LEVEL_3_CO %in% c("VNA", "MRQ"))))

mammalCurrVsPresNatBodySize <- ggplot(data = tdwg_final) + geom_point(aes(y = log(presNat_meanBodySize), x = log(curr_meanBodySize), color = THREEREALM)) + geom_abline(aes(intercept = 0, slope = 1))
ggsave(mammalCurrVsPresNatBodySize, file = file.path(fig.dir, "mammalCurrVsPresNatBodySize.pdf"))


ggplot(aes(x = log(curr_meanBodySize), y = meanFruitLengthFilled), data = tdwg_final) + geom_point() + geom_smooth(method = "lm") 
ggplot(aes(x = log(presNat_meanBodySize), y = log(meanFruitLengthFilled)), data = tdwg_final) + geom_point() + geom_smooth(method = "lm") 

hist(log(tdwg_final$meanFruitLengthFilled))
hist(tdwg_final$meanFruitLengthFilled)

# Calculate the proportion of mammals that are megafaunal, and the number of species of megafaunal fruits










## ISLAND TAXA
table(tdwg_env$GeologicalOrigin)
table(tdwg_env$Geology_3categ)
ggplot(aes(y = log(meanFruitLengthFilled), x=log(DistToCont_km), colour = Geology_3categ), data = subset(tdwg_final, Geology_3categ %in% c("continental", "mainland", "volcanic"))) + geom_point() + geom_smooth(method = "lm") + geom_label_repel(aes(y = log(meanFruitLengthFilled), x = log(DistToCont_km), label = LEVEL_3_CO))

library(ggrepel)

names(tdwg_env)
ggplot(aes(y = log(meanFruitLengthFilled), x=THREEREALM, color = Geology_3categ), data = tdwg_final) + geom_boxplot()
names(tdwg_env)

library(geosphere)

mat <- distm(tdwg_final[,c('LONG','LAT')], tdwg_final[,c('LONG','LAT')], fun=distGeo)
nearestTDWGlist <- tdwg_final$LEVEL_3_CO[unlist(apply(mat, MARGIN = 1, FUN = function(x){ which(min(x[x>0]) == x) }))]
distNearestTDWGlist <- unlist(apply(mat, MARGIN = 1, FUN = function(x){ min(x[x>0]) }))
tdwg_final$nearestTDWG <- nearestTDWGlist
tdwg_final$distNearestTDWG <- distNearestTDWGlist

ggplot(aes(y = log(meanFruitLengthFilled), x=distNearestTDWG, color = Geology_3categ), data = tdwg_final)
names(tdwg_env)

