################################################################################
#
# Structural equation modelling for ABMI data: WOODY
#
# May 2012, Daniel Kissling, with help from Jim Grace
#
################################################################################


rm(list=ls())
# Load the dataset
# load("C:\\ABMI & PSPs\\ABMI_Info\\Manuscript_Spatial variation of bird species diversity\\ABMI_03to10_Boreal.rdata")
load("D:/Data/wkissli1/My Documents/Projects/ABMI/1_BorealForest_GuildRichness/Analysis/2_Revision/ABMI/May_2012/SEM_lavaan_BIRDS_ABMI/ABMI_03to10_Boreal.rdata")


#########################
#Load packages
library(lavaan)
# library(ReadImages)

#########################
#The dataset
colnames(BR03to10.sites.Richness.Transformed)
SEM.dat<- BR03to10.sites.Richness.Transformed[,c(6:10,12)]    #with WOODY
str(SEM.dat)
colnames(SEM.dat)

#histograms
hist(SEM.dat$BIRDS)
hist(SEM.dat$TEMP)
hist(SEM.dat$PREC)
hist(SEM.dat$ELEV)
hist(SEM.dat$WOODY)
hist(SEM.dat$HUMAN)


#------------------ SEM analyses for ABMI data --------------------------------
# Criteria for fit measures:
# Model Chi-Square with its df and p-value: prefer p-value greater than 0.05
# Root Mean Square Error of Approximation (RMSEA): prefer lower 90%CI to be < 0.05
# Comparative Fit Index (CFI): prefer value greater than 0.90

# The SEM model with WOODY
# plot(read.jpeg("model_WOODY.jpg"))
model.1 <- 'BIRDS ~ ELEV + TEMP + PREC + HUMAN + WOODY
         TEMP ~ ELEV
         PREC ~ ELEV
			   WOODY ~ ELEV + TEMP + PREC + HUMAN'
sem.fit.1 <- sem(model.1, data=SEM.dat)
summary(sem.fit.1, stand=T, rsq=T, fit.measures=T)    #Chi-square p-value, 90% lower CI of RMSEA, and CFI not good
resid(sem.fit.1, type="standardized")                 #High correlation between TEMP and HUMAN
modindices(sem.fit.1)                                 #High for TEMP ~~ HUMAN  15.555 (also TEMP ~~  ELEV  26.913, but no need to include this as there is already an arrow)

# The same model, but adding TEMP ~~ HUMAN (error covariance)
model.2a <- 'BIRDS ~ ELEV + TEMP + PREC + HUMAN + WOODY
         TEMP ~ ELEV
         PREC ~ ELEV
			   WOODY ~ ELEV + TEMP + PREC + HUMAN
         TEMP ~~  HUMAN'
sem.fit.2a <- sem(model.2a, data=SEM.dat, fixed.x=FALSE)  
summary(sem.fit.2a, stand=T, rsq=T, fit.measures=T)   #Chi-square p-value, 90% lower CI of RMSEA, and CFI not good   
resid(sem.fit.2a, type="standardized")                #Highest correlation between PREC and WOODY
modindices(sem.fit.2a)                        #high modification index TEMP ~~  PREC 15.269       
anova(sem.fit.1,sem.fit.2a)                   #model.2a is better than model.1


# The same model, but adding TEMP ~~ PREC (error covariance)
model.2b <- 'BIRDS ~ ELEV + TEMP + PREC + HUMAN + WOODY
         TEMP ~ ELEV
         PREC ~ ELEV
			   WOODY ~ ELEV + TEMP + PREC + HUMAN
         TEMP ~~ HUMAN
         TEMP ~~ PREC'
sem.fit.2b <- sem(model.2b, data=SEM.dat, fixed.x=FALSE)  
summary(sem.fit.2b, stand=T, rsq=T, fit.measures=T) #Chi-square p-value, 90% lower CI of RMSEA, and CFI not good     
resid(sem.fit.2b, type="standardized")             #highest correlation for ELEV and HUMAN
modindices(sem.fit.2b)                             #modification index ELEV ~~ HUMAN 22.474 (also for PREC ~~ WOODY 54.515, but there is already a link in the model)
anova(sem.fit.2a,sem.fit.2b)                       #model.2b is better than model.2a


# The same model, but adding ELEV ~~ HUMAN (error covariance)
model.2c <- 'BIRDS ~ ELEV + TEMP + PREC + HUMAN + WOODY
         TEMP ~ ELEV
         PREC ~ ELEV
			   WOODY ~ ELEV + TEMP + PREC + HUMAN
         TEMP ~~ HUMAN
         TEMP ~~ PREC
         ELEV ~~ HUMAN'
sem.fit.2c <- sem(model.2c, data=SEM.dat, fixed.x=FALSE)  
summary(sem.fit.2c, stand=T, rsq=T, fit.measures=T)    #Chi-square p-value, 90% lower CI of RMSEA, and CFI now ok
resid(sem.fit.2c, type="standardized")         #now ok
modindices(sem.fit.2c)                               
anova(sem.fit.2b,sem.fit.2c)                   #model.2c is better than model.2b


# The same model, but now removing non-significant PREC from BIRDS model
model.2d <- 'BIRDS ~ ELEV + TEMP + HUMAN + WOODY
         TEMP ~ ELEV
         PREC ~ ELEV
			   WOODY ~ ELEV + TEMP + PREC + HUMAN
         TEMP ~~ HUMAN
         TEMP ~~ PREC
         ELEV ~~ HUMAN'
sem.fit.2d <- sem(model.2d, data=SEM.dat, fixed.x=FALSE)  
summary(sem.fit.2d, stand=T, rsq=T, fit.measures=T)    #looks good: Chi-square P-value, CFI, RMSEA lower 90%CI
resid(sem.fit.2d, type="standardized")         #ok
modindices(sem.fit.2d)                         #ok        
anova(sem.fit.2c,sem.fit.2d)                   #model.2d is similar to model.2c


# The same model, but now removing the most non-significant variable (ELEV) from WOODY model
model.2e <- 'BIRDS ~ ELEV + TEMP + HUMAN + WOODY
         TEMP ~ ELEV
         PREC ~ ELEV
			   WOODY ~ TEMP + PREC + HUMAN
         TEMP ~~ HUMAN
         TEMP ~~ PREC
         ELEV ~~ HUMAN'
sem.fit.2e <- sem(model.2e, data=SEM.dat, fixed.x=FALSE)  
summary(sem.fit.2e, stand=T, rsq=T, fit.measures=T)    #looks good
resid(sem.fit.2e, type="standardized")         #ok
modindices(sem.fit.2e)                         #ok       
anova(sem.fit.2d,sem.fit.2e)                   #model.2e is similar to model.2d


# The same model, but now removing the most non-significant variable (PREC) from WOODY model
model.2f <- 'BIRDS ~ ELEV + TEMP + HUMAN + WOODY
         TEMP ~ ELEV
         PREC ~ ELEV
			   WOODY ~ TEMP + HUMAN
         TEMP ~~ HUMAN
         TEMP ~~ PREC
         ELEV ~~ HUMAN'
sem.fit.2f <- sem(model.2f, data=SEM.dat, fixed.x=FALSE)  
summary(sem.fit.2f, stand=T, rsq=T, fit.measures=T)    #looks good
resid(sem.fit.2f, type="standardized")         #ok
modindices(sem.fit.2f)                         #ok      
anova(sem.fit.2e,sem.fit.2f)                   #model.2e is similar to model.2f

# The same model, but now removing the most non-significant variable (TEMP) from WOODY model
model.2g <- 'BIRDS ~ ELEV + TEMP + HUMAN + WOODY
         TEMP ~ ELEV
         PREC ~ ELEV
			   WOODY ~ HUMAN
         TEMP ~~ HUMAN
         TEMP ~~ PREC
         ELEV ~~ HUMAN'
sem.fit.2g <- sem(model.2g, data=SEM.dat, fixed.x=FALSE)  
summary(sem.fit.2g, stand=T, rsq=T, fit.measures=T)    #looks good, and WOODY ~ HUMAN is marginally significant (p = 0.05)
resid(sem.fit.2g, type="standardized")         #ok
modindices(sem.fit.2g)                       #ok       
anova(sem.fit.2f,sem.fit.2g)                  #model.2e is similar to model.2f

#Final model
plot(read.jpeg("model_WOODY_final_withStdCoeff.jpg"))