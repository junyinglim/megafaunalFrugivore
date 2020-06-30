## Estimate extinction rates ##
# Estimate extinction probabilities of extant mammals using observed changes in red list status over observed time frames
# Author: Jun Ying Lim
# Reference: Lim, J.Y., Svenning, J.-C., Göldel, B., Faurby, S. & Kissling, W.D. Past and future extinctions shape the body size - fruit size relationship between palms and mammalian frugivores.

## Packages and directories ========================
library(cowplot)
library(ggplot2)
rm(list = ls())
main.dir <- "/Users/junyinglim/Dropbox/Projects/2019/palms/projects/megafaunalFrugivore/src/R2_clean/"
data.dir <- file.path(main.dir, "data")
res.dir <- file.path(main.dir, "results")
fig.dir <- file.path(main.dir, "figs")

## Input empirical red list status change data ========================
set.seed(12345)

# Di Marco et al (2014) data was obtained by digitizing Table 1
dimarcodata <- matrix( c(262, 24, 15, 6, 3, 0,
                         3, 20, 18, 4, 1, 0,
                         1, 4, 52, 28, 3, 1,
                         1, 0, 2, 27, 10, 0,
                         0, 0, 1, 2, 9, 1,
                         0, 0, 0, 0, 0, 0), byrow = T, ncol = 6)

## Hoffmann et al (2011) data
# Number of species across categories at the start and end of each time point
# obtained by cross-referencing Table 1, Table S5 and Table S6 in the Supplementary
# Species identities of those that changed category are in Table S6
# Total number of species in each category for the two time periods are in Table S5

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

## Define model and likelihood function ========================
redlisttrans <- function(LC_NT, NT_LC, NT_VU, 
                         VU_EN, VU_NT, EN_CR,
                         EN_VU, CR_EX, CR_EN,
                         EX_CR, t){
  # Generates transition probabilities assuming a CTMC model
  # representing instantaneous rates of change between red list categories
  # Args:
  #     LC_NT = transition rate between LC and NT category
  #     NT_LC = transition rate between NT and LC category
  #     NT_VU = transition rate between NT and VU category
  #     VU_EN = transition rate between VU and EN category
  #     EN_CR = transition rate between EN and CR category
  #     EN_VU = transition rate between EN and VU category
  #     CR_EX = transition rate between CR and EX category
  #     CR_EN = transition rate between CR and EN category
  #     EX_CR = transition rate between EX and CR category
  #     t = amount of time
  #
  # Returns:
  #     matrix, transition probabilities for specified amount of time
  
  # Ensure rows sum to zero
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
  # Arguments:
  #    data: observed data
  #    par: parameter values to pass to redlisttrans
  # Returns:
  #   numeric, negative log-likelihood
  
  if(is.null(EX_CR)){
    P <- redlisttrans(par[1], par[2], par[3], par[4], par[5], par[6], par[7], par[8], par[9], par[10], ...)  
  } else {
    P <- redlisttrans(par[1], par[2], par[3], par[4], par[5], par[6], par[7], par[8], par[9], EX_CR = EX_CR, ...)  
  }
  
  P_s <- P / rowSums(P) # standardize rows
  
  return(-1 * sum(log(P_s[data > 0]) * data[data>0]))
}

## Estimate empirical rates of change in red list status and extinction probabilities ========================
# Fit CTMC model to empirical data using maximum likelihood (ML) to obtain ML-estimates of rates of transition between red list categories
init_par <- rep(0.01, 10)
dimarcoML <- optim(par = init_par[1:9], fn = redlistLik, method = "L-BFGS-B",
                   lower = 1E-6, upper = 1, data = dimarcodata, t = 33, EX_CR = 0)

hoffmannML <- optim(par = init_par, fn = redlistLik, method = "L-BFGS-B",
                    lower = 1E-6, upper = 1, data = hoffmanndata, t = 12)

# Using ML-estimated rates, compute transition probabilities assuming a time frame of 100 years
dimarcoML_pExt <- redlisttrans(dimarcoML$par[1], dimarcoML$par[2], dimarcoML$par[3],
                               dimarcoML$par[4], dimarcoML$par[5], dimarcoML$par[6],
                               dimarcoML$par[7], dimarcoML$par[8], dimarcoML$par[9],
                               EX_CR = 0, t= 100)

hoffmannML_pExt <- redlisttrans(hoffmannML$par[1], hoffmannML$par[2], hoffmannML$par[3],
                                hoffmannML$par[4], hoffmannML$par[5], hoffmannML$par[6],
                                hoffmannML$par[7], hoffmannML$par[8], hoffmannML$par[9],
                                hoffmannML$par[10], t= 100)

# Extinction probabilities for each red list category is simply the last column of the transition probability matrix
ctmc_pExt <- data.frame(IUCN.Status = c("LC", "NT", "VU", "EN", "CR"),
                        dimarco = dimarcoML_pExt[1:5,6],
                        hoffmann = hoffmannML_pExt[1:5,6],
                        davis = c(0.0017, 0.0141, 0.1000, 0.6723, 0.9990))

# Davis' rates uses the IUCN 2001 Criterion E definitions for threatened categories (CR, EN and VU)
# i.e., P(ext_CR, 10 years) = 0.5; P(ext_EN, 20 years) = 0.2 and P(ext_VU, 100 years) = 0.1
# Rearranging the equations P(r, t) = 1 - exp(-rt) , the rates for CR = -log(0.5)/10, EN = -log(0.8)/20, VU = -log(0.9) / 100
# However, Davis estimated them by assuming an exponential decrease in extinction rate across IUCN categories as an ordinal variable
# (i.e., CR = 1, EN = 2, VU = 3, NT = 4, LC = 5)
# Then assuming a constant rate through time, they calculated the extinction probability over a hundred year time frame
# as follows P(t) = 1 - exp(-rate * time)

write.csv(ctmc_pExt, file.path(data.dir, "ctmc_pExt.csv"),row.names = FALSE)

## Parametric bootstrapping of ML-estimated rates to test for robustness ========================

# Compute red list category transition probabilities under the same observed timeframe
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

# Simulate red list status changes over the same observed time frame and the same starting distribution of species across red lsit categories
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

# Estimate using ML the rates of transition between red list statuses using simulated red list changes under the empirical rates
dimarco_sim_MLpara <- list()
hoffmann_sim_MLpara <- list()
for(i in 1:1000){
  pb = txtProgressBar(min = 0, max = 1000, style = 3, initial = 0) 
  setTxtProgressBar(pb, i)
  dimarco_sim_MLpara[[i]] <- 
    optim(par = init_par[1:9], fn = redlistLik, method = "L-BFGS-B",
          lower = 1E-6, upper = 1, data = dimarco_sim[[i]], EX_CR = 0, t = 33)$par
  
  hoffmann_sim_MLpara[[i]] <- 
    optim(par = init_par, fn = redlistLik, method = "L-BFGS-B",
          lower = 1E-6, upper = 1, data = hoffmann_sim[[i]], t = 12)$par
}
close(pb)

# Clean up data and save a copy of the results
par_labels <- c("LC to NT", "NT to LC", "NT to VU", "VU to EN", "VU to NT", "EN to CR", "EN to VU", "CR to EX", "CR to EN", "EX to CR")
dimarco_par_sim <- data.frame( value = unlist(dimarco_sim_MLpara) ) 
dimarco_par_sim$par_label <- rep(par_labels[1:9], 1000)

dimarco_par_obs <- data.frame( value = dimarcoML$par)
dimarco_par_obs$par_label <- par_labels[1:9]

hoffmann_par_sim <- data.frame( value  = unlist(hoffmann_sim_MLpara) )
hoffmann_par_sim$par_label <- rep(par_labels[1:10], 1000)

hoffmann_par_obs <- data.frame(value = hoffmannML$par)
hoffmann_par_obs$par_label <- par_labels

dimarco_parboot <- list(dimarco_par_sim, dimarco_par_obs)
saveRDS(dimarco_parboot, file.path(res.dir, "dimarco_parboot.rds"))
hoffmann_parboot <- list(hoffmann_par_sim, hoffmann_par_obs)
saveRDS(hoffmann_parboot, file.path(res.dir, "hoffmann_parboot.rds"))

## Generate human readable files  ========================
dimarco_parboot <- readRDS(file.path(res.dir, "dimarco_parboot.rds"))[[1]]
dimarco_parboot$iteration <- rep(1:1000, each = 9)
write.csv(dimarco_parboot, file.path(res.dir, "dimarco_parboot.csv"), row.names= F)

dimarco_par <- readRDS(file.path(res.dir, "dimarco_parboot.rds"))[[2]]
write.csv(dimarco_par, file.path(res.dir, "dimarco_par.csv"), row.names= F)

hoffmann_parboot <- readRDS(file.path(res.dir, "hoffmann_parboot.rds"))[[1]]
hoffmann_parboot$iteration <- rep(1:1000, each = 10)
write.csv(hoffmann_parboot, file.path(res.dir, "hoffmann_parboot.csv"), row.names= F)

hoffmann_par <- readRDS(file.path(res.dir, "hoffmann_parboot.rds"))[[2]]
write.csv(hoffmann_par, file.path(res.dir, "hoffmann_par.csv"), row.names= F)