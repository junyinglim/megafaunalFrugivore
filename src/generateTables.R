## Generate tidy tables for result reporting
# Author: Jun Ying Lim

# Directories 
rm(list = ls())
main.dir <- "/Users/junyinglim/Dropbox/Projects/2019/palms/projects/megafaunalFrugivore"
data.dir <- file.path(main.dir, "data")
res.dir <- file.path(main.dir, "results")
src.dir <- file.path(main.dir, "src")
source(file.path(src.dir, "ggmodavg.R"))

# Define labels, and order of table
geographic_levels <- c("Global","Afrotropics","Neotropics","Indotropics")
med_levels <- c("(Intercept)", "curr_logMedBS_scl","pnat_logMedBS_scl", "pnat_logMedBS_cons_scl", "pnat_logMedBS_lib_scl","globalPC1_scl", "globalPC2_scl","globalPC3_scl", "regionalPC1_scl","regionalPC2_scl", "regionalPC3_scl","lgm_ens_Pano_scl", "lgm_ens_Tano_scl")
med_labels <- c("Intercept", "Log median body size", "Log median body size", "Log median body size", "Log median body size", "Climate PC1", "Climate PC2", "Climate PC3", "Climate PC1","Climate PC2", "Climate PC3","LGM Prec. Anom.", "LGM Temp. Anom.")
max_levels <- c("(Intercept)", "curr_logMax95BS_scl","pnat_logMax95BS_scl", "pnat_logMax95BS_cons_scl", "pnat_logMax95BS_lib_scl", "globalPC1_scl", "globalPC2_scl","globalPC3_scl", "regionalPC1_scl","regionalPC2_scl", "regionalPC3_scl","lgm_ens_Pano_scl", "lgm_ens_Tano_scl")
max_labels <- c("Intercept", "Log maximum body size", "Log maximum body size", "Log maximum body size", "Log maximum body size", "Climate PC1", "Climate PC2", "Climate PC3", "Climate PC1","Climate PC2", "Climate PC3","LGM Prec. Anom.", "LGM Temp. Anom.")

target_col <- c("Geographic.Scale", "Scenario", "coefficient", "fullAvgCoef", "fullAvgSE", "ConfidenceInterval","varexp", "totalR2")

# OLS Median
medBS_ols_modavg <- roundNumbers(read.csv(file.path(res.dir, "medBS_ols_modavg.csv")))
medBS_ols_modavg <- subset(medBS_ols_modavg, !coefficient == "(Intercept)" & Geographic.Scale %in% geographic_levels)
medBS_ols_modavg$Geographic.Scale <- factor(medBS_ols_modavg$Geographic.Scale,
                                            levels = c("Global", "Afrotropics", "Neotropics", "Indotropics"))
medBS_ols_modavg$coefficient <- factor(medBS_ols_modavg$coefficient,
                                       levels = med_levels,
                                       labels = med_labels)
medBS_ols_modavg2 <- with(medBS_ols_modavg, medBS_ols_modavg[order(Geographic.Scale, Scenario, coefficient),])
medBS_ols_modavg2 <- format(medBS_ols_modavg2)
medBS_ols_modavg2$ConfidenceInterval <- paste0("(", medBS_ols_modavg2$fulllower2.5, "; ",
                                               medBS_ols_modavg2$fullupper97.5, ")")
medBS_ols_modavg2$varexp[as.numeric(medBS_ols_modavg2$varexp) < 0.001] <- "< 0.001"
medBS_ols_modavg2$totalR2[!medBS_ols_modavg2$coefficient == "Log median body size"] <- ""
write.csv(subset(medBS_ols_modavg2, (!coefficient == "Intercept"))[target_col],
          file.path(res.dir, "medBS_ols_modavg_clean.csv"), row.names = F, quote = F)

# OLS Median Partial SD
medBS_ols_cade_modavg <- roundNumbers(read.csv(file.path(res.dir, "medBS_ols_cade_modavg.csv")))
medBS_ols_cade_modavg <- subset(medBS_ols_cade_modavg, !coefficient == "(Intercept)" & Geographic.Scale %in% geographic_levels)
medBS_ols_cade_modavg$Geographic.Scale <- factor(medBS_ols_cade_modavg$Geographic.Scale,
                                            levels = c("Global", "Afrotropics", "Neotropics", "Indotropics"))
medBS_ols_cade_modavg$coefficient <- factor(medBS_ols_cade_modavg$coefficient,
                                       levels = med_levels,
                                       labels = med_labels)
medBS_ols_cade_modavg2 <- with(medBS_ols_cade_modavg, medBS_ols_cade_modavg[order(Geographic.Scale, Scenario, coefficient),])
medBS_ols_cade_modavg2 <- format(medBS_ols_cade_modavg2)
medBS_ols_cade_modavg2$ConfidenceInterval <- paste0("(", medBS_ols_cade_modavg2$fulllower2.5, "; ",
                                                    medBS_ols_cade_modavg2$fullupper97.5, ")")
medBS_ols_cade_modavg2$varexp[as.numeric(medBS_ols_cade_modavg2$varexp) < 0.001] <- "< 0.001"
medBS_ols_cade_modavg2$totalR2[!medBS_ols_cade_modavg2$coefficient == "Log median body size"] <- ""
write.csv(medBS_ols_cade_modavg2[target_col],
          file.path(res.dir, "medBS_ols_cade_modavg_clean.csv"), row.names = F, quote = F)


# OLS Maximum BS
maxBS_ols_modavg <- roundNumbers(read.csv(file.path(res.dir, "maxBS_ols_modavg.csv")))
maxBS_ols_modavg <- subset(maxBS_ols_modavg, !coefficient == "(Intercept)" & Geographic.Scale %in% geographic_levels)
maxBS_ols_modavg$Geographic.Scale <- factor(maxBS_ols_modavg$Geographic.Scale,
                                            levels = geographic_levels)
maxBS_ols_modavg$coefficient <- factor(maxBS_ols_modavg$coefficient,
                                       levels = max_levels,
                                       labels = max_labels)
maxBS_ols_modavg2 <- with(maxBS_ols_modavg, maxBS_ols_modavg[order(Geographic.Scale, Scenario, coefficient),])
maxBS_ols_modavg2 <- format(maxBS_ols_modavg2)
maxBS_ols_modavg2$ConfidenceInterval <- paste0("(", maxBS_ols_modavg2$fulllower2.5, "; ",
                                               maxBS_ols_modavg2$fullupper97.5, ")")
maxBS_ols_modavg2$varexp[as.numeric(maxBS_ols_modavg2$varexp) < 0.001] <- "< 0.001"
maxBS_ols_modavg2$totalR2[!maxBS_ols_modavg2$coefficient == "Log maximum body size"] <- ""
write.csv(maxBS_ols_modavg2[target_col],
          file.path(res.dir, "maxBS_ols_modavg_clean.csv"), row.names = F, quote = F)

# OLS Maximum BS Partial SD
maxBS_ols_cade_modavg <- roundNumbers(read.csv(file.path(res.dir, "maxBS_ols_cade_modavg.csv")))
maxBS_ols_cade_modavg <- subset(maxBS_ols_cade_modavg, !coefficient == "(Intercept)" & Geographic.Scale %in% geographic_levels)
maxBS_ols_cade_modavg$Geographic.Scale <- factor(maxBS_ols_cade_modavg$Geographic.Scale,
                                            levels = geographic_levels)
maxBS_ols_cade_modavg$coefficient <- factor(maxBS_ols_cade_modavg$coefficient,
                                       levels = max_levels,
                                       labels = max_labels)
maxBS_ols_cade_modavg2 <- with(maxBS_ols_cade_modavg, maxBS_ols_cade_modavg[order(Geographic.Scale, Scenario, coefficient),])
maxBS_ols_cade_modavg2 <- format(maxBS_ols_cade_modavg2)
maxBS_ols_cade_modavg2$ConfidenceInterval <- paste0("(", maxBS_ols_cade_modavg2$fulllower2.5, "; ",
                                                    maxBS_ols_cade_modavg2$fullupper97.5, ")")
maxBS_ols_cade_modavg2$varexp[as.numeric(maxBS_ols_cade_modavg2$varexp) < 0.001] <- "< 0.001"
maxBS_ols_cade_modavg2$totalR2[!maxBS_ols_cade_modavg2$coefficient == "Log maximum body size"] <- ""
write.csv(maxBS_ols_cade_modavg2[target_col],
          file.path(res.dir, "maxBS_ols_cade_modavg_clean.csv"), row.names = F, quote = F)

# SAR results
maxBS_sar_modavg <- roundNumbers(read.csv(file.path(res.dir, "maxBS_sar_modavg.csv")) )
maxBS_sar_modavg <- subset(maxBS_sar_modavg, !coefficient %in% c("(Intercept)", "lambda"))
maxBS_sar_modavg$Geographic.Scale <- factor(maxBS_sar_modavg$Geographic.Scale,
                                            levels = geographic_levels)
maxBS_sar_modavg$coefficient <- factor(maxBS_sar_modavg$coefficient,
                                       levels = max_levels,
                                       labels = max_labels)
maxBS_sar_modavg2 <- with(maxBS_sar_modavg, maxBS_sar_modavg[order(Geographic.Scale, Scenario, coefficient),])
maxBS_sar_modavg2 <- format(maxBS_sar_modavg2) # prevents any changing of format (e.g., decimal places)
maxBS_sar_modavg2$ConfidenceInterval <- paste0("(", maxBS_sar_modavg2$fulllower2.5, "; ",
                                               maxBS_sar_modavg2$fullupper97.5, ")")
maxBS_sar_modavg2$varexp[as.numeric(maxBS_sar_modavg2$varexp) < 0.001] <- "< 0.001"
maxBS_sar_modavg2$totalR2[!maxBS_sar_modavg2$coefficient == "Log maximum body size"] <- ""
write.csv(maxBS_sar_modavg2[target_col],
          file.path(res.dir, "maxBS_sar_modavg_clean.csv"), row.names = F, quote = F)

medBS_sar_modavg <- roundNumbers(read.csv(file.path(res.dir, "medBS_sar_modavg.csv"))) 
medBS_sar_modavg <- subset(medBS_sar_modavg, !coefficient %in% c("(Intercept)", "lambda"))
medBS_sar_modavg$Geographic.Scale <- factor(medBS_sar_modavg$Geographic.Scale,
                                            levels = geographic_levels)
medBS_sar_modavg$coefficient <- factor(medBS_sar_modavg$coefficient,
                                       levels = med_levels,
                                       labels = med_labels)
medBS_sar_modavg2 <- with(medBS_sar_modavg, medBS_sar_modavg[order(Geographic.Scale, Scenario, coefficient),])
medBS_sar_modavg2 <- format(medBS_sar_modavg2)
medBS_sar_modavg2$ConfidenceInterval <- paste0("(", medBS_sar_modavg2$fulllower2.5, "; ",
                                               medBS_sar_modavg2$fullupper97.5, ")")
medBS_sar_modavg2$varexp[as.numeric(medBS_sar_modavg2$varexp) < 0.001] <- "< 0.001"
medBS_sar_modavg2$totalR2[!medBS_sar_modavg2$coefficient == "Log median body size"] <- ""
write.csv(medBS_sar_modavg2[target_col],
          file.path(res.dir, "medBS_sar_modavg_clean.csv"), row.names = F, quote = F)
