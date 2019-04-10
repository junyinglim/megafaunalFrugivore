## Clean up results for manuscript
main.dir <- "/Users/junyinglim/Dropbox/Projects/2019/palms/projects/megafaunalFrugivore"
data.dir <- file.path(main.dir, "data")
res.dir <- file.path(main.dir, "results")
src.dir <- file.path(main.dir, "src")
source(file.path(src.dir, "ggmodavg.R"))

# Define labels, and order of table
geographic_levels <- c("Global","Afrotropics","Neotropics","Indotropics")
med_levels <- c("(Intercept)", "curr_logMedBS_scl","pnat_logMedBS_scl","globalPC1_scl", "globalPC2_scl","globalPC3_scl", "regionalPC1_scl","regionalPC2_scl", "regionalPC3_scl","lgm_ens_Pano_scl", "lgm_ens_Tano_scl")
med_labels <- c("Intercept", "Log median body size", "Log median body size", "Climate PC1", "Climate PC2", "Climate PC3", "Climate PC1","Climate PC2", "Climate PC3","LGM Prec. Anom.", "LGM Temp. Anom.")
max_levels <- c("(Intercept)", "curr_logMax95BS_scl","pnat_logMax95BS_scl","globalPC1_scl", "globalPC2_scl","globalPC3_scl", "regionalPC1_scl","regionalPC2_scl", "regionalPC3_scl","lgm_ens_Pano_scl", "lgm_ens_Tano_scl")
max_labels <- c("Intercept", "Log maximum body size", "Log maximum body size", "Climate PC1", "Climate PC2", "Climate PC3", "Climate PC1","Climate PC2", "Climate PC3","LGM Prec. Anom.", "LGM Temp. Anom.")

# OLS Median
medBS_ols_modavg <- read.csv(file.path(res.dir, "medBS_ols_modavg.csv"))
medBS_ols_modavg$Geographic.Scale <- factor(medBS_ols_modavg$Geographic.Scale,
                                            levels = c("Global", "Afrotropics", "Neotropics", "Indotropics"))
medBS_ols_modavg$coefficient <- factor(medBS_ols_modavg$coefficient,
                                       levels = med_levels,
                                       labels = med_labels)
medBS_ols_modavg2 <- with(medBS_ols_modavg, medBS_ols_modavg[order(Geographic.Scale, Scenario, coefficient),])
medBS_ols_modavg2$ConfidenceInterval <- paste0("(", medBS_ols_modavg2$lower2.5, ", ",
                                               medBS_ols_modavg2$upper97.5, ")")
write.csv(subset(medBS_ols_modavg2, (!coefficient == "Intercept")),
          file.path(res.dir, "medBS_ols_modavg_clean.csv"), row.names = F)

# OLS Maximum BS
maxBS_ols_modavg <- read.csv(file.path(res.dir, "maxBS_ols_modavg.csv"))
maxBS_ols_modavg$Geographic.Scale <- factor(maxBS_ols_modavg$Geographic.Scale,
                                            levels = geographic_levels)
maxBS_ols_modavg$coefficient <- factor(maxBS_ols_modavg$coefficient,
                                       levels = max_levels,
                                       labels = max_labels)
maxBS_ols_modavg2 <- with(maxBS_ols_modavg, maxBS_ols_modavg[order(Geographic.Scale, Scenario, coefficient),])
maxBS_ols_modavg2$ConfidenceInterval <- paste0("(", maxBS_ols_modavg2$lower2.5, ", ",
                                               maxBS_ols_modavg2$upper97.5, ")")
write.csv(subset(maxBS_ols_modavg2, (!coefficient == "Intercept")),
          file.path(res.dir, "maxBS_ols_modavg_clean.csv"), row.names = F)

# SAR results
maxBS_sarknn_modavg <- roundNumbers(read.csv(file.path(res.dir, "maxBS_sarknn_modavg.csv")))
maxBS_sarknn_modavg$Geographic.scale <- factor(maxBS_sarknn_modavg$Geographic.scale, levels = geographic_levels)
maxBS_sarknn_modavg$coefficient <- factor(maxBS_sarknn_modavg$coefficient,
                                          levels = max_levels,
                                          labels = max_labels)
maxBS_sarknn_modavg2 <- with(maxBS_sarknn_modavg, maxBS_sarknn_modavg[order(Geographic.scale, Scenario, coefficient),])
maxBS_sarknn_modavg2$ConfidenceInterval <- paste0("(", maxBS_sarknn_modavg2$lower2.5, ", ",
                                                  maxBS_sarknn_modavg2$upper97.5, ")")
write.csv(subset(maxBS_sarknn_modavg2, (!coefficient == "Intercept") & Weighting == "No weighting"),
          file.path(res.dir, "maxBS_sarknn_nw_modavg_clean.csv"),
          row.names = F)
write.csv(subset(maxBS_sarknn_modavg2, (!coefficient == "Intercept") & Weighting == "Distance-weighting"),
          file.path(res.dir, "maxBS_sarknn_dist_modavg_clean.csv"),
          row.names = F)

maxBS_sarsoi_modavg <- roundNumbers(read.csv(file.path(res.dir, "maxBS_sarsoi_modavg.csv")))
maxBS_sarsoi_modavg$Geographic.scale <- factor(maxBS_sarsoi_modavg$Geographic.scale, levels = geographic_levels)
maxBS_sarsoi_modavg$coefficient <- factor(maxBS_sarsoi_modavg$coefficient,
                                          levels = max_levels,
                                          labels = max_labels)
maxBS_sarsoi_modavg2 <- with(maxBS_sarsoi_modavg, maxBS_sarsoi_modavg[order(Geographic.scale, Scenario, coefficient),])
maxBS_sarsoi_modavg2$ConfidenceInterval <- paste0("(", maxBS_sarsoi_modavg2$lower2.5, ", ",
                                                  maxBS_sarsoi_modavg2$upper97.5, ")")
write.csv(subset(maxBS_sarsoi_modavg2, (!coefficient == "Intercept") & Weighting == "No weighting"),
          file.path(res.dir, "maxBS_sarsoi_nw_modavg_clean.csv"),
          row.names = F)
write.csv(subset(maxBS_sarsoi_modavg2, (!coefficient == "Intercept") & Weighting == "Distance-weighting"),
          file.path(res.dir, "maxBS_sarsoi_dist_modavg_clean.csv"),
          row.names = F)

medBS_sarknn_modavg <- roundNumbers(read.csv(file.path(res.dir, "medBS_sarknn_modavg.csv")))
medBS_sarknn_modavg$Geographic.scale <- factor(medBS_sarknn_modavg$Geographic.scale, levels = geographic_levels)
medBS_sarknn_modavg$coefficient <- factor(medBS_sarknn_modavg$coefficient,
                                          levels = med_levels,
                                          labels = med_labels)
medBS_sarknn_modavg2 <- with(medBS_sarknn_modavg, medBS_sarknn_modavg[order(Geographic.scale, Scenario, coefficient),])
medBS_sarknn_modavg2$ConfidenceInterval <- paste0("(", medBS_sarknn_modavg2$lower2.5, ", ",
                                                  medBS_sarknn_modavg2$upper97.5, ")")
write.csv(subset(medBS_sarknn_modavg2, (!coefficient == "Intercept") & Weighting == "No weighting"),
          file.path(res.dir, "medBS_sarknn_nw_modavg_clean.csv"),
          row.names = F)
write.csv(subset(medBS_sarknn_modavg2, (!coefficient == "Intercept") & Weighting == "Distance-weighting"),
          file.path(res.dir, "medBS_sarknn_dist_modavg_clean.csv"),
          row.names = F)

medBS_sarsoi_modavg <- roundNumbers(read.csv(file.path(res.dir, "medBS_sarsoi_modavg.csv")))
medBS_sarsoi_modavg$Geographic.scale <- factor(medBS_sarsoi_modavg$Geographic.scale, levels = geographic_levels)
medBS_sarsoi_modavg$coefficient <- factor(medBS_sarsoi_modavg$coefficient,
                                          levels = med_levels,
                                          labels = med_labels)
medBS_sarsoi_modavg2 <- with(medBS_sarsoi_modavg, medBS_sarsoi_modavg[order(Geographic.scale, Scenario, coefficient),])
medBS_sarsoi_modavg2$ConfidenceInterval <- paste0("(", medBS_sarsoi_modavg2$lower2.5, ", ",
                                                  medBS_sarsoi_modavg2$upper97.5, ")")
write.csv(subset(medBS_sarsoi_modavg2, (!coefficient == "Intercept") & Weighting == "No weighting"),
          file.path(res.dir, "medBS_sarsoi_nw_modavg_clean.csv"),
          row.names = F)
write.csv(subset(medBS_sarsoi_modavg2, (!coefficient == "Intercept") & Weighting == "Distance-weighting"),
          file.path(res.dir, "medBS_sarsoi_dist_modavg_clean.csv"),
          row.names = F)

