require(ggplot2)

summarizeLM <- function(mod, scale = "global", scenario = "current"){
  
  mod_summary <- summary(mod)
  mod_coeff <- as.data.frame(mod_summary$coefficients)
  n_coeff <- nrow(mod_coeff)
  varnames <- gsub(rownames(mod_coeff), pattern = "log\\(", replacement = "log_")
  varnames <- gsub(varnames, pattern = "scale\\(|\\)|\\(", replacement = "")
  adj.r.squared <- mod_summary$adj.r.squared
  mod_coeff_est <- paste(round(mod_coeff$Estimate, 3), "±", round(mod_coeff$`Std. Error`, 3))
  mod_coeff_p <- mod_coeff$`Pr(>|t|)`
  mod_coeff_p_sym <- ifelse(mod_coeff_p < 0.001, "***", ifelse(mod_coeff_p < 0.01, "**", ifelse(mod_coeff_p < 0.05, "*", "n.s.") ))
  mod_aic <- round(AIC(mod), 3)
  data.frame("Geographic Scale" = c(scale, rep("", n_coeff-1)),
             "Scenario" = c(scenario, rep("", n_coeff-1)),
             "Variable" = varnames,
             "Coefficient" = mod_coeff_est,
             "P-value" = mod_coeff_p_sym,
             "Adjusted R sq." = c(round(adj.r.squared, digits = 3), rep("", n_coeff-1)),
             "AIC" = c(mod_aic, rep("", n_coeff-1)))
}



summarizeRelImportance <- function(x, plotIntercept = FALSE){
  # Plots relative
  # x = "averaging" class from MuMIn package
  #x = fullmod_avg
  summaryStats <- data.frame(coefficient = colnames(x$coefficients), fullAvgCoef = x$coefficients[1,], condAvgCoef = x$coefficients[2,])
  
  relimportStats <- data.frame(importance = x$importance, coefficient = names(x$importance))
  confintStats <- data.frame(confint(x))
  names(confintStats) <- c("lower2.5", "upper97.5")
  summaryStats <- cbind(summaryStats, confintStats)
  if(!plotIntercept){
    summaryStats <- summaryStats[-1,]    
  }
  summaryStats <- merge(summaryStats, relimportStats, by = "coefficient")
  return(summaryStats)
}

plotRelImportance <- function(x){
  # Plots basic diagnostic output of the summarizeRelImportance function
  ggplot(data = x) + geom_point(aes(y = fullAvgCoef, x= coefficient, size = importance)) + geom_segment(aes(y = lower2.5, yend = upper97.5, x = coefficient, xend = coefficient)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
}


# x$coefficients
# 
# names(summaryStats)
# ?geom_segment
# str(fullmod_avg, max.level = 1)
# fullmod_avg$coefficients # full average (full) and conditional average (subset)
# 
# fullmod_avg$importance
# vif(fullmod)
# ?vif
# fullmod_avg$importance
# fullmod_avg$formula
# ?model.avg
