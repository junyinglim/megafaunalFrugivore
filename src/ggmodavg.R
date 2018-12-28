require(ggplot2)

plotRelImportance <- function(x, plotIntercept = FALSE){
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
  
  #ggplot(data = summaryStats) + geom_point(aes(y = fullAvgCoef, x= coefficient, size = importance)) + geom_segment(aes(y = lower2.5, yend = upper97.5, x = coefficient, xend = coefficient)) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
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
