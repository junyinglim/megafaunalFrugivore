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

require(MuMIn)

# gpa_data <- data.frame(gpa = c(1.97, 2.74, 2.19, 2.60, 2.98, 1.65, 1.89, 2.38, 2.66, 1.96, 3.14, 1.96, 2.20, 3.90, 2.02, 3.61, 3.07, 2.63, 3.11, 3.20), 
#                        sat_math = c(321, 718, 358, 403, 640, 237, 270, 418, 443, 359, 669, 409, 582, 750, 451, 645, 791, 521, 594, 653),
#                        sat_verb = c(247, 436, 578, 447, 563, 342, 472, 356, 327, 385, 664, 518, 364, 632, 435, 704, 341, 483, 665, 606),
#                        hs_math = c(2.3, 3.8, 2.98, 3.58, 3.38, 1.48, 1.67, 3.73, 3.09, 1.54, 3.21, 2.77, 1.47, 3.14, 1.54, 3.50, 3.20, 3.59, 3.42, 3.69),
#                        hs_eng = c(2.63, 3.57, 2.57, 2.21, 3.48, 2.14, 2.64, 2.52, 3.2, 3.46, 3.37, 2.60, 2.90, 3.49, 3.20, 3.74, 2.93, 3.32, 2.70, 3.52))
dredge2 <- function(full_mod){
  # Extract terms from model object
  full_pred_terms <- attr(full_mod$terms, "term.labels")
  lhs_term <- attr(full_mod$terms, "variables")[[2]]
  nterms <- length(full_pred_terms)
  # Create all possible combinations of predictor variables. Remove first row as no predictor variables are included
  mod_list <- expand.grid(rep(list(0:1), nterms))[-1,] 
  names(mod_list) <- full_pred_terms
  
  # Calculate standard deviations of predictors
  pred_sd <- sapply(full_pred_terms, FUN = function(x){ sd(full_mod$model[[x]]) } )
  
  # Create empty lists to hold results
  vif_list <- list()
  mod_AICc <- vector()
  vif_list <- list()
  sd_list <- list()
  pSD_list <- list()
  stdCoef_list <- list()
  coef_list <- list()
  SE_list <- list()
  stdSE_list <- list()

  for(i in 1:nrow(mod_list)){
    rhs_term <- names(mod_list)[mod_list[i,] == 1]
    mod_formula <- as.formula( paste(lhs_term, "~",
                                     paste(rhs_term, collapse = " + ")) )
    mod_obj[[i]] <- lm(mod_formula, data = full_mod$model)
    
    nobs <- nobs(mod_obj[[i]])
    npred <- length(rhs_term)
    npara <- length(coef(mod_obj[[i]]))+1# includes intercept and standard dev. in AIC calculations
    
    # Calculate model AICc
    mod_AICc[i] <- AIC(mod_obj[[i]]) + ((2*npara ^2 + 2*npara) / (nobs-npara-1))
    
    sd_list[[i]] <- pred_sd[rhs_term] # sd of the variable
    coef_list[[i]] <- coef( mod_obj[[i]] )[rhs_term]
    
    if(length(rhs_term) > 1){
      vif_list[[i]] <- vif(mod_obj[[i]])
      pSD_list[[i]] <- sd_list[[i]] * vif_list[[i]]^-0.5 * ((nobs - 1) / (nobs - npred))^0.5
      stdCoef_list[[i]] <- coef( mod_obj[[i]] )[rhs_term] * pSD_list[[i]]
      SE_list[[i]] <- summary(mod_obj[[i]])$coefficients[,2][rhs_term]
      stdSE_list[[i]] <- (SE_list[[i]]^2 * pSD_list[[i]]^2)^0.5
    } else {
      vif_list[[i]] <- setNames(1, rhs_term)
      pSD_list[[i]] <- NA # Partial SD undefined for single variable models
      stdCoef_list[[i]] <- coef( mod_obj[[i]] )[rhs_term] * sd_list[[i]]
      SE_list[[i]] <- summary(mod_obj[[i]])$coefficients[,2][rhs_term]
      stdSE_list[[i]] <- (SE_list[[i]]^2 * sd_list[[i]]^2)^0.5
    }
  }
  mod_dAICc <- mod_AICc - min(mod_AICc)
  mod_wt <- exp(-0.5 * mod_dAICc) / sum(exp(-0.5 * mod_dAICc))
  list("Models" = mod_list,
              "Model obj" = mod_obj,
              "AICc" = mod_AICc,
              "dAICc" = mod_dAICc,
              "Weights" = mod_wt,
              "SD" = sd_list,
              "VIF" = vif_list,
              "pSD" = pSD_list,
              "Coef" = coef_list,
              "stdCoef" = stdCoef_list, 
              "SE" = SE_list,
              "stdSE" = stdSE_list)
}


model.avg2 <- function(model_list){
  # Generate model averaged coefficients
  model_list <- full_moddr
  predVars <- names(model_list$Models)
  
  avgStdCoef <- vector()
  avgUnstdCoef <- vector()
  avgStdSE <- vector()
  avgUnstdCoef <- vector()
  var_impt <- vector()
  for(i in 1:length(predVars)){
    
    modInd <- which(model_list$Models[predVars[i]] == 1) # models with target pred. variable
    weights <- model_list$Weights[modInd]
    weight_sum <- sum(weights)
    
    var_impt[i] <- weight_sum
    
    StdCoefs <- unlist(lapply(model_list$stdCoef[modInd], function(x) { x[predVars[i]]} ))
    avgStdCoef[i] <- sum(StdCoefs * (weights/weight_sum)) 
    
    UnstdCoefs <- unlist(lapply(model_list$Coef[modInd], function(x) { x[predVars[i]]} ))
    avgUnstdCoef[i] <- sum(UnstdCoefs * (weights/weight_sum)) 
    
    # See Burnham & Anderson 2002, p180 for SE calculation
    StdSEs <- unlist(lapply(model_list$stdSE[modInd], function(x) { x[predVars[i]]} ))
    avgStdSE[i] <- sum(( (StdCoefs - avgStdCoef[i])^2 + (StdSEs)^2)^0.5 * (weights / weight_sum))
    
    UnstdSEs <- unlist(lapply(model_list$SE[modInd], function(x) { x[predVars[i]]} ))
    avgUnstdSE[i] <- sum(( (UnstdCoefs - avgUnstdCoef[i])^2 + (UnstdSEs)^2)^0.5 * (weights / weight_sum))
    
  }
  data.frame("Variable" = predVars,
             avgStdCoef = avgStdCoef,
             avgStdSE = avgStdSE,
             avgUnstdCoef = avgUnstdCoef,
             avgUnstdSE = avgUnstdSE,
             var_impt = var_impt)
}

full_mod <- lm(gpa ~ sat_math + sat_verb + hs_math + hs_eng, data = gpa_data, na.action = "na.fail")
full_moddr <- dredge2(full_mod)
full_modavg <- model.avg2(full_moddr)

full_moddr2 <- dredge(full_mod)
full_modavg2 <- model.avg(full_moddr2)
