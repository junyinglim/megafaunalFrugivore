
createFormula <- function(df){
  lhsList <- as.vector(unique(df$lhs))
  formula <- ""
  for(i in 1:length(lhsList)){
    rhs <- paste(subset(df, lhs == lhsList[i])$rhs, collapse = " + ")
    if(formula == ""){
      formula <- paste0(formula,paste0(lhs[i], " ~ ", rhs))
    } else {
      formula <- paste0(formula, "\n ", paste0(lhs[i], " ~ ", rhs))
    }
  }
  return(formula)
}

fitUpdateSEM <- function(formula, data){
  formula_clean <- gsub(formula, pattern = " ", replacement = "")
  subformula_clean <- strsplit(x = formula_clean, split = "\n")[[1]]
  subformulaPara <- list()
  for(i in 1:length(subformula_clean)){
    subformula_split <- strsplit(subformula_clean[i], split = "~")[[1]]
    subformula_lhs <- subformula_split[1]
    subformula_rhs <- strsplit(subformula_split[2], split = "\\+")[[1]]
    subformulaPara[[i]] <- data.frame("lhs" = subformula_lhs, "rhs" = subformula_rhs)
  }
  formulaParaDF <- do.call("rbind", subformulaPara)
  
  mod <- sem(formula, data  = data) # initialize model
  nvar <- nrow(formulaParaDF)
  lrt <- 1
  while(lrt > 0.05){
    pvals <- round(pnorm((-abs(mod@ParTable$est[1:nvar]) / mod@ParTable$se[1:nvar]))*2, digits = 4)
    criteria <- which.max(pvals) & which(pvals > 0.05)
    varNS <- paste0(formulaParaDF$lhs[criteria]," ~ ",formulaParaDF$rhs[criteria])
    print(paste0("Removing this relationship:", varNS))
    
    formulaParaDF <- formulaParaDF[-criteria,]
    
    print("Fitting new model ... ")
    newmod <- sem(model = createFormula(formulaParaDF), data = data)
    
    lrt <- lavTestLRT(mod, newmod)$'Pr(>Chisq)'[2]
    print(paste0("Likelihood ratio p-value is", lrt))
    
    if(lrt > 0.05){
      mod <- newmod
    }
    
  }
  return(mod)
}
