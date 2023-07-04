
cox_reg = function(time, event, covariates, weights = NULL){

  if(!is.null(weights)){
    data <- cbind(time, event, covariates, weights)
    data$set <- 1:nrow(data)
    fmla <- as.formula(paste("Surv(time, event) ~ ", paste(names(covariates), collapse= "+"), "+ cluster(set)"))
    coxph_reg <- coxph(fmla, data = data, weights = weights)
    
  } else {
    data <- cbind(time, event, covariates)
    fmla <- as.formula(paste("Surv(time, event) ~ ", paste(names(covariates), collapse= "+")))
    coxph_reg <- coxph(fmla, data = data)
    }
  
  
  beta <- coxph_reg$coefficients[1]
  exp_beta <- exp(coxph_reg$coefficients[1])
  rob_se <- sqrt(coxph_reg$var[1, 1])
  names(beta) <- names(exp_beta) <- names(rob_se) <- NULL
  ciL <- exp(beta - 1.96 * rob_se)
  ciU <- exp(beta + 1.96 * rob_se)
  wald <- (beta / rob_se) ^ 2
  pval <- pchisq(wald, df = 1, lower.tail = FALSE)
  
  
  return(list(beta = beta , exp_beta = exp_beta,
              ciL = ciL, ciU = ciU))
}

logit <- function(x){
  res <- log(x/ (1 - x))
  return(res)
}


