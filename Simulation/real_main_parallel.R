rm(list=ls())

setwd("~/Documents/GitHub/Impact_Non_Informative_Censoring_MHR/Simulation")


source("gen_surv_data.R")
source("beta_bissection.R")
source("marginal_beta.R")
source("censor_data.R")
source("cox_reg.R")
library(survival)
library(doParallel)
library(doRNG)
library(rngtools)

seed = 2355
n <- 1000 #Sample size
n_back = n
rep <- 100 #Number of iterations


rng <- RNGseq(rep, seed)
cl <- makeCluster(1) #28 usually, 1 for local
registerDoParallel(cl)
registerDoRNG(seed)


censoring_list = c(0.6, 0.7, 0.8, 0.9)
censoring_list = c(0.8)
target_list = c(2)
counterfactual_list = c(TRUE)
type_censoring_list = c('uniform')
#counterfactual = TRUE
#censoring_p = 0.7
#target_beta = 2
#type_censoring = 'uniform'

for(counterfactual in counterfactual_list){
for(censoring_p in censoring_list){
  for(target_beta in target_list){
    for(type_censoring in type_censoring_list){
        set.seed(seed) 
      n = n_back      
counterfactual_text = 'notcounterfactual'
if(counterfactual == TRUE){
  n <- n_back/2
  counterfactual_text = 'counterfactual'
}
print(paste0(counterfactual_text,
             '_n', n_back, 'rep', rep,
             '_target', target_beta,
             '_censoring', censoring_p,
             '_censTYPE',type_censoring,
             '_.RData'))
#List to store all generated data
CHR_est_uncens_res <- CHR_est_res <- MHR_est_unweighted_uncens_res <- IPTW_uncens_res <- IPTW_res<- IPTPEW_res <- matrix(nrow = rep, ncol = 5)
IPTWPE2_res = IPTPEW_res
MATCH_res = IPTPEW_res
MATCH_PE_res = IPTPEW_res
MATCH_PE2_res = IPTPEW_res
NAIVE_res = IPTPEW_res
IPTW_real_pc = IPTPEW_res
MATCH_real_pc = IPTPEW_res
IPTW2_real_pc = IPTPEW_res
MATCH2_real_pc = IPTPEW_res

chr <- numeric(rep)

# This is weird, but sometimes not calling the functions again here gives an error.
source("gen_surv_data.R")
source("beta_bissection.R")
source("marginal_beta.R")
source("censor_data.R")
source("cox_reg.R")

sim_results = list(CHR_est_uncens_res = matrix(nrow = rep, ncol = 5),
                   CHR_est_res= matrix(nrow = rep, ncol = 5),
                   MHR_est_unweighted_uncens_res = matrix(nrow = rep, ncol = 5),
                   IPTW_uncens_res = matrix(nrow = rep, ncol = 5),
                   IPTW_res = matrix(nrow = rep, ncol = 5),
                   IPTPEW_res = matrix(nrow = rep, ncol = 5),
                   IPTPEW2_res = matrix(nrow = rep, ncol = 5),
                   IPTW_real_pc = matrix(nrow = rep, ncol = 5),
                   IPTW_real_pc2 = matrix(nrow = rep, ncol = 5),
                   MATCH_res  = matrix(nrow = rep, ncol = 5),
                   MATCH_PE_res  = matrix(nrow = rep, ncol = 5),
                   MATCH_PE2_res  = matrix(nrow = rep, ncol = 5),
                   MATCH_real_pc = matrix(nrow = rep, ncol = 5),
                   MATCH_real_pc2 = matrix(nrow = rep, ncol = 5),
                   NAIVE_res = matrix(nrow = rep, ncol = 5),
                   mhr = rep(NA,rep),
                   chr = rep(NA,rep))
reslist_n_rep = foreach(i = 1:rep, .packages = c('survival', 'MatchIt', 'Matching'),.export = 'rng')%dopar%{
  #print('I started the loop!')
  
  source("gen_surv_data_ATE_220823.R")
  source("beta_bissection_ATE_JH.R")
  source("marginal_beta_ATE_JH.R")
  source("censor_data.R")
  source("cox_reg_220823.R")
  ########################
  eta = 2
  lambda = 0.00002
  datalist <- gen_surv_data_ATE(n = n, counterfactual = counterfactual, lambda = lambda, eta = eta, scenario = 'A',m = 1, k = 1, 
                                     marginal_hazard_ratio = target_beta, censoring_p = censoring_p, type_censoring = type_censoring)
  #print('I generated data!')
  
  mhr <- datalist$marginal_hazard_ratio
  chr <- datalist$conditional_hazard_ratio
  data <- datalist$data
  
  ###########################
  #Estimate true prob of censoring 
  ###########################
  
  
  theta = datalist$theta
  lambda_i = datalist$lambda_i
  
  
  if(type_censoring == 'uniform'){
    true_cens_p = (lambda_i/(eta*theta))*pgamma((theta/lambda_i)^eta,1/eta)*gamma(1/eta)
  }
  if(type_censoring == 'weibull'){
    true_cens_p = 1/(1 + ((theta/lambda_i)^eta))
  }
  true_event_p = 1-true_cens_p
  mean(data$event)
  mean(true_event_p)
  data$true_cens_p = true_cens_p
  #print('I calculated real event probs!')
  
  #ggplot2::qplot(data$true_cens_p,data$Y_true)
  cor(true_event_p,data$Y_true)
  cor(true_event_p,data$Treatment)
  cor(true_event_p,data$event)

  
  
  
  ###########################
  #Estimate propensity score 
  ###########################
  #Logistic regression and correct model specification
  ps_est <- glm(Treatment ~ X1 + X2 + X3 + X4 + X5 + X6 + X7+ X8 +
                  X9 + X10, data = data, family = 'binomial')
  data$ps_est <- ps_est$fitted.values
 
  ##########################
  #Calculate ATE ps-weights
  ##########################
  data$pt_weights <- ifelse(data$Treatment == 1, 1/data$ps_est, 1/(1-data$ps_est))
  
  ###########################
  #Store uncensored data 
  ###########################
  data_uncens <- subset(data,  select = -Y)
  data_uncens <- cbind(data_uncens[, which(colnames(data_uncens) == 'Y_true')], data_uncens[, -which(colnames(data_uncens) == 'Y_true')])
  colnames(data_uncens)[1] <- c("Y")
  data_uncens$event <- 1
  
  ##############################################
  #UNCENSORED: Calculate CHR estimate
  ##############################################
  time <- data_uncens$Y
  event <- data_uncens$event
  covariates <- data_uncens[, 3:13]
  CHR_est_uncens <- cox_reg(time, event, covariates, weights = NULL)
  
  
  ##############################################
  #UNCENSORED: Calculate unweighted MHR estimate
  ##############################################
  covariates <- data.frame(Treatment = data_uncens[, 3])
  MHR_est_unweighted_uncens <- cox_reg(time, event, covariates, weights = NULL)
  
  
  
  ##################################
  #UNCENSORED: Calculate common IPTW
  ##################################
  weights <- data_uncens$pt_weights
  IPTW_uncens <- cox_reg(time, event, covariates, weights)
  
  ##############################################
  #CENSORED: Calculate CHR estimate
  ##############################################
  time <- data$Y
  event <- data$event
  covariates <- data[, 3:13]
  CHR_est <- cox_reg(time, event, covariates, weights = NULL)
  
  ##########################
  #CENSORED: Calculate naive IPTW
  ##########################
  covariates <- data.frame(Treatment = data[, 3])
  weights <- data$pt_weights
  IPTW_naive <- cox_reg(time, event, covariates, weights = rep(1,nrow(data)))
  
  ##########################
  #CENSORED: Calculate common IPTW
  ##########################
  covariates <- data.frame(Treatment = data[, 3])
  weights <- data$pt_weights
  IPTW <- cox_reg(time, event, covariates, weights)

  ###################################
  #Estimate probability of censoring 
  ###################################
  #Unweighted logistic regression including both treatment and covariates in the model
  data$censored <- ifelse(data$event == 1, 0, 1)
  pc_est <-glm(censored ~ Treatment + X1 + X2 + X3 + X4 + X5 + X6 + X7+ X8 +
                 X9 + X10, data = data, family = 'binomial')
  data$pc_est <- pc_est$fitted.values
 
  ####################################################
  #CENSORED: Calculate IPTW with unweighted PE weights
  ####################################################
  weights <- data$pt_weights * (1-data$pc_est)
  IPTPEW <- cox_reg(time, event, covariates, weights)
  
  weights_pc = ifelse(data$event == 1, data$pc_est, (1-data$pc_est))
  weights <- data$pt_weights * weights_pc
  
  ## Remember to comment out
  ##weights <- data$pt_weights * (data$pc_est)
  
  IPTPEW2 <- cox_reg(time, event, covariates, weights)
  
  
  ####################################################
  #CENSORED: Calculate IPTW with TRUE PE weights
  ####################################################
  weights <- data$pt_weights * (1-data$true_cens_p)
  IPTW_real_pc_line <- cox_reg(time, event, covariates, weights)
  
  weights_pc = ifelse(data$event == 1, data$true_cens_p, (1-data$true_cens_p))
  weights <- data$pt_weights * weights_pc
  IPTW_real_pc_line2 <- cox_reg(time, event, covariates, weights)
  
  
  ##########################
  #CENSORED: Calculate common full matching
  ##########################
  
  match.object_ATE = suppressWarnings(MatchIt::matchit(Treatment~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10, data = data, method = 'full', caliper = 0.2, 
                                              estimand = 'ATE', distance = logit(data$ps_est)))
  matched.data_ATE = MatchIt::match.data(match.object_ATE)
  matched.data_ATE$set = matched.data_ATE$subclass
  matched_rows = as.numeric(row.names(matched.data_ATE))
  
  fmla <- paste("Surv(time = Y, event = event) ~ Treatment + cluster(set)")
  
  result_match_ATE = coxph(as.formula(fmla), data = matched.data_ATE, weights = matched.data_ATE$weights,
                           control = coxph.control(timefix = FALSE))
  beta = result_match_ATE$coefficients[1]
  rob_se <- sqrt(result_match_ATE$var[1, 1])
  MATCH = list(beta = beta , exp_beta = exp(beta),
       ciL = exp(beta - 1.96 * rob_se), ciU = exp(beta + 1.96 * rob_se))
  
  ##########################
  #CENSORED: Calculate  PE full matching
  ##########################
  
  result_match_PE = coxph(as.formula(fmla), data = matched.data_ATE, 
                          weights = matched.data_ATE$weights*(1-data[matched_rows,]$pc_est))
  beta = result_match_PE$coefficients[1]
  rob_se <- sqrt(result_match_PE$var[1, 1])
  MATCH_PE = list(beta = beta , exp_beta = exp(beta),
               ciL = exp(beta - 1.96 * rob_se), ciU = exp(beta + 1.96 * rob_se))
  
  weights_pc = ifelse(data$event == 0, data$pc_est, (1-data$pc_est))
  result_match_PE2 = coxph(as.formula(fmla), data = matched.data_ATE, 
                          weights = matched.data_ATE$weights*(weights_pc[matched_rows]))
  beta = result_match_PE2$coefficients[1]
  rob_se <- sqrt(result_match_PE2$var[1, 1])
  MATCH_PE2 = list(beta = beta , exp_beta = exp(beta),
                  ciL = exp(beta - 1.96 * rob_se), ciU = exp(beta + 1.96 * rob_se))
  
  ##########################
  #CENSORED: Calculate true PE full matching
  ##########################
  

  result_match_PE = coxph(as.formula(fmla), data = matched.data_ATE, 
                          weights = matched.data_ATE$weights*(1-data[matched_rows,]$true_cens_p))
  beta = result_match_PE$coefficients[1]
  rob_se <- sqrt(result_match_PE$var[1, 1])
  MATCH_real_pc_line = list(beta = beta , exp_beta = exp(beta),
                  ciL = exp(beta - 1.96 * rob_se), ciU = exp(beta + 1.96 * rob_se))
  
  weights_pc = ifelse(data$event == 0, data$true_cens_p, (1-data$true_cens_p))
  result_match_PE2 = coxph(as.formula(fmla), data = matched.data_ATE, 
                           weights = matched.data_ATE$weights*(weights_pc[matched_rows]))
  beta = result_match_PE2$coefficients[1]
  rob_se <- sqrt(result_match_PE2$var[1, 1])
  MATCH_real_pc_line2 = list(beta = beta , exp_beta = exp(beta),
                   ciL = exp(beta - 1.96 * rob_se), ciU = exp(beta + 1.96 * rob_se))
  print('I did all calculations!')
  
  ##########################
  #Saving results
  ##########################
  
  CHR_est_uncens_res <- c(CHR_est_uncens$exp_beta, CHR_est_uncens$exp_beta- chr, CHR_est_uncens$ciL,
                           CHR_est_uncens$ciU, (CHR_est_uncens$exp_beta- chr)/chr)
  CHR_est_res <- c(CHR_est$exp_beta, CHR_est$exp_beta- chr, CHR_est$ciL,
                           CHR_est$ciU, (CHR_est$exp_beta- chr)/chr)
  MHR_est_unweighted_uncens_res <- c(MHR_est_unweighted_uncens$exp_beta, MHR_est_unweighted_uncens$exp_beta- mhr, MHR_est_unweighted_uncens$ciL,
                                      MHR_est_unweighted_uncens$ciU, (MHR_est_unweighted_uncens$exp_beta- mhr)/mhr)
  IPTW_uncens_res <- c(IPTW_uncens$exp_beta, IPTW_uncens$exp_beta- mhr, IPTW_uncens$ciL,
                        IPTW_uncens$ciU, (IPTW_uncens$exp_beta- mhr)/mhr)
  IPTW_res <- c(IPTW$exp_beta, IPTW$exp_beta- mhr, IPTW$ciL,
                 IPTW$ciU, (IPTW$exp_beta- mhr)/mhr)
  IPTPEW_res <- c(IPTPEW$exp_beta, IPTPEW$exp_beta- mhr, IPTPEW$ciL,
                   IPTPEW$ciU, (IPTPEW$exp_beta- mhr)/mhr)
  IPTPEW2_res <- c(IPTPEW2$exp_beta, IPTPEW2$exp_beta- mhr, IPTPEW2$ciL,
                  IPTPEW2$ciU, (IPTPEW2$exp_beta- mhr)/mhr)
  IPTW_real_pc <- c(IPTW_real_pc_line$exp_beta, IPTW_real_pc_line$exp_beta- mhr, IPTW_real_pc_line$ciL,
                    IPTW_real_pc_line$ciU, (IPTW_real_pc_line$exp_beta- mhr)/mhr)
  IPTW_real_pc2 <- c(IPTW_real_pc_line2$exp_beta, IPTW_real_pc_line2$exp_beta- mhr, IPTW_real_pc_line2$ciL,
                    IPTW_real_pc_line2$ciU, (IPTW_real_pc_line2$exp_beta- mhr)/mhr)
  MATCH_res <- c(MATCH$exp_beta, MATCH$exp_beta- mhr, MATCH$ciL,
                      MATCH$ciU, (MATCH$exp_beta- mhr)/mhr)
  MATCH_PE_res <- c(MATCH_PE$exp_beta, MATCH_PE$exp_beta- mhr, MATCH_PE$ciL,
                         MATCH_PE$ciU, (MATCH_PE$exp_beta- mhr)/mhr)
  MATCH_PE2_res <- c(MATCH_PE2$exp_beta, MATCH_PE2$exp_beta- mhr, MATCH_PE2$ciL,
                    MATCH_PE2$ciU, (MATCH_PE2$exp_beta- mhr)/mhr)
  MATCH_real_pc <- c(MATCH_real_pc_line$exp_beta, MATCH_real_pc_line$exp_beta- mhr, MATCH_real_pc_line$ciL,
                    MATCH_real_pc_line$ciU, (MATCH_real_pc_line$exp_beta- mhr)/mhr)
  MATCH_real_pc2 <- c(MATCH_real_pc_line2$exp_beta, MATCH_real_pc_line2$exp_beta- mhr, MATCH_real_pc_line2$ciL,
                     MATCH_real_pc_line2$ciU, (MATCH_real_pc_line2$exp_beta- mhr)/mhr)
  
  NAIVE_res <- c(IPTW_naive$exp_beta, IPTW_naive$exp_beta- mhr, IPTW_naive$ciL,
                      IPTW_naive$ciU, (IPTW_naive$exp_beta- mhr)/mhr)
  
  sim_results_line1 <- list(CHR_est_uncens_res = CHR_est_uncens_res, CHR_est_res = CHR_est_res, MHR_est_unweighted_uncens_res = MHR_est_unweighted_uncens_res,
                      IPTW_uncens_res = IPTW_uncens_res, IPTW_res = IPTW_res, IPTPEW_res = IPTPEW_res,
                      IPTPEW2_res = IPTPEW2_res,
                      IPTW_real_pc = IPTW_real_pc,IPTW_real_pc2 = IPTW_real_pc2,
                      MATCH_res = MATCH_res, 
                      MATCH_PE_res = MATCH_PE_res, MATCH_PE2_res = MATCH_PE2_res, 
                      MATCH_real_pc = MATCH_real_pc, MATCH_real_pc2 = MATCH_real_pc2,
                      NAIVE_res = NAIVE_res,
                      chr = chr, mhr = mhr, censoring_p = censoring_p, n = n, counterfactual = counterfactual)
  sim_results_line1
  #print('I finished a loop!')
  
}

      for(j in 1:rep){
        
        sim_results_line = reslist_n_rep[[j]]
        sim_results$CHR_est_uncens_res[j,] = sim_results_line$CHR_est_uncens_res
        sim_results$CHR_est_res[j,] = sim_results_line$CHR_est_res
        sim_results$MHR_est_unweighted_uncens_res[j,] = sim_results_line$MHR_est_unweighted_uncens_res
        sim_results$IPTW_uncens_res[j,] = sim_results_line$IPTW_uncens_res
        sim_results$IPTW_res[j,] = sim_results_line$IPTW_res
        sim_results$IPTPEW_res[j,] = sim_results_line$IPTPEW_res
        sim_results$IPTPEW2_res[j,] = sim_results_line$IPTPEW2_res
        sim_results$IPTW_real_pc[j,] = sim_results_line$IPTW_real_pc
        sim_results$IPTW_real_pc2[j,] = sim_results_line$IPTW_real_pc2
        sim_results$MATCH_res[j,] = sim_results_line$MATCH_res
        sim_results$MATCH_PE_res[j,] = sim_results_line$MATCH_PE_res
        sim_results$MATCH_PE2_res[j,] = sim_results_line$MATCH_PE2_res
        sim_results$MATCH_real_pc[j,] = sim_results_line$MATCH_real_pc
        sim_results$MATCH_real_pc2[j,] = sim_results_line$MATCH_real_pc2
        sim_results$NAIVE_res[j,] = sim_results_line$NAIVE_res
        sim_results$chr[j] = sim_results_line$chr
        sim_results$mhr[j] = sim_results_line$mhr
        
        
        
      }
      
      bias <- c(mean(sim_results$CHR_est_uncens_res[, 2]), mean(sim_results$CHR_est_res[, 2]), mean(sim_results$MHR_est_unweighted_uncens_res[, 2]), 
                mean(sim_results$IPTW_uncens_res[, 2]), mean(sim_results$IPTW_res[, 2]), 
                mean(sim_results$IPTPEW_res[, 2]), mean(sim_results$IPTPEW2_res[, 2]),
                mean(sim_results$IPTW_real_pc[, 2]), mean(sim_results$IPTW_real_pc2[, 2]), 
                mean(sim_results$MATCH_res[, 2]), mean(sim_results$MATCH_PE_res[, 2]),
                mean(sim_results$MATCH_PE2_res[, 2]),
                mean(sim_results$MATCH_real_pc[, 2]), mean(sim_results$MATCH_real_pc2[, 2]),
                mean(sim_results$NAIVE_res[, 2]))
      
      sd <- c(sd(sim_results$CHR_est_uncens_res[, 1]), sd(sim_results$CHR_est_res[, 1]), sd(sim_results$MHR_est_unweighted_uncens_res[, 1]), 
              sd(sim_results$IPTW_uncens_res[, 1]), sd(sim_results$IPTW_res[, 1]), 
              sd(sim_results$IPTPEW_res[, 1]),sd(sim_results$IPTPEW2_res[, 1]),
              sd(sim_results$IPTW_real_pc[, 1]), sd(sim_results$IPTW_real_pc2[, 1]), 
              sd(sim_results$MATCH_res[, 1]), 
              sd(sim_results$MATCH_PE_res[, 1]),sd(sim_results$MATCH_PE2_res[, 1]),
              sd(sim_results$MATCH_real_pc[, 1]), sd(sim_results$MATCH_real_pc2[, 1]), 
              sd(sim_results$NAIVE_res[, 1]))
      
      rmse <- sqrt(bias^2 + sd^2)
      
      rel_bias <- c(mean(sim_results$CHR_est_uncens_res[, 5]), mean(sim_results$CHR_est_res[, 5]), mean(sim_results$MHR_est_unweighted_uncens_res[, 5]), 
                    mean(sim_results$IPTW_uncens_res[, 5]), mean(sim_results$IPTW_res[, 5]), 
                    mean(sim_results$IPTPEW_res[, 5]),mean(sim_results$IPTPEW2_res[, 5]),
                    mean(sim_results$IPTW_real_pc[, 5]), mean(sim_results$IPTW_real_pc2[, 5]), 
                    mean(sim_results$MATCH_res[, 5]), 
                    mean(sim_results$MATCH_PE_res[, 5]),mean(sim_results$MATCH_PE2_res[, 5]),
                    mean(sim_results$MATCH_real_pc[, 5]), mean(sim_results$MATCH_real_pc2[, 5]),
                    mean(sim_results$NAIVE_res[, 5]))
      
      ciL <- c(mean(sim_results$CHR_est_uncens_res[, 3]), mean(sim_results$CHR_est_res[, 3]), mean(sim_results$MHR_est_unweighted_uncens_res[, 3]), 
               mean(sim_results$IPTW_uncens_res[, 3]), mean(sim_results$IPTW_res[, 3]), 
               mean(sim_results$IPTPEW_res[, 3]),mean(sim_results$IPTPEW2_res[, 3]),
               mean(sim_results$IPTW_real_pc[, 3]), mean(sim_results$IPTW_real_pc2[, 3]), 
               mean(sim_results$MATCH_res[, 3]), 
               mean(sim_results$MATCH_PE_res[, 3]),mean(sim_results$MATCH_PE2_res[, 3]),
               mean(sim_results$MATCH_real_pc[, 3]),  mean(sim_results$MATCH_real_pc2[, 3]), 
               mean(sim_results$NAIVE_res[, 3]))
      
      
      ciU <- c(mean(sim_results$CHR_est_uncens_res[, 4]), mean(sim_results$CHR_est_res[, 4]), mean(sim_results$MHR_est_unweighted_uncens_res[, 4]), 
               mean(sim_results$IPTW_uncens_res[, 4]), mean(sim_results$IPTW_res[, 4]), 
               mean(sim_results$IPTPEW_res[, 4]),mean(sim_results$IPTPEW2_res[, 4]),
               mean(sim_results$IPTW_real_pc[, 4]),  mean(sim_results$IPTW_real_pc2[, 4]),
               mean(sim_results$MATCH_res[, 4]), 
               mean(sim_results$MATCH_PE_res[, 4]),mean(sim_results$MATCH_PE2_res[, 4]),
               mean(sim_results$MATCH_real_pc[, 4]), mean(sim_results$MATCH_real_pc2[, 4]), 
               mean(sim_results$NAIVE_res[, 4]))
      
      
      cp_CHR_est_uncens_res <- cp_CHR_est_res <- numeric(length(sim_results$chr))
      
      for(j in 1:rep){
        cp_CHR_est_uncens_res[j] <- ifelse(findInterval(sim_results$chr[j],c(sim_results$CHR_est_uncens_res[j,3], sim_results$CHR_est_uncens_res[j,4]), rightmost.closed = TRUE)==1, 1, 0)
        cp_CHR_est_res[j] <-  ifelse(findInterval(sim_results$chr[j],c(sim_results$CHR_est_res[j,3], sim_results$CHR_est_res[j,4]), rightmost.closed = TRUE)==1, 1, 0)
      }
      
      cp <- c(mean(cp_CHR_est_uncens_res), mean(cp_CHR_est_res), 
              mean(apply(sim_results$MHR_est_unweighted_uncens_res, 1, function(x) ifelse(findInterval(sim_results$mhr,c(x[3], x[4]), rightmost.closed = TRUE)==1, 1, 0))),
              mean(apply(sim_results$IPTW_uncens_res, 1, function(x) ifelse(findInterval(sim_results$mhr,c(x[3], x[4]), rightmost.closed = TRUE)==1, 1, 0))),
              mean(apply(sim_results$IPTW_res, 1, function(x) ifelse(findInterval(sim_results$mhr,c(x[3], x[4]), rightmost.closed = TRUE)==1, 1, 0))),
              mean(apply(sim_results$IPTPEW_res, 1, function(x) ifelse(findInterval(sim_results$mhr,c(x[3], x[4]), rightmost.closed = TRUE)==1, 1, 0))),
              mean(apply(sim_results$IPTPEW2_res, 1, function(x) ifelse(findInterval(sim_results$mhr,c(x[3], x[4]), rightmost.closed = TRUE)==1, 1, 0))),
              mean(apply(sim_results$IPTW_real_pc, 1, function(x) ifelse(findInterval(sim_results$mhr,c(x[3], x[4]), rightmost.closed = TRUE)==1, 1, 0))),
              mean(apply(sim_results$IPTW_real_pc2, 1, function(x) ifelse(findInterval(sim_results$mhr,c(x[3], x[4]), rightmost.closed = TRUE)==1, 1, 0))),
              mean(apply(sim_results$MATCH_res, 1, function(x) ifelse(findInterval(sim_results$mhr,c(x[3], x[4]), rightmost.closed = TRUE)==1, 1, 0))),
              mean(apply(sim_results$MATCH_PE_res, 1, function(x) ifelse(findInterval(sim_results$mhr,c(x[3], x[4]), rightmost.closed = TRUE)==1, 1, 0))),
              mean(apply(sim_results$MATCH_PE2_res, 1, function(x) ifelse(findInterval(sim_results$mhr,c(x[3], x[4]), rightmost.closed = TRUE)==1, 1, 0))),
              mean(apply(sim_results$MATCH_real_pc, 1, function(x) ifelse(findInterval(sim_results$mhr,c(x[3], x[4]), rightmost.closed = TRUE)==1, 1, 0))),
              mean(apply(sim_results$MATCH_real_pc2, 1, function(x) ifelse(findInterval(sim_results$mhr,c(x[3], x[4]), rightmost.closed = TRUE)==1, 1, 0))),
              mean(apply(sim_results$NAIVE_res, 1, function(x) ifelse(findInterval(sim_results$mhr,c(x[3], x[4]), rightmost.closed = TRUE)==1, 1, 0)))
      )
      
      
      res_table <- cbind(bias, sd, rmse, rel_bias, ciL, ciU, cp)
      colnames(res_table) <- c("Bias", "SD", "RMSE", "Rel.bias", "CIL", "CIU", "CP")
      rownames(res_table) <- c("CHR_est_uncens", "CHR_est", "MHR_est_unweighted_uncens",
                               "IPTW_uncens", "IPTW",  "IPTPEW",  "IPTPEW2", 
                               "IPTW_true_PE","IPTW_true_PE2",
                               "MATCH", "MATCH_PE","MATCH_PE2",
                               "MATCH_true_PE","MATCH_true_PE2", 'NAIVE')
      
      
      
      final_list = list(sim_results = sim_results, res_table = res_table)
      save(final_list, file = paste0('Results/',counterfactual_text,
                                     '_n', n_back, 'rep', rep,
                                     '_target', target_beta,
                                     '_censoring', censoring_p,
                                     '_censTYPE',type_censoring,
                                     '_.RData'))
    }
  }
}
}
round(res_table,3)

