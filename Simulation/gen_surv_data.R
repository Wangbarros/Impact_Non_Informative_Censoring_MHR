#' @title Generate survival data with a specified marginal hazard ratio and
#'   censoring rate
#'
#' @description Function used to generate survival data, mostly for testing and
#'   example purposes. The data is generated with a specified marginal hazard
#'   ratio using a bissection method according to aa(2122) and covariate
#'   scenarios and transformations according to bb(21212). The data generated
#'   can also be censored with a uniform or weibull distribution with a
#'   specified censoring rate, this process is an implementation of ccc(21212).
#'   The time to event outcome follows
#'   \deqn{Y_i = \left( \frac{-\text{log}(u_i)}{\lambda e^{LP_i}} \right)^{1/\eta}}
#'   Where LP is the linear predictor of a subject, with the effects of the covariates and the treament on survival time
#'   and the treatment assignment is a logit transformation on the covariates effects on treatment.
#'   For more details see ddd(12121)
#'
#' @param n A numeric value of how many observations/subjects the dataset should have.
#' @param lambda A numeric value representing the scale parameter.
#' @param eta A numeric value representing the shape parameter.
#' @param m A numeric value representing the coeffecient of covariate effect on survival time.
#' @param k A numeric value representing the coeffecient of covariate effect on treatment probability.
#' @param marginal_hazard_ratio A numeric value which is the pre determined marginal hazard ratio
#' of the treatment.
#' @param censoring_p A numeric value which is the pre determined censoring proportion of the data.
#' @param type_censoring The distribution that should be used for the censoring. 'uniform' or 'weibull'.
#'
#' @return A list of containing: \itemize{
#'   \item data, a dataframe with n observations with survival time, censoring/event indicator,
#'   treatment indicator and covariates.
#'   \item conditional_hazard_ratio, conditional hazard ratio used to generate the data
#'   with the specified marginal hazard ratio.
#'   \item prob_T, the treatment probability of each individual.
#'   \item conditional_hazard_ratio, the pre determined marginal hazard ratio.

#' }
#' @export
#' @examples
#'   sim_res = gen_surv_data(n = 1000, lambda = 0.00002, eta = 2,s
#'   scenario = 'A', m = 1, k = 1, marginal_hazard_ratio = 0.8,
#'   censoring_p = 0.1, type_censoring = 'uniform')
#'   data = sim_res$data
gen_surv_data_ATE = function(n = 1000, counterfactual = FALSE, lambda = 0.00002, eta = 2,
                             scenario = 'A',
                             m = 1, k = 1, marginal_hazard_ratio = 0.8,
                             censoring_p = 0, type_censoring = c('uniform', 'weibull')) {
  event = rep(1,n)
  X_1 = rbinom(n = n, size = 1, prob = 0.5)
  X_3 = rbinom(n = n, size = 1, prob = 0.5)
  X_5 = rbinom(n = n, size = 1, prob = 0.5)
  X_6 = rbinom(n = n, size = 1, prob = 0.5)
  X_8 = rbinom(n = n, size = 1, prob = 0.5)
  X_9 = rbinom(n = n, size = 1, prob = 0.5)

  X_2 = rnorm(n = n)
  X_4 = rnorm(n = n)
  X_7 = rnorm(n = n)
  X_10 = rnorm(n = n)
  
  
  X_22 = X_2^2
  X_44 = X_4^2
  X_77 = X_7^2

  X_13 = X_1*X_3
  X_24 = X_2*X_4
  X_35 = X_3*X_5
  X_46 = X_4*X_6
  X_57 = X_5*X_7
  X_16 = X_1*X_6
  X_23 = X_2*X_3
  X_34 = X_3*X_4
  X_45 = X_4*X_5
  X_56 = X_5*X_6

  b_1 = 0.8 * k
  b_2 = -0.25 * k
  b_3 = 0.6 * k
  b_4 = -0.4 * k
  b_5 = -0.8 * k
  b_6 = -0.5 * k
  b_7 = 0.7 * k

  a_1 = 0.3 * m
  a_2 = -0.36 * m
  a_3 = -0.73 * m
  a_4 = -0.2 * m
  a_5 = 0.71 * m
  a_6 = -0.19 * m
  a_7 = 0.26 * m

  a_h_1 = 0.5 * m
  a_h_2 = 0.5 * m
  b_h_1 = 0.5 * k
  b_h_2 = 0.5 * k

  X_h_1 = rnorm(n = n)
  X_h_2 = rnorm(n = n)

  # Scenario A
  if(scenario == 'A'){
    logit_pr_T = b_1*X_1 + b_2*X_2 + b_3*X_3 + b_4*X_4 + b_5*X_5 + b_6*X_6 + b_7*X_7
  }
  # Scenario B
  if(scenario == 'B'){
    logit_pr_T = b_1*X_1 + b_2*X_2 + b_3*X_3 + b_4*X_4 + b_5*X_5 + b_6*X_6 + b_7*X_7 +
      b_2*X_22
  }
  # Scenario C
  if(scenario == 'C'){
    logit_pr_T = b_1*X_1 + b_2*X_2 + b_3*X_3 + b_4*X_4 + b_5*X_5 + b_6*X_6 + b_7*X_7 +
      b_2*X_22 + b_4*X_44 + b_7*X_77
  }
  # Scenario D
  if(scenario == 'D'){
    logit_pr_T = b_1*X_1 + b_2*X_2 + b_3*X_3 + b_4*X_4 + b_5*X_5 + b_6*X_6 + b_7*X_7 +
      b_1*0.5*X_13 + b_2*0.7*X_24 + b_4*0.5*X_45 + b_5*0.5*X_56
  }
  # Scenario E
  if(scenario == 'E'){
    logit_pr_T = b_1*X_1 + b_2*X_2 + b_3*X_3 + b_4*X_4 + b_5*X_5 + b_6*X_6 + b_7*X_7 +
      b_2*X_22 +
      b_1*0.5*X_13 + b_2*0.7*X_24 + b_4*0.5*X_45 + b_5*0.5*X_56
  }
  # Scenario F
  if(scenario == 'F'){
    logit_pr_T = b_1*X_1 + b_2*X_2 + b_3*X_3 + b_4*X_4 + b_5*X_5 + b_6*X_6 + b_7*X_7 +
      b_1*0.5*X_13 + b_2*0.7*X_24 + b_3*0.5*X_35 + b_4*0.7*X_46 + b_5*0.5*X_57 +
      b_1*0.5*X_16 + b_2*0.7*X_23 + b_3*0.5*X_34 + b_4*0.5*X_45 + b_5*0.5*X_56
  }
  # Scenario G
  if(scenario == 'G'){
    logit_pr_T = b_1*X_1 + b_2*X_2 + b_3*X_3 + b_4*X_4 + b_5*X_5 + b_6*X_6 + b_7*X_7 +
      b_2*X_22 + b_4*X_44 + b_7*X_77 +
      b_1*0.5*X_13 + b_2*0.7*X_24 + b_3*0.5*X_35 + b_4*0.7*X_46 + b_5*0.5*X_57 +
      b_1*0.5*X_16 + b_2*0.7*X_23 + b_3*0.5*X_34 + b_4*0.5*X_45 + b_5*0.5*X_56
  }
  prob_T <- exp(logit_pr_T) / (1 + exp(logit_pr_T))
  Treatment <- rbinom(n, 1, p = prob_T)

  U <- runif(n)
  U = log(U)
  X = as.matrix(cbind(X_1,X_2,X_3,X_4,X_5,X_6,X_7,X_8,X_9,X_10))
  
  p = ncol(X)
 
  
  LP = a_1*X_1 + a_2*X_2 + a_3*X_3 + a_4*X_4 + a_5*X_8 + a_6*X_9 + a_7*X_10

  # suppressWarnings({betaT = sbwsur:::beta_bissection(U=U,Treatment = Treatment, target = marginal_hazard_ratio,
  #                                           LP = LP, lambda = lambda, eta = eta)$beta_T})
  suppressWarnings({betaT = beta_bissection_ATE(U=U,Treatment = Treatment, target = marginal_hazard_ratio,
                                                     LP = LP, lambda = lambda, eta = eta)$beta_T})
  #latent survival time
  Y1 <- ((- U)/(lambda*exp((betaT) + LP))) ^ (1/eta)
  Y0 <- ((- U)/(lambda*exp(LP))) ^ (1/eta)
  if(counterfactual == TRUE){
    Y <- c(Y1, Y0)
    Treatment = c(rep(1, n), rep(0, n))
    LP = c(LP,LP)
    data_out = data.frame(cbind(Y, event = rep(event, 2), Treatment = Treatment, rbind(X, X)))
  } else {
    Y <- (Y1 * Treatment) + (Y0 * (1 - Treatment))
    data_out = data.frame(cbind(Y, event, Treatment, X))
    
  }
  
  data_out$Y_true <- Y
 
  if (censoring_p>0){
    # censor_values = sbwsur:::censoring_data(censoring_p = censoring_p, type = type_censoring,
    #                                betaT = betaT, Treatment = Treatment,
    #                                LP = LP, eta = eta, lambda = lambda)
    
  censoring_res <- censoring_data(censoring_p = censoring_p, type = type_censoring,
                   betaT = betaT, Treatment = Treatment,
                   LP = LP, eta = eta, lambda = lambda)
  censor_values = censoring_res$censor_values
  lambda_i = censoring_res$lambda_i
  theta = censoring_res$theta
  
  
  Ycens <- pmin(data_out$Y, censor_values)
  event <- as.numeric(data_out$Y <= censor_values)
  data_out$event <- event
  data_out$Y <- Ycens
  
  
  } 
  names(data_out) = c('Y', 'event', 'Treatment', sprintf("X%s",seq(1:p)), 'Y_true')
  
  return(list(data = data_out, conditional_hazard_ratio = exp(betaT), prob_T = prob_T, 
              marginal_hazard_ratio = marginal_hazard_ratio, LP = LP, theta = theta,
              lambda_i = lambda_i))
}
