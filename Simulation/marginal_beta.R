#' Marginal Hazard Ratio
#'
#' This is an auxiliary function to the other functions of the package. It calculates the current
#' marginal hazard ratio using the observed value and its cunterfactural.
#'
#' @param betaT The conditional hazard ratio.
#' @param U A vector with value drawn from the uniform distribution to generate the survival time.
#' @param Treatment A vector with the treatment status of each observation.
#' @param lambda A numeric value representing the scale parameter.
#' @param eta  A numeric value representing the shape parameter.
#' @return The conditional hazard ratio for the desired marginal hazard ratio.
#' @keywords internal
#' @noRd
marginal_beta_ATE = function(betaT, U, Treatment, LP, lambda, eta){
  n = length(U)
  T1 <- rep(1,n)
  T2 <- rep(0,n)
  Y1 <- ((- U)/(lambda*exp(betaT * T1 + LP))) ^ (1/eta)
  Y2 <- ((- U)/(lambda*exp(betaT * T2 + LP))) ^ (1/eta)
  # Y1 = Y1[which(Treatment==1)]
  # Y2 = Y2[which(Treatment==1)]
  # T1 = T1[which(Treatment==1)]
  # T2 = T2[which(Treatment==1)]
  event = rep(1, 2 * length(Y1))
  data = data.frame(Y = c(Y1,Y2), T_final = c(T1,T2), event = event)
  fmla <- paste("survival::Surv(time = Y, event = event) ~ T_final ")
  cph <- survival::coxph(as.formula(fmla), data = data)
  return(exp(cph$coefficients))
}
