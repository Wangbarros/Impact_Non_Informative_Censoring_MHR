#' Marginal hazard ratio bissection
#'
#' This is an auxiliary function to the other functions of the package. It applies the bissection method
#' to find the conditional hazard ratio that corresponds to the pre determined marginal hazard ratio.
#'
#' @param target The desired marginal hazard ratio.
#' @param U A vector with value drawn from the uniform distribution to generate the survival time.
#' @param Treatment A vector with the treatment status of each observation.
#' @param lambda A numeric value representing the scale parameter.
#' @param eta  A numeric value representing the shape parameter.
#' @return The marginal hazard ratio.
#' @keywords internal
#' @noRd
beta_bissection_ATE = function(U, Treatment, target, LP, lambda, eta){
  a_b = c(-50, 40)
  f_a_b = c(Inf,Inf)
  i=0
  convergence = FALSE
  max_iter = 200
  while (i<max_iter && !convergence){
    #print(a_b)
    c = (a_b[1]+a_b[2])/2
    f_a = marginal_beta_ATE(betaT=a_b[1], U=U, Treatment=Treatment,
                        LP=LP, lambda=lambda, eta=eta)
    f_b = marginal_beta_ATE(betaT=a_b[2], U=U, Treatment=Treatment,
                        LP=LP, lambda=lambda, eta=eta)
    f_c = marginal_beta_ATE(betaT=c, U=U, Treatment=Treatment,
                        LP=LP, lambda=lambda, eta=eta)
    f_a = f_a - target
    f_b = f_b - target
    f_c = f_c - target

    if (sign(f_c) == sign(f_a)){
      a_b[1] = c
    }
    if (sign(f_c) == sign(f_b)){
      a_b[2] = c
    }
    if(abs(a_b[1]-a_b[2])< 10^-6){
      convergence = TRUE
    }
    i = i+1
    if(i==max_iter-1){
      stop('Warning: beta did not converge')
    }
  }

  beta_T = c
  return (list(beta_T = beta_T, a_b = a_b))
}
