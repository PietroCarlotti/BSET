#' Monte Carlo Computation of the Parameter \eqn{\theta} from Carlotti and Parast (2026)
#'
#' This function implements a Monte Carlo approach to estimate the parameter
#' \eqn{\theta} from Carlotti and Parast (2026). This parameter represents the
#' difference in treatment effects between the primary and surrogate outcomes,
#' both measured using the probability that the treated outcome is larger
#' than the control outcome.
#' 
#' @details The function processes data from a chosen data generating process,
#' computing the sample probabilities for both the primary outcome \eqn{Y} and
#' the surrogate \eqn{S}:
#' \deqn{\hat{V}_Y = \frac{1}{n} \sum\limits^{n}_{i=1} I(Y_{1i} > Y_{0i}),}
#' \deqn{\hat{V}_S = \frac{1}{n} \sum\limits^{n}_{i=1} I(S_{1i} > S_{0i}).}
#' Then, it calculates
#' \deqn{\hat{\theta} = \hat{V}_Y - \hat{V}_S.}
#' This function is generally not intended to be called directly by the user
#' and is instead used internally within \code{BSET_no_X} and \code{BSET_X}.
#'
#' @param MC_data A list containing:
#' \itemize{
#'   \item \code{P}: A matrix or data frame of potential outcomes with columns 
#'     "Y1", "Y0", "S1", and "S0".
#' }
#'
#' @return A list containing the true values:
#' \itemize{
#'   \item \code{V_Y}: The Monte Carlo estimate of \eqn{P(Y_{1i} > Y_{0i})} computed on \code{P}.
#'   \item \code{V_S}: The Monte Carlo estimate of \eqn{P(S_{1i} > S_{0i})} computed on \code{P}.
#'   \item \code{theta}: The difference \code{V_Y} - \code{V_S}.
#' }
#' 
#' @export
compute_theta <- function(MC_data) {
  
  # Calculate true probabilities from potential outcomes
  # V_Y = P(Y1 > Y0)
  V_Y <- mean(MC_data$P[, "Y1"] > MC_data$P[, "Y0"])
  
  # V_S = P(S1 > S0)
  V_S <- mean(MC_data$P[, "S1"] > MC_data$P[, "S0"])
  
  # Theta (True discrepancy)
  theta <- V_Y - V_S
  
  # Output
  output <- list(
    V_Y = V_Y,
    V_S = V_S,
    theta = theta
  )
  
  return(output)
}
