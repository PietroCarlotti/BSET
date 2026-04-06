#' Frequentist Test from Parast et al. (2024)
#'
#' This function performs a frequentist test for surrogate evaluation based on 
#' the rank-based methodology of Parast et al. (2024). It calculates confidence
#' intervals, the \eqn{\varepsilon} threshold used, and determines if the surrogate is
#' valid.
#'
#' @param P_observed Observed outcomes \eqn{(Y, S)} corresponding to the assigned treatment \eqn{Z}.
#' @param Z Treatment assignment vector.
#' @param alpha Numeric. Significance level for the confidence interval (default is 0.05).
#' @param beta Numeric. Type II error rate (default is 0.2).
#' @param delta_true Numeric (optional). The true value of \eqn{\delta}, used to calculate 
#'   frequentist coverage during simulations.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{CI}: The calculated confidence interval for \eqn{\delta}.
#'   \item \code{epsilon}: The \eqn{\varepsilon} threshold value used in the test.
#'   \item \code{coverage}: Logical indicating if \code{delta_true} falls within the \code{CI} 
#'     (if \code{delta_true} is provided).
#'   \item \code{power}: Logical indicating if the upper bound of \code{CI} is below \code{epsilon}, which indicates that the test identifies the surrogate as valid.
#' }
#' 
#' @importFrom SurrogateRank test.surrogate
#' @export
frequentist_test <- function(P_observed, Z, alpha = 0.05, beta = 0.2, delta_true = NULL) {
  
  # Surrogacy test from Parast et al. (2024)
  Parast_test <- SurrogateRank::test.surrogate(
    yone = P_observed[Z == 1, "Y"],
    yzero = P_observed[Z == 0, "Y"],
    sone = P_observed[Z == 1, "S"],
    szero = P_observed[Z == 0, "S"],
    power.want.s = 1 - beta
  )
  
  # Calculations for coverage
  if (!is.null(delta_true)) {
    # Check if true delta is below the upper bound of the CI
    frequentist_coverage <- delta_true < Parast_test$ci.delta[2]
  } else {
    frequentist_coverage <- NA
  }
  
  # Return as a list
  list(
    delta = Parast_test$delta.estimate, 
    CI = Parast_test$ci.delta,
    epsilon = Parast_test$epsilon.used,
    coverage = frequentist_coverage,
    power = Parast_test$is.surrogate
  )
}
