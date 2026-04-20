#' Data Generating Process without Baseline Covariates
#'
#' This function generates potential outcomes from the simulation settings
#' described in Parast et al. (2024). It creates a dataset of potential outcomes
#' \deqn{P = (Y_1, S_1, Y_0, S_0)} and observed outcomes 
#' \deqn{P_{observed} = (Y, S)} based on a random treatment assignment \eqn{Z}.
#' 
#' @details The function supports two types of data-generating processes:
#' \itemize{
#'   \item \strong{Gaussian model:} Potential outcomes are drawn from a 
#'   multivariate normal distribution:
#'   \deqn{P \sim \mathcal{N}_{4}(\mu^{*}, \Sigma^{*}).}
#'   \item \strong{Non-linear model:} Potential outcomes for the surrogate are generated from a non-Gaussian
#'   distribution, and the potential outcomes for the primary outcome are generated from a non-linear function of the surrogate plus non-Gaussian noise.
#' }
#'
#' @param n Integer. Total sample size.
#' @param p Numeric. Probability of being assigned to the treatment group (Z=1).
#' @param mu_star Numeric vector. The mean vector for \eqn{P}. Required if \code{model = "Gaussian"}.
#' @param Sigma_star Matrix. The covariance matrix for \eqn{P}. Required if \code{model = "Gaussian"}.
#' @param model Character. The type of data generation: \code{"Gaussian"} 
#'   or \code{"misspecified"}.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{Z}: Treatment assignment vector.
#'   \item \code{n1}: Number of treated units.
#'   \item \code{n0}: Number of control units.
#'   \item \code{P}: Full matrix of potential outcomes.
#'   \item \code{P_observed}: Observed outcomes \eqn{(Y, S)} corresponding to the assigned treatment \eqn{Z}.
#'   \item \code{P_unobserved}: Counterfactual outcomes under the opposite treatment.
#' }
#' This function is not a primary user-facing function of the package and
#' does not include examples.
#'
#' @importFrom Rdpack reprompt
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats rbinom rnorm
#' @export
DGP_no_X <- function(n, p, mu_star = NULL, Sigma_star = NULL, model = c("Gaussian", "misspecified")) {
  
  model <- match.arg(model)
  
  # Treatment assignment
  Z <- stats::rbinom(n = n, size = 1, prob = p)
  n1 <- sum(Z)
  n0 <- n - n1
  
  # Potential outcomes
  if (model == "Gaussian") {
    if (is.null(mu_star) || is.null(Sigma_star)) {
      stop("For Gaussian DGP, mu_star and Sigma_star must be provided.")
    }
    P <- mvtnorm::rmvnorm(n = n, mean = mu_star, sigma = Sigma_star)
    colnames(P) <- c("Y1", "S1", "Y0", "S0")
  } else if (model == "misspecified") {
    P <- matrix(data = NA, nrow = n, ncol = 4)
    colnames(P) <- c("Y1", "S1", "Y0", "S0")
    for (i in 1:n) {
      P[i, "S1"] <- exp(stats::rnorm(1, 2.5, 1.5))
      P[i, "S0"] <- exp(stats::rnorm(1, 0.5, 1.5))
      P[i, "Y1"] <- 2.0 + 1.2*sqrt(P[i, "S1"]) + 0.3*exp(P[i, "S1"]/500) + exp(stats::rnorm(1, 0, sqrt(0.3)))
      P[i, "Y0"] <- 0.8*sqrt(P[i, "S0"]) + 0.2*exp(P[i, "S0"]/50) + exp(stats::rnorm(1, 0, sqrt(0.3)))
    }
  }
  
  # Observed and unobserved outcomes (Vectorized for speed)
  P_observed <- matrix(NA, nrow = n, ncol = 2, dimnames = list(NULL, c("Y", "S")))
  P_unobserved <- matrix(NA, nrow = n, ncol = 2, dimnames = list(NULL, c("Y", "S")))
  
  # If Z=1, observed is (Y1, S1) [cols 1:2], unobserved is (Y0, S0) [cols 3:4]
  P_observed[Z == 1, ]   <- P[Z == 1, 1:2]
  P_unobserved[Z == 1, ] <- P[Z == 1, 3:4]
  
  # If Z=0, observed is (Y0, S0) [cols 3:4], unobserved is (Y1, S1) [cols 1:2]
  P_observed[Z == 0, ]   <- P[Z == 0, 3:4]
  P_unobserved[Z == 0, ] <- P[Z == 0, 1:2]
  
  return(list(
    Z = Z,
    n1 = n1,
    n0 = n0,
    P = P,
    P_observed = P_observed,
    P_unobserved = P_unobserved
  ))
}
