#' Data Generating Process with a Binary Covariate
#'
#' This function generates potential outcomes from a data generating process
#' similar to the one described in Parast et al. (2024), but with the addition
#' of a binary covariate X. It creates a dataset of potential outcomes
#' \deqn{P = (Y_1, S_1, Y_0, S_0)} and observed outcomes
#' \deqn{P_{observed} = (Y, S)} based on a random treatment assignment \eqn{Z}.
#' 
#' @details The potential outcomes are generated from multivariate normal
#' distributions with different mean vectors and covariance matrices depending
#' on the value of \eqn{X}. Specifically:
#' \deqn{P \mid X = 0 \sim \mathcal{N}_{4}(\mu^{0}, \Sigma^{0}), \\
#' P \mid X = 1 \sim \mathcal{N}_{4}(\mu^{1}, \Sigma^{1}).}
#'
#' @param n Integer. Total sample size.
#' @param p Numeric. Probability of being assigned to the treatment group \eqn{(Z=1)}.
#' @param q Numeric. Probability of the binary covariate X being 1.
#' @param mu_0 Numeric vector. Mean vector for \eqn{P} when \eqn{X=0}.
#' @param mu_1 Numeric vector. Mean vector for \eqn{P} when \eqn{X=1}.
#' @param Sigma_0 Matrix. Covariance matrix for \eqn{P} when \eqn{X=0}.
#' @param Sigma_1 Matrix. Covariance matrix for \eqn{P} when \eqn{X=1}.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{X}: The binary covariate vector.
#'   \item \code{n_X1}: Number of units with \eqn{X=1}.
#'   \item \code{n_X0}: Number of units with \eqn{X=0}.
#'   \item \code{Z}: Treatment assignment vector.
#'   \item \code{n1}: Number of treated units.
#'   \item \code{n0}: Number of control units.
#'   \item \code{P}: Full matrix of potential outcomes.
#'   \item \code{P_observed}: Observed outcomes \eqn{(Y, S)} corresponding to the assigned treatment \eqn{Z}.
#'   \item \code{P_unobserved}: Counterfactual outcomes under the opposite treatment.
#' }
#' 
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats rbinom
#' @export
DGP_X_binary <- function(n, p, q, mu_0, mu_1, Sigma_0, Sigma_1) {
  
  # Binary covariate
  X <- stats::rbinom(n = n, size = 1, prob = q)
  n_X1 <- sum(X)
  n_X0 <- n - n_X1
  
  # Potential outcomes generation
  P_X0 <- mvtnorm::rmvnorm(n = n, mean = mu_0, sigma = Sigma_0)
  P_X1 <- mvtnorm::rmvnorm(n = n, mean = mu_1, sigma = Sigma_1)
  
  colnames(P_X0) <- colnames(P_X1) <- c("Y1", "S1", "Y0", "S0")
  
  # Composite potential outcomes matrix
  P <- X * P_X1 + (1 - X) * P_X0
  
  # Treatment assignment
  Z <- stats::rbinom(n = n, size = 1, prob = p)
  n1 <- sum(Z)
  n0 <- n - n1
  
  # Observed and Unobserved (Vectorized logic)
  # When Z=1, we observe cols 1-2 (Y1, S1). When Z=0, we observe cols 3-4 (Y0, S0).
  P_observed <- Z * P[, c("Y1", "S1")] + (1 - Z) * P[, c("Y0", "S0")]
  P_unobserved <- Z * P[, c("Y0", "S0")] + (1 - Z) * P[, c("Y1", "S1")]
  
  colnames(P_observed) <- colnames(P_unobserved) <- c("Y", "S")
  
  return(list(
    X = X,
    n_X1 = n_X1,
    n_X0 = n_X0,
    Z = Z,
    n1 = n1,
    n0 = n0,
    P = P,
    P_observed = P_observed,
    P_unobserved = P_unobserved
  ))
}

#' Data Generating Process with a Continuous Covariate
#' 
#' This function generates potential outcomes from a data generating process
#' similar to the one described in Parast et al. (2024), but with the addition
#' of a continuous covariate X. It creates a dataset of potential outcomes
#' \deqn{P = (Y_1, S_1, Y_0, S_0)} and observed outcomes
#' \deqn{P_{observed} = (Y, S)} based on a random treatment assignment \eqn{Z}.
#' 
#' @details The potential outcomes are generated from multivariate normal
#' distributions with mean vector and covariance matrix that depend on the
#' value of \eqn{X}. Specifically, the mean vector is a linear function of \eqn{X}:
#' \deqn{\mu(X) = x \cdot (\beta_{Y1}, \beta_{S1}, \beta_{Y0}, \beta_{S0})^T,}
#' and the covariance matrix is constant across values of \eqn{X}:
#' \deqn{\Sigma(X) = \Sigma.}
#' 
#' @param n Integer. Total sample size.
#' @param p Numeric. Probability of being assigned to the treatment group \eqn{(Z=1)}.
#' @param beta Numeric vector. Coefficients for the linear function of the mean vector.
#' @param Sigma Matrix. Covariance matrix for the potential outcomes.
#' @param m Numeric. Mean of the continuous covariate X.
#' @param s Numeric. Standard deviation of the continuous covariate X.
#' 
#' @return A list containing:
#' \itemize{
#'  \item \code{X}: The continuous covariate vector.
#'  \item \code{Z}: Treatment assignment vector.
#'  \item \code{n1}: Number of treated units.
#'  \item \code{n0}: Number of control units.
#'  \item \code{P}: Full matrix of potential outcomes.
#'  \item \code{P_observed}: Observed outcomes \eqn{(Y, S)} corresponding to the assigned treatment \eqn{Z}.
#'  \item \code{P_unobserved}: Counterfactual outcomes under the opposite treatment.
#' }
#' 
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats rbinom rnorm
#' @export
DGP_X_Gaussian <- function(n, p, beta, Sigma, m, s) {
  
  # Continuous covariate
  X <- stats::rnorm(n = n, mean = m, sd = s)
  
  # Potential outcomes generation
  mu_X <- X %*% t(beta)
  P <- matrix(data = NA, nrow = n, ncol = 4)
  colnames(P) <- c("Y1", "S1", "Y0", "S0")
  
  for (i in 1:n) {
    P[i, ] <- mvtnorm::rmvnorm(n = 1, mean = mu_X[i, ], sigma = Sigma)
  }
  
  # Treatment assignment
  Z <- stats::rbinom(n = n, size = 1, prob = p)
  n1 <- sum(Z)
  n0 <- n - n1
  
  # Observed and Unobserved (Vectorized logic)
  P_observed <- Z * P[, c("Y1", "S1")] + (1 - Z) * P[, c("Y0", "S0")]
  P_unobserved <- Z * P[, c("Y0", "S0")] + (1 - Z) * P[, c("Y1", "S1")]
  
  colnames(P_observed) <- colnames(P_unobserved) <- c("Y", "S")
  
  return(list(
    X = X,
    Z = Z,
    n1 = n1,
    n0 = n0,
    P = P,
    P_observed = P_observed,
    P_unobserved = P_unobserved
  ))
}
