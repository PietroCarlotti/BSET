#' Bayesian Test from Carlotti and Parast (2026)
#'
#' This function performs a Bayesian test for surrogate evaluation based on the
#' imputation-based methodology of \insertCite{carlotti2026bayesian;textual}{BSET}. It calculates
#' credible intervals, the \eqn{\eta} threshold used, and determines if the
#' surrogate is valid.
#' This function is generally not intended to be called directly by the user
#' and is instead used internally within \code{BSET_no_X} and \code{BSET_X}.
#'
#' @param P_MCMC A three-dimensional array of potential outcomes obtained via MCMC sampling with dimensions \code{[n_subjects, variables, n_samples]}. The variables should correspond to \eqn{(Y_1, S_1, Y_0, S_0)}.
#' @param alpha Numeric. Significance level for the credible interval (default is 0.05).
#' @param V_S_star Numeric. The value of \eqn{v_S} that satisfies the power constraint for the surrogate validation test.
#' @param theta_true Numeric (optional). The true value of \eqn{\eta}, used to calculate frequentist coverage during simulations.
#' @param V_Y_true Numeric (optional). The true value of \eqn{V_Y}, used to calculate the surrogate validation threshold \eqn{\eta} during simulations. If not provided, it will be estimated from the MCMC samples.
#' @return A list containing:
#' \itemize{
#'   \item \code{V_Y_MCMC}: Posterior draws for \eqn{V_Y}.
#'   \item \code{V_S_MCMC}: Posterior draws for \eqn{V_S}.
#'   \item \code{theta_MCMC}: Posterior draws for \eqn{\theta = V_Y - V_S}.
#'   \item \code{CI}: The calculated credible interval for \eqn{\theta}.
#'   \item \code{threshold}: The \eqn{\eta} threshold value used in the test.
#'   \item \code{coverage}: Logical indicating if \code{theta_true} falls within the \code{CI} (if \code{theta_true} is provided).
#'   \item \code{power}: Logical indicating if the upper bound of \code{CI} is below \code{threshold}, which indicates that the test identifies the surrogate as valid.
#' }
#' @references
#' \insertRef{carlotti2026bayesian}{BSET}
#' @export
Bayesian_test <- function(P_MCMC, alpha = 0.05, V_S_star, theta_true = NULL, V_Y_true = NULL) {
  
  # Determine number of samples iterations
  n_samples <- dim(P_MCMC)[3]
  
  # Set column names for clarity
  colnames(P_MCMC) <- c("Y1", "S1", "Y0", "S0")
  
  # Initialize vectors to store posterior draws
  V_Y_MCMC <- numeric(n_samples)
  V_S_MCMC <- numeric(n_samples)
  theta_MCMC <- numeric(n_samples)
  
  # Loop over samples
  for (k in 1:n_samples) {
    V_Y_MCMC[k] <- mean(P_MCMC[, "Y1", k] > P_MCMC[, "Y0", k])
    V_S_MCMC[k] <- mean(P_MCMC[, "S1", k] > P_MCMC[, "S0", k])
    theta_MCMC[k] <- V_Y_MCMC[k] - V_S_MCMC[k]
  }
  
  # Bayesian credible interval
  Bayes_CI <- c(-1, stats::quantile(theta_MCMC, probs = 1 - alpha))
  
  # Calculations
  Bayes_threshold <- max(mean(V_Y_MCMC) - V_S_star, 0)
  Bayes_power <- Bayes_CI[2] < Bayes_threshold
  
  if (!is.null(theta_true)) {
    Bayes_coverage <- theta_true < Bayes_CI[2]
  } else {
    Bayes_coverage <- NA
  }
  
  # Return as a list
  list(
    V_Y_MCMC = V_Y_MCMC,
    V_S_MCMC = V_S_MCMC,
    theta_MCMC = theta_MCMC,
    CI = Bayes_CI,
    threshold = Bayes_threshold,
    coverage = Bayes_coverage,
    power = Bayes_power
  )
}
