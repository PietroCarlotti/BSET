#' True estimands for Parast et al. (2024)
#'
#' A dataset containing the Monte Carlo estimates of the parameters of interest
#' for the simulation settings considered in Parast et al. (2024).
#'
#' @format A data frame with 4 rows and 7 columns:
#' \describe{
#'   \item{setting}{The index of the simulation setting.}
#'   \item{U_Y_MC}{Numeric vector of Monte Carlo estimates.}
#'   \item{U_S_MC}{Numeric vector of Monte Carlo estimates.}
#'   \item{delta_MC}{Numeric vector of Monte Carlo estimates.}
#'   \item{V_Y_MC}{Numeric vector of Monte Carlo estimates.}
#'   \item{V_S_MC}{Numeric vector of Monte Carlo estimates.}
#'   \item{theta_MC}{Numeric vector of Monte Carlo estimates.}
#' }
#' @source Computed using the function \code{compute_estimands_Parast_et_al_2024(1000000)}.
"estimands_Parast_et_al_2024"

#' True estimands for Carlotti and Parast (2026)
#'
#' A dataset containing the Monte Carlo estimates of the parameters of interest
#' for the simulation settings considered in Carlotti and Parast (2026).
#'
#' @format A data frame with 1 row and 6 columns:
#' \describe{
#'   \item{U_Y_MC}{Numeric vector of Monte Carlo estimates.}
#'   \item{U_S_MC}{Numeric vector of Monte Carlo estimates.}
#'   \item{delta_MC}{Numeric vector of Monte Carlo estimates.}
#'   \item{V_Y_MC}{Numeric vector of Monte Carlo estimates.}
#'   \item{V_S_MC}{Numeric vector of Monte Carlo estimates.}
#'   \item{theta_MC}{Numeric vector of Monte Carlo estimates.}
#' }
#' @source Computed using the function \code{compute_estimands_Carlotti_and_Parast_2026(1000000)}.
"estimands_Carlotti_and_Parast_2026"

#' Simulation grid for Parast et al. (2024)
#'
#' A dataset containing the grid of simulation results for Parast et al. (2024).
#'
#' @format A data frame with multiple columns:
#' \describe{
#'   \item{setting}{The index of the simulation setting.}
#'   \item{simulation}{The index of the simulation run.}
#'   \item{n_sample}{The sample size used.}
#'   \item{n_chains}{The number of MCMC chains used.}
#'   \item{n_iterations}{The number of MCMC iterations.}
#'   \item{n_simulations}{The total number of simulations.}
#'   \item{timestamp}{The timestamp of execution.}
#'   \item{seed}{The random seed used.}
#'   \item{BF_alternative}{The alternative hypothesis used for the Bayes factor.}
#'   \item{Bayesian_epsilon}{The threshold for the Bayesian test.}
#'   \item{frequentist_epsilon}{The threshold for the frequentist test.}
#'   \item{Bayesian_CI_upper}{The upper bound of the Bayesian credible interval.}
#'   \item{frequentist_CI_upper}{The upper bound of the frequentist confidence interval.}
#'   \item{Bayesian_coverage}{Logical. Indicates if the true value is in the credible interval.}
#'   \item{frequentist_coverage}{Logical. Indicates if the true value is in the confidence interval.}
#'   \item{Bayesian_power}{Logical. Indicates if the Bayesian test rejected the null.}
#'   \item{frequentist_power}{Logical. Indicates if the frequentist test rejected the null.}
#' }
#' @source Computed using
#'   \code{Parast_et_al_2024_simulations(seed = 124, n_simulations = 50)}.
"Parast_et_al_2024_simulations_grid"

#' Simulation grid for Carlotti and Parast (2026)
#'
#' A dataset containing the grid of simulation results for Carlotti and
#' Parast (2026).
#'
#' @format A data frame with multiple columns:
#' \describe{
#'   \item{setting}{The index of the simulation setting.}
#'   \item{simulation}{The index of the simulation run.}
#'   \item{n_sample}{The sample size used.}
#'   \item{n_chains}{The number of MCMC chains used.}
#'   \item{n_iterations}{The number of MCMC iterations.}
#'   \item{n_simulations}{The total number of simulations.}
#'   \item{timestamp}{The timestamp of execution.}
#'   \item{seed}{The random seed used.}
#'   \item{BF_alternative}{The alternative hypothesis used for the Bayes factor.}
#'   \item{Bayesian_epsilon}{The threshold for the Bayesian test.}
#'   \item{frequentist_epsilon}{The threshold for the frequentist test.}
#'   \item{Bayesian_CI_upper}{The upper bound of the Bayesian credible interval.}
#'   \item{frequentist_CI_upper}{The upper bound of the frequentist confidence interval.}
#'   \item{Bayesian_coverage}{Logical. Indicates if the true value is in the credible interval.}
#'   \item{frequentist_coverage}{Logical. Indicates if the true value is in the confidence interval.}
#'   \item{Bayesian_power}{Logical. Indicates if the Bayesian test rejected the null.}
#'   \item{frequentist_power}{Logical. Indicates if the frequentist test rejected the null.}
#' }
#' @source Computed using
#'   \code{Carlotti_and_Parast_2026_simulations(seed = 123,
#'   n_simulations = 500)}.
"Carlotti_and_Parast_2026_simulations_grid"
