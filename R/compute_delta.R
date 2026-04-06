#' Monte Carlo Computation of the Parameter \eqn{\delta} from Parast et al. (2024)
#'
#' This function implements a Monte Carlo approach to estimate the parameter
#' \eqn{\delta} from Parast et al. (2024) . This parameter represents the
#' difference in treatment effects between the primary and surrogate outcomes,
#' both measured using the Mann-Whitney statistic.
#' 
#' @details The function processes data from a chosen data generating process,
#' computing the Mann-Whitney U statistic for both the primary outcome \eqn{Y}
#' and the surrogate \eqn{S}:
#' \deqn{\hat{U}_Y = \frac{1}{n_1 n_0} \sum\limits_{i:Z_i=1} \sum\limits_{j:Z_j=0} I(Y_i > Y_j),}
#' \deqn{\hat{U}_S = \frac{1}{n_1 n_0} \sum\limits_{i:Z_i=1} \sum\limits_{j:Z_j=0} I(S_i > S_j).}
#' Then, it calculates
#' \deqn{\hat{\delta} = \hat{U}_Y - \hat{U}_S.}
#'
#' @param MC_data A list containing:
#' \itemize{
#'   \item \code{P_observed}: A data frame or matrix with columns "Y" and "S".
#'   \item \code{Z}: Treatment assignment vector.
#'   \item \code{n1}: Number of treated units.
#'   \item \code{n0}: Number of control units.
#' }
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{U_Y}: Mann-Whitney U statistic for the primary outcome Y computed on \code{P_observed}.
#'   \item \code{U_S}: Mann-Whitney U statistic for the surrogate S computed on \code{P_observed}.
#'   \item \code{delta}: The difference \code{U_Y} - \code{U_S}.
#' }
#' 
#' @importFrom stats wilcox.test
#' @export
compute_delta <- function(MC_data) {
  
  # Standardized Wilcoxon rank-sum statistics for Y
  U_Y <- stats::wilcox.test(
    x = MC_data$P_observed[MC_data$Z == 1, "Y"],
    y = MC_data$P_observed[MC_data$Z == 0, "Y"]
  )$statistic / (MC_data$n1 * MC_data$n0)
  
  # Standardized Wilcoxon rank-sum statistics for S
  U_S <- stats::wilcox.test(
    x = MC_data$P_observed[MC_data$Z == 1, "S"],
    y = MC_data$P_observed[MC_data$Z == 0, "S"]
  )$statistic / (MC_data$n1 * MC_data$n0)
  
  # Delta calculation
  delta <- U_Y - U_S
  
  # Output as a list
  output <- list(
    U_Y = as.numeric(U_Y),
    U_S = as.numeric(U_S),
    delta = as.numeric(delta)
  )
  
  return(output)
}
