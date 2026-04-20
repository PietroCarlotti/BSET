#' Compute \eqn{v_S} from Carlotti and Parast (2026)
#'
#' This function determines the value \eqn{v_S} that is used to compute the
#' surrogate validation threshold \eqn{\eta} from Carlotti and Parast (2026):
#' \deqn{\eta = \max \{v_Y - v_S, 0\},}
#' where \eqn{v_Y} is the hypothesized value of the treatment effect on the
#' primary outcome (typically set equal to the estimate computed on the
#' available data) and \eqn{v_S} is the value that satisfies the following power
#' constraint:
#' \deqn{P(\text{BF}_n \geq \text{BF}_{n, \alpha} \; | \; V_S = v_S) = 1 - \beta,}
#' where \eqn{\text{BF}_{n, \alpha}} is the \eqn{(1 - \alpha)} quantile of the
#' Bayes factor distribution under the null hypothesis \eqn{V_S = V^0_{S}}, and \eqn{1 - \beta} is the desired power of the
#' test. The function computes the distribution of the Bayes factor under the null
#' hypothesis, derives the critical value \eqn{\text{BF}_{n, \alpha}}, and then
#' uses a root-finding algorithm to solve for the value of \eqn{v_S} that
#' satisfies the power constraint.
#' This function is generally not intended to be called directly by the user
#' and is instead used internally within \code{BSET_no_X} and \code{BSET_X}.
#'
#' @param n Integer. Sample size.
#' @param alpha Numeric. Type I error rate (default is 0.05).
#' @param beta Numeric. Type II error rate (default is 0.2).
#' @param V_S_zero Numeric. The hypothesized value of the surrogate's treatment effect under the null hypothesis (default is 0.5).
#' @param a Numeric. First shape parameter alpha for the Beta prior (default is 1).
#' @param b Numeric. Second shape parameter beta for the Beta prior (default is 1).
#' @param BF_alternative Character. The type of alternative hypothesis: either \code{"two_sided"} or \code{"greater"}.
#' @param root_tolerance Numeric. Tolerance level for the root-finding algorithm (default is 1e-16).
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{BF_alpha}: The critical value of the Bayes factor corresponding to the specified alpha level.
#'   \item \code{V_S_star}: The value of \eqn{v_S} that satisfies the power constraint for the surrogate validation test.
#' }
#' 
#' @importFrom stats uniroot
#' @importFrom dplyr %>% filter slice pull
#' @importFrom rlang .data
#' @export
compute_V_S_star <- function(n, alpha = 0.05, beta = 0.2, V_S_zero = 0.5, 
                             a = 1, b = 1, BF_alternative = "greater", root_tolerance = 1e-16) {
  
  # Find BF_alpha using the distribution function we formatted earlier
  BF_alpha <- compute_BF_distribution(
    n = n,
    V_S_true = V_S_zero,
    V_S_zero = V_S_zero,
    a = a,
    b = b,
    BF_alternative = BF_alternative
  ) %>%
    dplyr::filter(.data$BF_CDF >= 1 - alpha) %>%
    dplyr::slice(1) %>%
    dplyr::pull(.data$BF_values)
  
  # Objective function for uniroot
  objective_function <- function(V) {
    BF_dist <- compute_BF_distribution(
      n = n,
      V_S_true = V,
      V_S_zero = V_S_zero,
      a = a,
      b = b,
      BF_alternative = BF_alternative
    )
    
    # Output: difference between current power and target power
    output <- beta - BF_dist$BF_CDF[BF_dist$BF_values == BF_alpha]
    return(output)
  }
  
  # Define search interval based on logic
  if (BF_alternative == "two_sided") {
    interval <- c(0.5, 1) 
  } else if (BF_alternative == "greater") {
    interval <- if (V_S_zero < 0.5) c(0.5, 1) else c(V_S_zero, 1)
  }
    
  # Solve for V_S_star
  sol <- stats::uniroot(
    f = objective_function,
    interval = interval,
    tol = root_tolerance
  )
  
  # Output requested alternative
  output <- list(
    BF_alpha = BF_alpha,
    V_S_star = sol$root
  )
  
  return(output)
}
