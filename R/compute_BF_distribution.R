#' Compute the Distribution of the Bayes Factor from Carlotti and Parast (2026)
#'
#' This function calculates the probability mass function and cumulative distribution 
#' function of the Bayes factor defined in \insertCite{carlotti2026bayesian;textual}{BSET} for the
#' following hypothesis test:
#' \deqn{\begin{cases}
#' H_0: V_S = V_S^{0} \\
#' H_1: V_S > V_S^{0}
#' \end{cases}}
#' where \eqn{V_S} is the surrogate's treatment effect on \eqn{S} measured as
#' the probability
#' \deqn{V_S = P(S_{1i} > S_{0i})}
#' and \eqn{V_S^{0}} is a hypothesized value under the null hypothesis. These
#' hypotheses can be tested by fitting the following Beta-binomial model to the
#' data:
#' \deqn{\hat{V}_S \mid V_S \sim \text{Binomial} (n, V_S)}
#' \deqn{p(V_S) = \frac{\text{Beta}(V_S \mid a, b)}{\int_{1/2}^{1} \text{Beta}(v \mid a, b) \, dv}, \quad V_S \in \left(\tfrac{1}{2}, 1\right),}
#' where \eqn{\hat{V}_S} is the sample estimate of the surrogate's treatment
#' effect on \eqn{S} computed as
#' \deqn{\hat{V}_S = \frac{1}{n} \sum\limits^{n}_{i=1} I(S_{1i} > S_{0i}).}
#' The Bayes factor is then computed as the ratio of the marginal likelihoods
#' under the alternatives:
#' \deqn{BF_{n} = \frac{B(a + n \hat{V}_S,\, b + n - n \hat{V}_S) \left(1 - F_{\text{Beta}(a + n \hat{V}_S,\, b + n - n \hat{V}_S)}\!\left(\frac{1}{2}\right)\right)}{B(a, b) \left(1 - F_{\text{Beta}(a, b)}\!\left(\frac{1}{2}\right)\right)} \cdot 2^n,}
#' where \eqn{B(a, b)} is the Beta function, \eqn{F_{\text{Beta}(a,b)}} is the
#' cumulative distribution function of the \eqn{\text{Beta}(a, b)} distribution,
#' and \eqn{a} and \eqn{b} are the shape parameters of the Beta prior.
#' Given the true value of \eqn{V_S}, the distribution of the Bayes factor can
#' be computed by evaluating \eqn{BF_n} for all possible values of
#' \eqn{\hat{V}_S} and their corresponding probabilities under the Binomial
#' distribution with parameters \eqn{n} and the true value of \eqn{V_S}.
#' This function is generally not intended to be called directly by the user
#' and is instead used internally within \code{BSET_no_X} and \code{BSET_X}.
#'
#' @param n Integer. The sample size.
#' @param V_S_true Numeric. The true value of treatment effect on the surrogate.
#' @param V_S_zero Numeric. The hypothesized value of the surrogate's treatment effect under the null hypothesis (default is 0.5).
#' @param a Numeric. First shape parameter alpha for the Beta prior (default is 1).
#' @param b Numeric. Second shape parameter beta for the Beta prior (default is 1).
#' @param BF_alternative Character. The type of alternative hypothesis: either \code{"two_sided"} or \code{"greater"}.
#'
#' @return A data frame containing:
#' \itemize{
#'   \item \code{BF_values}: The possible values of the Bayes Factor.
#'   \item \code{BF_PMF}: The probability mass function for the Bayes Factor.
#'   \item \code{BF_CDF}: The cumulative distribution function for the Bayes Factor.
#' }
#' 
#' @references
#' \insertRef{carlotti2026bayesian}{BSET}
#' @importFrom stats pbeta dbinom
#' @importFrom dplyr %>% group_by summarise arrange mutate
#' @importFrom rlang .data
#' @export
compute_BF_distribution <- function(n, V_S_true, V_S_zero = 0.5, a = 1, b = 1, BF_alternative) {
  
  # Values of the Bayes Factor
  BF_values <- numeric(n + 1)
  
  for (k in 0:n) {
    # Log beta-binomial marginal likelihood under the alternative hypothesis
    if (BF_alternative == "two_sided") {
      log_marginal_H1 <- lbeta(a + k, b + n - k) - lbeta(a, b)
    } else if (BF_alternative == "greater") {
      log_marginal_H1 <- lbeta(a + k, b + n - k) - lbeta(a, b) + 
        log1p(- stats::pbeta(V_S_zero, a + k, b + n - k)) - 
        log1p(- stats::pbeta(V_S_zero, a, b))
    }
    
    # Log binomial marginal likelihood under the null hypothesis
    log_likelihood_H0 <- k * log(V_S_zero) + (n - k) * log(1 - V_S_zero)
    
    # Bayes Factor
    BF_values[k + 1] <- exp(log_marginal_H1 - log_likelihood_H0)
  }
  
  # Probabilities of observing each Bayes Factor value
  BF_probs <- stats::dbinom(
    x = 0:n,
    size = n,
    prob = V_S_true
  )
  
  # Tabulate the distribution of the Bayes Factor
  # Note: Using .data$ to avoid 'no visible binding for global variable' notes in check
  BF_distribution <- data.frame(
    BF_values = BF_values,
    BF_probs = BF_probs) %>%
    dplyr::group_by(.data$BF_values) %>% 
    dplyr::summarise(BF_PMF = sum(.data$BF_probs), .groups = "drop") %>% 
    dplyr::arrange(.data$BF_values) %>% 
    dplyr::mutate(BF_CDF = cumsum(.data$BF_PMF))
  
  return(BF_distribution)
}
