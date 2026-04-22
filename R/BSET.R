#######################
# FINAL USER FUNCTION #
#######################

#' Bayesian Surrogate Evaluation Test from Carlotti and Parast (2026)
#'
#' This function implements the Bayesian Surrogate Evaluation Test (BSET)
#' as proposed by \insertCite{carlotti2026bayesian;textual}{BSET}. When \code{X}
#' is \code{NULL} (the default), the model is fit without covariates; when \code{X}
#' is provided, a multivariate regression model is used to adjust for baseline
#' covariates. In both cases the function fits a Bayesian model using Stan to
#' generate posterior samples for the parameters of interest, which are then
#' used to conduct a Bayesian hypothesis test for evaluating the validity of
#' the surrogate marker. A frequentist test is also performed for comparison.
#'
#' \strong{Without covariates} (\code{X = NULL}), the Bayesian model is:
#' \deqn{(Y_i, S_i) = \begin{cases}
#' P_{1i}, & \text{if } Z_i = 1 \\
#' P_{0i}, & \text{if } Z_i = 0
#' \end{cases}, \quad i = 1, \ldots, n,}
#' \deqn{P_i \overset{\text{ind.}}{\sim} \mathcal{N}_4(\mu, \Sigma), \quad i = 1, \ldots, n,}
#' \deqn{\mu \sim \mathcal{N}_4(\mu_0, \Sigma_0),}
#' \deqn{\Sigma = \text{diag}(\sigma_{1:4}) \, \Omega \, \text{diag}(\sigma_{1:4}),}
#' \deqn{\sigma_k \sim \text{Half-Normal}(0, s_k), \quad k = 1, \ldots, 4,}
#' \deqn{\Omega \sim \text{LKJ}(\tau).}
#'
#' \strong{With covariates} (\code{X} provided), the Bayesian model is:
#' \deqn{P_i \overset{\text{ind.}}{\sim} \mathcal{N}_4(B X_i, \Sigma), \quad i = 1, \ldots, n,}
#' \deqn{B = \begin{bmatrix}
#' \beta_{1} & \beta_{2} & \beta_{3} & \beta_{4} \\
#' \end{bmatrix}^\top,}
#' \deqn{\beta_{k} \sim \mathcal{N}_d (\mu_{\beta}, \Sigma_{\beta}), \quad k = 1, \ldots, 4,}
#' with \eqn{\Sigma}, \eqn{\sigma_k}, and \eqn{\Omega} as above.
#'
#' This is a primary user-facing function of the package and includes working examples below.
#'
#' @param data A data frame containing the observed data.
#' @param Y Character. Name of the outcome variable.
#' @param S Character. Name of the surrogate variable.
#' @param Z Character. Name of the treatment assignment variable.
#' @param X Character vector. Names of the covariate columns to include in the model.
#'   When \code{NULL} (the default), the model is fit without covariates and the
#'   \code{mu_0} / \code{Sigma_0} prior parameters are used. When provided, the
#'   covariate-adjusted model is fit and the \code{intercept} / \code{mu_beta} /
#'   \code{Sigma_beta} prior parameters are used instead.
#' @param delta_true The true value of delta, used to calculate frequentist coverage during simulations (optional).
#' @param theta_true The true value of theta, used to calculate Bayesian coverage during simulations (optional).
#' @param seed Random seed for reproducibility (optional. If not provided, a random seed will be generated).
#' @param n_chains Number of MCMC chains to run (default is 4).
#' @param n_iter Number of iterations per MCMC chain (default is 2000).
#' @param burn_in_ratio Proportion of iterations to discard as burn-in (default is 0.25).
#' @param a,b Prior parameters for prior Beta distribution on V_S (default is a = 1, b = 1).
#' @param alpha Type I error rate for the test (default is 0.05).
#' @param beta Type II error rate for the test (default is 0.2).
#' @param V_S_zero Value of V_S under the null hypothesis (default is 0.5).
#' @param BF_alternative Alternative hypothesis for the Bayes factor test ("greater", "less" or "two.sided").
#' @param root_tolerance Numerical tolerance for root-finding algorithms (default is 1e-16).
#' @param mu_0 Prior mean vector for the mean parameters. Used only when \code{X = NULL}
#'   (default is a vector of zeros of length 4).
#' @param Sigma_0 Prior covariance matrix for the mean vector. Used only when \code{X = NULL}
#'   (default is a 4x4 identity matrix).
#' @param intercept Logical. Whether to include an intercept in the regression model.
#'   Used only when \code{X} is provided (default is \code{TRUE}).
#' @param mu_beta Prior mean vector for the regression coefficients. Used only when \code{X}
#'   is provided (default is a vector of zeros with length equal to the number of covariates,
#'   including the intercept if \code{intercept = TRUE}).
#' @param Sigma_beta Prior covariance matrix for the regression coefficients. Used only when
#'   \code{X} is provided (default is a diagonal matrix with variance 10, with dimensions equal
#'   to the number of covariates, including the intercept if \code{intercept = TRUE}).
#' @param s Prior scale parameters for the error variance (default is a vector of ones of length 4).
#' @param tau Prior parameter for the LKJ correlation distribution (default is 1).
#' @param plot Logical. Whether to plot the posterior distribution of theta (default is FALSE).
#' @param mute Logical. Whether to suppress Stan output during model fitting (default is TRUE).
#' @param parallel Logical. Whether to use parallel processing for MCMC sampling (default is TRUE).
#'
#' @return A list containing:
#' \itemize{
#'  \item \code{stan_cores}: The number of CPU cores used for MCMC sampling.
#'  \item \code{Bayesian_test}: A list containing the results of the Bayesian test, including posterior samples and credible intervals.
#'  \item \code{frequentist_test}: A list containing the results of the frequentist test, including point estimates and confidence intervals.
#'  \item \code{theta_posterior_plot}: A \code{ggplot} object showing the posterior distribution of \eqn{\theta}, with vertical lines indicating the credible interval, the \eqn{\eta} threshold, and the true values of \eqn{\delta} and \eqn{\theta} (if provided).
#'  }
#' @examples
#' # Generate data from the perfect surrogate setting of Parast et al. (2024)
#' set.seed(123)
#' data_no_X <- DGP_no_X(
#'   n = 100,
#'   p = 0.5,
#'   mu_star = c(6, 6, 2.5, 2.5),
#'   Sigma_star = kronecker(diag(2), matrix(c(3, 3, 3, 3.1), 2, 2)),
#'   model = "Gaussian"
#' )
#'
#' # Prepare the data frame
#' df_no_X <- data.frame(
#'   Y = data_no_X$P_observed[, "Y"],
#'   S = data_no_X$P_observed[, "S"],
#'   Z = data_no_X$Z
#' )
#'
#' # Run BSET without covariates (computationally intensive)
#' \donttest{
#' result_no_X <- BSET(
#'   data = df_no_X,
#'   Y = "Y",
#'   S = "S",
#'   Z = "Z",
#'   seed = 123,
#'   n_chains = 2,
#'   n_iter = 500
#' )
#' }
#'
#' # Generate data from the setting of Carlotti and Parast (2026) with a binary covariate
#' set.seed(123)
#' data_X <- DGP_X_binary(
#'   n = 100,
#'   p = 0.5,
#'   q = 0.5,
#'   mu_0 = c(5, 5, 0, 0),
#'   mu_1 = c(5, -5, 0, -10),
#'   Sigma_0 = kronecker(diag(2), matrix(c(1, 1, 1, 2), 2, 2)),
#'   Sigma_1 = kronecker(diag(2), matrix(c(1, 1, 1, 2), 2, 2))
#' )
#'
#' # Prepare the data frame
#' df_X <- data.frame(
#'   Y = data_X$P_observed[, "Y"],
#'   S = data_X$P_observed[, "S"],
#'   Z = data_X$Z,
#'   X = data_X$X
#' )
#'
#' # Run BSET with covariates (computationally intensive)
#' \donttest{
#' result_X <- BSET(
#'   data = df_X,
#'   Y = "Y",
#'   S = "S",
#'   Z = "Z",
#'   X = "X",
#'   seed = 123,
#'   n_chains = 2,
#'   n_iter = 500
#' )
#' }
#' @references
#' \insertRef{carlotti2026bayesian}{BSET}
#'
#' \insertRef{parast2024rank}{BSET}
#' @importFrom rlang .data
#' @export
BSET <- function(
    data, Y, S, Z, X = NULL, delta_true = NULL, theta_true = NULL,
    seed = NULL, n_chains = 4, n_iter = 2000, burn_in_ratio = 0.25,
    a = 1, b = 1, alpha = 0.05, beta = 0.2, V_S_zero = 0.5, BF_alternative = "greater", root_tolerance = 1e-16,
    mu_0 = rep(0, 4), Sigma_0 = diag(4),
    intercept = TRUE, mu_beta = NULL, Sigma_beta = NULL,
    s = rep(1, 4), tau = 1,
    plot = FALSE, mute = TRUE, parallel = TRUE
) {

  use_X <- !is.null(X)

  # Ensure all requested columns actually exist in the data
  required_cols <- c(Y, S, Z)
  if (use_X) required_cols <- c(required_cols, X)
  missing_cols <- setdiff(required_cols, names(data))

  if (length(missing_cols) > 0) {
    stop(paste("The following columns are missing from the data frame:", paste(missing_cols, collapse = ", ")))
  }

  # Extract variables from the data frame
  Y <- data[[Y]]
  S <- data[[S]]
  Z <- data[[Z]]

  # Check for NAs
  if (anyNA(data[, required_cols])) {
    warning("Data contains missing values. Filter them first.")
  }

  # Sample size
  n <- nrow(data)

  # Covariate-specific processing
  if (use_X) {
    X <- as.matrix(data[X])
    d <- ncol(X)
    if (intercept) {
      X <- cbind(1, X)
      d <- d + 1
    }
    if (is.null(mu_beta))    mu_beta    <- rep(0, d)
    if (is.null(Sigma_beta)) Sigma_beta <- 10 * diag(d)
  }

  # Burn-in iterations
  burn_in <- floor(n_iter * burn_in_ratio)

  # Total number of posterior samples
  n_samples <- n_chains * (n_iter - burn_in)

  # Compute V_S_star
  V_S_star <- compute_V_S_star(
    n = n,
    alpha = alpha,
    beta = beta,
    V_S_zero = V_S_zero,
    a = a,
    b = b,
    BF_alternative = BF_alternative,
    root_tolerance = root_tolerance
  )$V_S_star

  # Parallel vs sequential execution
  if (parallel) {
    # Detect available cores (cap at 2 during CRAN checks)
    available_cores <- if (isTRUE(as.logical(Sys.getenv("_R_CHECK_LIMIT_CORES_", "FALSE")))) {
      2L
    } else {
      parallel::detectCores() - 1
    }

    # Maximum allowed PSOCK connections. R default is 128
    max_connections <- getOption("future.connections.max", 128) - 4

    # Number of cores to use
    stan_cores <- min(available_cores, max_connections, n_chains)

    # Plan parallel execution
    future::plan(future::multisession, workers = stan_cores)

    # Shut down workers on exit
    on.exit(future::plan(future::sequential), add = TRUE)
  } else {
    stan_cores <- 1
  }

  # Ensure compiled models are saved and reused
  rstan::rstan_options(auto_write = TRUE)

  # Set the seed for reproducibility
  stan_seed <- if (!is.null(seed)) seed else sample.int(.Machine$integer.max, 1)

  # Select Stan model and data based on whether covariates are provided
  if (use_X) {
    stan_file <- system.file("stan", "BSET_X.stan", package = "BSET")
    stan_data <- list(
      n = n,
      d = d,
      P_observed = cbind(Y, S),
      Z = Z,
      X = X,
      mu_beta = mu_beta,
      Sigma_beta = Sigma_beta,
      s = s,
      tau = tau
    )
  } else {
    stan_file <- system.file("stan", "BSET_no_X.stan", package = "BSET")
    stan_data <- list(
      n = n,
      P_observed = cbind(Y, S),
      Z = Z,
      mu_0 = mu_0,
      Sigma_0 = Sigma_0,
      s = s,
      tau = tau
    )
  }

  # Compile the model if it doesn't exist, otherwise load the pre-compiled model
  stan_model <- rstan::stan_model(file = stan_file)

  # Fit Stan model
  if (mute) {
    stan_fit <- suppressMessages(
      suppressWarnings(
        rstan::sampling(
          object = stan_model,
          data = stan_data,
          chains = n_chains,
          cores = stan_cores,
          iter = n_iter,
          warmup = burn_in,
          seed = stan_seed,
          control = list(adapt_delta = 0.999, max_treedepth = 15),
          refresh = 0
        )
      )
    )
  } else {
    stan_fit <- rstan::sampling(
      object = stan_model,
      data = stan_data,
      chains = n_chains,
      cores = stan_cores,
      iter = n_iter,
      warmup = burn_in,
      seed = stan_seed,
      control = list(adapt_delta = 0.999, max_treedepth = 15),
      refresh = 0
    )
  }

  # Extract posterior samples
  HMC_samples <- rstan::extract(stan_fit)

  # Reorder posterior draws: [n_subjects, variables, n_samples]
  P_HMC <- aperm(HMC_samples$P, c(2, 3, 1))

  # Bayesian test
  Bayesian_test_results <- Bayesian_test(
    P_MCMC = P_HMC,
    alpha = alpha,
    V_S_star = V_S_star,
    theta_true = if (!is.null(theta_true)) theta_true else NULL
  )

  # Frequentist test
  frequentist_test_results <- frequentist_test(
    P_observed = cbind(Y, S),
    Z = Z,
    alpha = alpha,
    beta = beta,
    delta_true = if (!is.null(delta_true)) delta_true else NULL
  )

  # Posterior plot for theta
  if (plot) {
    theta_posterior_plot <- ggplot2::ggplot(data.frame(theta = Bayesian_test_results$theta_MCMC), mapping = ggplot2::aes(x = .data$theta)) +
      ggplot2::geom_histogram(
        binwidth = 0.02,
        fill = "lightblue",
        color = "black",
        alpha = 0.8
      ) +
      ggplot2::geom_vline(
        xintercept = Bayesian_test_results$CI[2],
        color = "blue",
        linetype = "solid",
        linewidth = 1
      ) +
      ggplot2::geom_vline(
        xintercept = Bayesian_test_results$epsilon,
        color = "darkgreen",
        linetype = "solid",
        linewidth = 1
      ) +
      ggplot2::coord_cartesian(xlim = c(-1, 1)) +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        x = expression(theta),
        y = "Frequency"
      )

    if (!is.null(delta_true)) {
      theta_posterior_plot <- theta_posterior_plot +
        ggplot2::geom_vline(
          xintercept = delta_true,
          color = "orange",
          linetype = "solid",
          linewidth = 1
        )
    }

    if (!is.null(theta_true)) {
      theta_posterior_plot <- theta_posterior_plot +
        ggplot2::geom_vline(
          xintercept = theta_true,
          color = "red",
          linetype = "solid",
          linewidth = 1
        )
    }
  }

  # Output
  output <- list(
    stan_cores = stan_cores,
    Bayesian_test = Bayesian_test_results,
    frequentist_test = frequentist_test_results,
    theta_posterior_plot = if (plot) theta_posterior_plot else NULL
  )

  return(output)
}
