#' Run Simulations for the Settings in Carlotti and Parast (2026)
#'
#' This function executes the simulation study described in Carlotti and Parast
#' (2026) across two covariate settings: a binary covariate setting (setting 1,
#' \code{X_binary}) and a Gaussian covariate setting (setting 2,
#' \code{X_Gaussian}). It performs the frequentist and Bayesian surrogate
#' evaluation tests from Parast et al. (2024) and Carlotti and Parast (2026),
#' respectively, and returns a structured data frame of the results. The
#' execution is parallelized using the \code{future} framework.
#'
#' @details
#' The function runs \code{n_simulations} simulations for each of the two
#' settings, generating data using \code{DGP_X_binary} (setting 1) or
#' \code{DGP_X_Gaussian} (setting 2) and processes the results using
#' \code{BSET}. The simulation utilizes MCMC sampling via \code{rstan} for the
#' Bayesian estimation components. Note that it relies on an external object
#' \code{estimands_Carlotti_and_Parast_2026} that contains the true values of
#' the parameters computed using the function
#' \code{compute_estimands_Carlotti_and_Parast_2026}.
#'
#' @param seed Numeric. A random seed for reproducibility of the simulations.
#' @param n_simulations Numeric. The number of simulation runs to execute per setting.
#' @param parallel Logical. Whether to run the simulations in parallel or sequentially (default is TRUE).
#' @return A data frame containing:
#' \itemize{
#'   \item \code{setting}: The index of the simulation setting (1 for binary covariate, 2 for Gaussian covariate).
#'   \item \code{simulation}: The individual simulation run ID.
#'   \item \code{Bayesian_epsilon}, \code{frequentist_epsilon}: Numeric. Surrogate validation thresholds for the Bayesian and frequentist tests, respectively.
#'   \item \code{Bayesian_coverage}, \code{frequentist_coverage}: Boolean. Indicates whether the true surrogate effect is within the credible or confidence interval for the Bayesian and frequentist tests, respectively.
#'   \item \code{Bayesian_power}, \code{frequentist_power}: Boolean. Indicates whether the upper bound of the credible or confidence interval is below the threshold for the Bayesian and frequentist tests, respectively.
#'   \item \code{n_sample}, \code{n_simulations}, \code{n_chains}, \code{n_iterations}: Numeric. Metadata about the simulation parameters.
#'   \item \code{timestamp}: Character. Date and time when the simulation was completed.
#'   \item \code{seed}: Numeric. The random seed used for reproducibility.
#' }
#' This function is not a primary user-facing function of the package and
#' does not include examples.
#'
#' @import dplyr
#' @importFrom future plan multisession
#' @importFrom future.apply future_lapply
#' @importFrom rstan rstan_options
#' @importFrom parallel detectCores
#' @importFrom tibble tibble
#' @export
Carlotti_and_Parast_2026_simulations <- function(seed, n_simulations, parallel = TRUE) {
  ###################################
  # Random seed for reproducibility #
  ###################################
  set.seed(seed)
  
  ###########################
  # Data generating process #
  ###########################
  
  # Sample size
  n <- 50
  
  # Treatment assignment probability
  p <- 0.5
  
  # Possible settings
  n_settings <- 2
  settings <- vector("list", n_settings)
  
  for (setting in 1:n_settings) {
    
    # Choice of the setting
    if (setting == 1) {
      ############
      # X binary #
      ############
      
      # Setting name
      setting_name <- "X_binary"
      
      # Setting label
      setting_label <- "X binary"
      
      # Binary covariate probability
      q <- 0.5
      
      # Number of covariates (including the intercept)
      d <- 2
      
      # Probability of X = 1
      q <- 0.5
      
      # Mean vector for potential outcomes when X = 0
      mu_0 <- c(5, 5, 0, 0)
      
      # Mean vector for potential outcomes when X = 1
      mu_1 <- c(5, -5, 0, -10)
      
      # Covariance matrix for potential outcomes when X = 0
      Sigma_0 <- kronecker(diag(2), matrix(
        data = c(1, 1,
                 1, 2),
        nrow = 2,
        ncol = 2))
      
      # Covariance matrix for potential outcomes when X = 1
      Sigma_1 <- kronecker(diag(2), matrix(
        data = c(1, 1,
                 1, 2),
        nrow = 2,
        ncol = 2))
      
      # Store setting parameters in the list
      settings[[setting]] <- list(
        setting_name = setting_name,
        setting_label = setting_label,
        q = q,
        d = d,
        mu_0 = mu_0,
        mu_1 = mu_1,
        Sigma_0 = Sigma_0,
        Sigma_1 = Sigma_1
      )
    } else if (setting == 2) {
      ##############
      # X Gaussian #
      ##############
      
      # Setting name
      setting_name <- "X_gaussian"
      
      # Setting label
      setting_label <- "X Gaussian"
      
      # Number of covariates (including the intercept)
      d <- 2
      
      # Vector of coefficients for the covariates in the potential outcomes model
      beta <- c(1, 7, 0, 6)
      
      # Covariance matrix of the potential outcomes
      Sigma <- 0.5 * diag(4)
      
      # Mean vector the Gaussian covariate X
      m <- 3
      
      # Standard deviation of the Gaussian covariate X
      s <- 1
      
      # Store setting parameters in the list
      settings[[setting]] <- list(
        setting_name = setting_name,
        setting_label = setting_label,
        d = d,
        beta = beta,
        Sigma = Sigma,
        m = m,
        s = s
      )
    }
  }
  
  ######################
  # Testing parameters #
  ######################
  
  # Value of epsilon under the null hypothesis for the Bayesian test
  V_S_zero <- 0.5
  
  # Parameters for the Beta distribution used in the Bayesian test
  a <- 1
  b <- 1
  
  # Bayes factor alternative hypothesis
  BF_alternative <- "greater"
  
  # Probability of Type I error
  alpha <- 0.05
  
  # Probability of Type II error
  beta <- 0.2
  
  ####################
  # Prior parameters #
  ####################
  
  # Prior parameters for the prior distributions
  s <- rep(2, 4)
  tau <- 1
  
  ######################
  # Posterior sampling #
  ######################
  
  # Number of MCMC chains
  n_chains <- 2
  
  # Number of MCMC iterations per chain
  n_iter <- 500
  
  ###################
  # Run simulations #
  ###################
  
  # Grid of simulations
  simulations_grid <- base::expand.grid(
    setting = 1:n_settings,
    simulation = 1:n_simulations
  )
  
  # Function to run a single simulation for Parast et al. (2024) settings
  Carlotti_and_Parast_2026_single_simulation <- function(simulation_id, simulation_seed) {
    # Set the seed for reproducibility
    if (!is.null(simulation_seed)) {
      set.seed(simulation_seed)
    }
    
    # Extract setting and simulation number
    simulation_setting <- simulation_id$setting
    simulation_number <- simulation_id$simulation
    
    # Generate data
    if (simulation_setting == 1) {
      data <- DGP_X_binary(
        n = n,
        p = p,
        q = settings[[simulation_setting]]$q,
        mu_0 = settings[[simulation_setting]]$mu_0,
        mu_1 = settings[[simulation_setting]]$mu_1,
        Sigma_0 = settings[[simulation_setting]]$Sigma_0,
        Sigma_1 = settings[[simulation_setting]]$Sigma_1
      )
    } else if (simulation_setting == 2) {
      data <- DGP_X_Gaussian(
        n = n,
        p = p,
        beta = settings[[simulation_setting]]$beta,
        Sigma = settings[[simulation_setting]]$Sigma,
        m = settings[[simulation_setting]]$m,
        s = settings[[simulation_setting]]$s
      )
    }
    
    observed_data <- data.frame(
      Y = data$P_observed[, "Y"],
      S = data$P_observed[, "S"],
      Z = data$Z,
      X = data$X
    )
    
    # Prior parameters for the prior distributions
    mu_beta <- array(0, dim = c(settings[[simulation_setting]]$d))
    Sigma_beta <- 10*diag(settings[[simulation_setting]]$d)
    
    BSET_results <- BSET(
      data = observed_data,
      Y = "Y",
      S = "S",
      Z = "Z",
      X = "X",
      delta_true = estimands_Carlotti_and_Parast_2026$delta_MC[simulation_setting],
      theta_true = estimands_Carlotti_and_Parast_2026$theta_MC[simulation_setting],
      seed = seed + simulation_seed,
      n_chains = n_chains,
      n_iter = n_iter,
      burn_in_ratio = 0.25,
      a = a,
      b = b,
      alpha = alpha,
      beta = beta,
      V_S_zero = V_S_zero,
      BF_alternative = BF_alternative,
      root_tolerance = 1e-16,
      mu_beta = mu_beta,
      Sigma_beta = Sigma_beta,
      s = s,
      tau = tau,
      plot = TRUE,
      mute = TRUE,
      parallel = FALSE
    )

    return(BSET_results)
  }
  
  # Run all simulations
  if (parallel) {
    # Detect available cores
    available_cores <- parallel::detectCores() - 1
    
    # Maximum allowed PSOCK connections. R default is 128
    max_connections <- getOption("future.connections.max", 128) - 100
    
    # Number of cores to use
    workers <- min(available_cores, max_connections)
    
    # Plan parallel execution
    future::plan(future::multisession, workers = workers)
    
    # Run all simulations in parallel
    Carlotti_and_Parast_2026_all_simulations <- future.apply::future_lapply(
      X = 1:(n_settings*n_simulations),
      FUN = function(i) {
        Carlotti_and_Parast_2026_single_simulation(
          simulation_id = simulations_grid[i, ],
          simulation_seed = i
        )
      },
      future.seed = TRUE
    )
    
    # Return to sequential processing for any subsequent steps
    future::plan(future::sequential)
  } else {
    # Run all simulations sequentially
    Carlotti_and_Parast_2026_all_simulations <- lapply(
      X = 1:(n_settings*n_simulations),
      FUN = function(i) {
        Carlotti_and_Parast_2026_single_simulation(
          simulation_id = simulations_grid[i, ],
          simulation_seed = i
        )
      }
    )
  }
  
  # Remove names from the list of results for easier processing
  names(Carlotti_and_Parast_2026_all_simulations) <- NULL
  
  # Extract results into a dataframe
  Carlotti_and_Parast_2026_all_simulations_df <- dplyr::bind_rows(lapply(Carlotti_and_Parast_2026_all_simulations, function(res) {
    tibble::tibble(
      Bayesian_epsilon = res$Bayesian_test$epsilon,
      Bayesian_CI_upper = res$Bayesian_test$CI[2],
      Bayesian_coverage = as.logical(res$Bayesian_test$coverage),
      Bayesian_power = as.logical(res$Bayesian_test$power),
      frequentist_epsilon = res$frequentist_test$epsilon,
      frequentist_CI_upper = res$frequentist_test$CI[2],
      frequentist_coverage = as.logical(res$frequentist_test$coverage),
      frequentist_power = as.logical(res$frequentist_test$power)
    )
  }))

  Carlotti_and_Parast_2026_simulations_grid <-
    bind_cols(
      simulations_grid,
      Carlotti_and_Parast_2026_all_simulations_df
    )
  
  Carlotti_and_Parast_2026_simulations_grid <-
    Carlotti_and_Parast_2026_simulations_grid %>%
    mutate(
      n_sample = n,
      n_chains = n_chains,
      n_iterations = n_iter,
      n_simulations = n_simulations,
      BF_alternative = BF_alternative,
      timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      seed = seed
    )
  
  ##########
  # Output #
  ##########
  
  return(Carlotti_and_Parast_2026_simulations_grid)
}
