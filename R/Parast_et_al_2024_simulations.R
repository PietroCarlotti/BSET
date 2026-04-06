#' Run Simulations for the Settings in Parast et al. (2024)
#'
#' This function executes the simulation study described in Parast et al. (2024)
#' across four different settings:
#' \enumerate{
#'  \item Useless surrogate.
#'  \item Perfect surrogate.
#'  \item Imperfect surrogate.
#'   \item Misspecified model.
#'}
#' It performs the frequentist and Bayesian surrogate evaluation tests from
#' Parast et al. (2024) and Carlotti and Parast (2026), respectively, and
#' returns a structured data frame of the results. The execution is parallelized
#' using the \code{future} framework.
#'
#' @details
#' The function runs a total of 500 simulations for each of the four settings,
#' generating data using \code{DGP_no_X} and processes the results using
#' \code{BSET_no_X}. The simulation utilizes MCMC sampling via \code{rstan} for
#' the Bayesian estimation components. Note that it relies on an external object 
#' \code{estimands_Parast_et_al_2024} for true parameter values, which are
#' computed using the function \code{compute_estimands_Parast_et_al_2024}.
#' 
#' @param seed Numeric. A random seed for reproducibility of the simulations.
#' @param n_simulations Numeric. The number of simulations to run for each setting.
#' @param parallel Logical. Whether to run the simulations in parallel or sequentially (default is TRUE).
#' @return A data frame containing:
#' \itemize{
#'   \item \code{setting}: The index of the simulation setting (1 to 4).
#'   \item \code{simulation}: The individual simulation run ID.
#'   \item \code{Bayesian_epsilon}, \code{frequentist_epsilon}: Numeric. Surrogate validation thresholds for the Bayesian and frequentist tests, respectively.
#'   \item \code{Bayesian_coverage}, \code{frequentist_coverage}: Boolean. Indicates whether the true surrogate effect is within the credible or confidence interval for the Bayesian and frequentist tests, respectively.
#'   \item \code{Bayesian_power}, \code{frequentist_power}: Boolean. Indicates whether the upper bound of the credible or confidence interval is below the threshold for the Bayesian and frequentist tests, respectively.
#'   \item \code{n_sample}, \code{n_simulations}, \code{n_chains}, \code{n_iterations}: Numeric. Metadata about the simulation parameters.
#'   \item \code{timestamp}: Character. Date and time when the simulation was completed.
#'   \item \code{seed}: Numeric. The random seed used for reproducibility.
#' }
#' 
#' @import dplyr
#' @importFrom purrr map_dfr
#' @importFrom future plan multisession
#' @importFrom future.apply future_lapply
#' @importFrom rstan rstan_options
#' @importFrom parallel detectCores
#' @importFrom tibble tibble
#' @export
Parast_et_al_2024_simulations <- function(seed, n_simulations, parallel) {
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
  n_settings <- 4
  settings <- vector("list", n_settings)
  
  for (setting in 1:n_settings) {
    
    # Choice of the setting
    if (setting == 1) {
      
      #####################
      # Useless surrogate #
      #####################
      
      # Setting name
      setting_name <- "useless_surrogate"
      
      # Setting label
      setting_label <- "Useless Surrogate"
      
      # Model type
      model <- "Gaussian"
      
      # True parameters
      mu_star <- c(5, 2/3, 3, 2/3)
      Sigma_star <- diag(c(3, 1, 3, 1))
      
    } else if (setting == 2) {
      
      #####################
      # Perfect surrogate #
      #####################
      
      # Setting name
      setting_name <- "perfect_surrogate"
      
      # Setting label
      setting_label <- "Perfect Surrogate"
      
      # Model type
      model <- "Gaussian"
      
      # True parameters for Y and epsilon
      mu_Y <- c(6, 2.5)
      Sigma_Y <- 3*diag(2)
      
      mu_epsilon <- rep(0, 2)
      Sigma_epsilon <- 0.1*diag(2)
      
      mu_Y_epsilon <- c(mu_Y, mu_epsilon)
      Sigma_Y_epsilon <- rbind(
        cbind(Sigma_Y, matrix(0, 2, 2)),
        cbind(matrix(0, 2, 2), Sigma_epsilon)
      )
      
      # Mixing matrix
      A <- matrix(
        data = c(1, 0, 0, 0,
                 1, 0, 1, 0,
                 0, 1, 0, 0,
                 0, 1, 0, 1),
        nrow = 4,
        ncol = 4,
        byrow = TRUE
      )
      
      # True parameters
      mu_star <- A %*% mu_Y_epsilon
      Sigma_star <- A %*% Sigma_Y_epsilon %*% t(A)
      
    } else if (setting == 3) {
      
      #######################
      # Imperfect surrogate #
      #######################
      
      # Setting name
      setting_name <- "imperfect_surrogate"
      
      # Setting label
      setting_label <- "Imperfect Surrogate"
      
      # Model type
      model <- "Gaussian"
      
      # True parameters for S and epsilon
      mu_S <- c(5, 3)
      Sigma_S <- 3*diag(2)
      
      mu_epsilon <- c(1.5, 0)
      Sigma_epsilon <- 0.6*diag(2)
      
      mu_S_epsilon <- c(mu_S, mu_epsilon)
      Sigma_S_epsilon <- rbind(
        cbind(Sigma_S, matrix(0, 2, 2)),
        cbind(matrix(0, 2, 2), Sigma_epsilon)
      )
      
      # Mixing matrix
      A <- matrix(
        data = c(1, 0, 1, 0,
                 1, 0, 0, 0,
                 0, 1, 0, 1,
                 0, 1, 0, 0),
        nrow = 4,
        ncol = 4,
        byrow = TRUE
      )
      
      # True parameters
      mu_star <- A %*% mu_S_epsilon
      Sigma_star <- A %*% Sigma_S_epsilon %*% t(A)
      
    } else if (setting == 4) {
      
      ######################
      # Misspecified model #
      ######################
      
      # Setting name
      setting_name <- "misspecified"
      
      # Setting label
      setting_label <- "Misspecified Model"
      
      # Model type
      model <- "misspecified"
      
      # True parameters
      mu_star <- NULL
      Sigma_star <- NULL
      
    }
    
    settings[[setting]] <- list(
      setting_name = setting_name,
      setting_label = setting_label,
      model = model,
      mu_star = mu_star,
      Sigma_star = Sigma_star
    )
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
  mu_0 <- rep(0, 4)
  Sigma_0 <- 10*diag(4)
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
  simulations_grid <- expand.grid(
    setting = 1:n_settings,
    simulation = 1:n_simulations
  )
  
  # Ensure compiled models are saved and reused
  rstan_options(auto_write = TRUE)
  
  # Function to run a single simulation for Parast et al. (2024) settings
  Parast_et_al_2024_single_simulation <- function(simulation_id, simulation_seed) {
    # Set the seed for reproducibility
    if (!is.null(simulation_seed)) {
      set.seed(simulation_seed)
    }
    
    # Extract setting and simulation number
    simulation_setting <- simulation_id$setting
    simulation_number <- simulation_id$simulation
    
    # PRINT PROGRESS (DELETE AFTER YOU ARE DONE DEBUGGING)
    message(paste0("Running simulation ", simulation_number, " for setting ", simulation_setting, " (", settings[[simulation_setting]]$setting_label, ")"))
    
    # Generate data
    data <- DGP_no_X(
      n = n,
      p = p,
      model = settings[[simulation_setting]]$model,
      mu_star = settings[[simulation_setting]]$mu_star,
      Sigma_star = settings[[simulation_setting]]$Sigma_star
    )
    
    observed_data <- data.frame(
      Y = data$P_observed[, "Y"],
      S = data$P_observed[, "S"],
      Z = data$Z
    )
    
    BSET_no_X_results <- BSET_no_X(
      data = observed_data,
      Y = "Y",
      S = "S",
      Z = "Z",
      delta_true = estimands_Parast_et_al_2024$delta_MC[simulation_setting],
      theta_true = estimands_Parast_et_al_2024$theta_MC[simulation_setting],
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
      mu_0 = mu_0,
      Sigma_0 = Sigma_0,
      s = s,
      tau = tau,
      plot = TRUE,
      mute = TRUE,
      parallel = FALSE
    )
    
    return(BSET_no_X_results)
  }
  
  # Run all simulations
  if (parallel) {
    # Detect available cores
    available_cores <- parallel::detectCores() - 1
    
    # Maximum allowed PSOCK connections. R default is 128
    max_connections <- getOption("future.connections.max", 128) - 4
    
    # Number of cores to use
    workers <- min(available_cores, max_connections)
    
    # Plan parallel execution
    future::plan(multisession, workers = workers)
    
    # Run all simulations in parallel
    Parast_et_al_2024_all_simulations <- future_lapply(
      X = 1:(n_settings*n_simulations),
      FUN = function(i) {
        Parast_et_al_2024_single_simulation(
          simulation_id = simulations_grid[i, ],
          simulation_seed = i
        )
      },
      future.seed = TRUE
    )
    
    # Close workers after finishing
    future::plan(future::sequential)
  } else {
    # Run all simulations sequentially
    Parast_et_al_2024_all_simulations <- lapply(
      X = 1:(n_settings*n_simulations),
      FUN = function(i) {
        Parast_et_al_2024_single_simulation(
          simulation_id = simulations_grid[i, ],
          simulation_seed = i
        )
      }
    )
  }
  
  # Remove names from the list of results for easier processing
  names(Parast_et_al_2024_all_simulations) <- NULL
  
  # Extract results into a dataframe
  Parast_et_al_2024_all_simulations_df <- map_dfr(Parast_et_al_2024_all_simulations, function(res) {
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
  })
  
  Parast_et_al_2024_simulations_grid <-
    bind_cols(
      simulations_grid,
      Parast_et_al_2024_all_simulations_df
    )
  
  Parast_et_al_2024_simulations_grid <-
    Parast_et_al_2024_simulations_grid %>%
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
  
  return(Parast_et_al_2024_simulations_grid)
}
