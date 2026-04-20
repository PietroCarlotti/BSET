#' Monte Carlo Computation of the Estimands for the Simulation Study in Parast et al. (2024)
#'
#' This function iterates through the four simulation settings defined in Parast
#' et al. (2024) and estimates the true values of \eqn{U_Y}, \eqn{U_S},
#' \eqn{\delta}, \eqn{V_Y}, \eqn{V_S}, and \eqn{\theta} using a Monte Carlo
#' dataset generated according to the specified data-generating processes. 
#' 
#' @details
#' The settings are defined as follows:
#' \itemize{
#'   \item Setting 1: \strong{Useless surrogate (Gaussian model)}
#'     \deqn{Y_1 \sim \mathcal{N}(5, 3), \quad Y_0 \sim \mathcal{N}(3, 3),}
#'     \deqn{S_1 \sim \mathcal{N}(2/3, 1), \quad S_0 \sim \mathcal{N}(2/3, 1).}
#'   \item Setting 2: \strong{Perfect surrogate (Gaussian model)}
#'     \deqn{Y_1 \sim \mathcal{N}(6, 3), \quad Y_0 \sim \mathcal{N}(5/2, 3),}
#'     \deqn{S_1 = Y_1 + \mathcal{N}(0, 1/10), \quad S_0 = Y_0 + \mathcal{N}(0, 1/10).}
#'   \item Setting 3: \strong{Imperfect surrogate (Gaussian model)}
#'     \deqn{S_1 \sim \mathcal{N}(5, 3), \quad S_0 \sim \mathcal{N}(3, 3),}
#'     \deqn{Y_1 = S_1 + \mathcal{N}(1.5, 0.6), \quad Y_0 = S_0 + \mathcal{N}(0, 0.6).}
#'   \item Setting 4: \strong{Misspecified model (non-Gaussian model)}
#'     \deqn{S_1 \sim \exp(\mathcal{N}(2.5, 1.5)), \quad S_0 \sim \exp(\mathcal{N}(0.5, 1.5)),}
#'     \deqn{Y_1 = 2 + 6/5 \sqrt{S_1} + 3/10 \exp(S_1 / 500) + \exp(\mathcal{N}(0, 0.3)),}
#'     \deqn{Y_0 = 4/5 \sqrt{S_0} + 1/5 \exp(S_0 / 50) + \exp(\mathcal{N}(0, 0.3)).}
#' }
#' This function is not a primary user-facing function of the package and
#' does not include examples.
#'
#' @param MC_samples Integer. The number of Monte Carlo samples to generate per setting.
#'
#' @return A data frame containing the Monte Carlo estimates for each setting:
#' \itemize{
#'   \item \code{setting}: The index of the simulation setting.
#'   \item \code{U_Y_MC}, \code{U_S_MC}, \code{delta_MC}: Parameters of interest from Parast et al. (2024).
#'   \item \code{V_Y_MC}, \code{V_S_MC}, \code{theta_MC}: Parameters of interest from Carlotti and Parast (2026).
#' }
#'
#' @export
compute_estimands_Parast_et_al_2024 <- function(MC_samples) {
  # Treatment assignment probability
  p = 0.5
  
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
  
  # Initialize estimands dataframe
  estimands <- data.frame(
    setting = 1:length(settings),
    U_Y_MC = as.numeric(NA),
    U_S_MC = as.numeric(NA),
    delta_MC = as.numeric(NA),
    V_Y_MC = as.numeric(NA),
    V_S_MC = as.numeric(NA),
    theta_MC = as.numeric(NA)
  )
  
  # Loop over settings
  for (i in 1:n_settings) {
    # Extract setting parameters
    model_type <- settings[[i]]$model
    mu_vals    <- settings[[i]]$mu_star
    sigma_vals <- settings[[i]]$Sigma_star
    
    # Generate Monte Carlo data
    MC_data <- DGP_no_X(
      n = MC_samples,
      p = p,
      mu_star = mu_vals,
      Sigma_star = sigma_vals,
      model = model_type
    )
    
    # Compute estimands
    delta_vals <- compute_delta(MC_data)
    theta_vals <- compute_theta(MC_data)
      
    # Store results
    estimands$U_Y_MC[i] <- delta_vals$U_Y
    estimands$U_S_MC[i] <- delta_vals$U_S
    estimands$delta_MC[i] <- delta_vals$delta
    estimands$V_Y_MC[i] <- theta_vals$V_Y
    estimands$V_S_MC[i] <- theta_vals$V_S
    estimands$theta_MC[i] <- theta_vals$theta
  }

  return(estimands)
}
