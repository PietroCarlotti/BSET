#' Monte Carlo Computation of the Estimands for the Simulation Study in Carlotti and Parast (2026)
#'
#' This function iterates through the simulation settings defined in Carlotti
#' and Parast (2026) and estimates the true values of \eqn{U_Y}, \eqn{U_S},
#' \eqn{\delta}, \eqn{V_Y}, \eqn{V_S}, and \eqn{\theta} using a Monte Carlo
#' dataset generated according to the specified data-generating processes. The
#' settings are defined as follows:
#' \itemize{
#'   \item Setting 1: \strong{X binary}
#'   \itemize{
#'     \item If \eqn{X = 1}:
#'     \deqn{Y_1 \sim \mathcal{N}(5, 1), \quad Y_0 \sim \mathcal{N}(0, 1)}
#'     \deqn{S_1 = Y_1 + \mathcal{N}(0, 1), \quad S_0 = Y_0 + \mathcal{N}(0, 1)}
#'     \item If \eqn{X = 0}:
#'     \deqn{Y_1 \sim \mathcal{N}(5, 1), \quad Y_0 \sim \mathcal{N}(0, 1)}
#'     \deqn{S_1 = Y_1 + \mathcal{N}(-10, 1), \quad S_0 = Y_0 + \mathcal{N}(-10, 1)}
#'   }
#'   \item Setting 2: \strong{X Gaussian}
#'   \deqn{()}
#' }
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
compute_estimands_Carlotti_and_Parast_2026 <- function(MC_samples) {
  # Treatment assignment probability
  p = 0.5
  
  # Possible settings
  n_settings <- 1
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
    setting_name <- settings[[i]]$setting_name
    setting_label <- settings[[i]]$setting_label
    d <- settings[[i]]$d
    
    if (setting_name == "X_binary") {
      # Binary covariate parameters
      q <- settings[[i]]$q
      mu_0 <- settings[[i]]$mu_0
      mu_1 <- settings[[i]]$mu_1
      Sigma_0 <- settings[[i]]$Sigma_0
      Sigma_1 <- settings[[i]]$Sigma_1
      
      # Generate large-scale Monte Carlo data
      MC_data <- DGP_X_binary(
        n = MC_samples,
        p = p,
        q = q,
        mu_0 = mu_0,
        mu_1 = mu_1,
        Sigma_0 = Sigma_0,
        Sigma_1 = Sigma_1
      )
    } else if (setting_name == "X_gaussian") {
      # Gaussian covariate parameters
      beta <- settings[[i]]$beta
      Sigma <- settings[[i]]$Sigma
      m <- settings[[i]]$m
      s <- settings[[i]]$s
      
      # Generate large-scale Monte Carlo data
      MC_data <- DGP_X_Gaussian(
        n = MC_samples,
        p = p,
        beta = beta,
        Sigma = Sigma,
        m = m,
        s = s
      )
    }
    
    # Compute estimands
    delta_estimands <- compute_delta(MC_data)
    theta_estimands <- compute_theta(MC_data)
    
    # Store results
    estimands$U_Y_MC[i] <- delta_estimands$U_Y
    estimands$U_S_MC[i] <- delta_estimands$U_S
    estimands$delta_MC[i] <- delta_estimands$delta
    estimands$V_Y_MC[i] <- theta_estimands$V_Y
    estimands$V_S_MC[i] <- theta_estimands$V_S
    estimands$theta_MC[i] <- theta_estimands$theta
  }
  
  return(estimands)
}
