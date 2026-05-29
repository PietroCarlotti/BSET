## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----install-BSET-------------------------------------------------------------
# install.packages("pak")
# pak::pak("PietroCarlotti/BSET")

## ----packages, message=FALSE, warning=FALSE-----------------------------------
library(ggplot2)
library(dplyr)

## ----quick-start-data---------------------------------------------------------
# Set the random seed for reproducibility
set.seed(1234)

# Sample size
n <- 50

# Binary covariate
X <- rbinom(n, 1, 0.5)

# Treatment assignment
Z <- rbinom(n, 1, 0.5)

# Primary outcome
Y <- 2 * Z -3 * X + rnorm(n)

# Surrogate outcome
S <- Y + rnorm(n)

df <- data.frame(
  Y = Y,
  S = S,
  Z = Z,
  X = X
)

## ----quick-start-no-X, eval=FALSE---------------------------------------------
# result_no_X <- BSET::BSET(
#   data = df,
#   Y = "Y",
#   S = "S",
#   Z = "Z",
#   seed = 123,
#   plot = TRUE
# )
# result_no_X$theta_posterior_plot

## ----quick-start-no-X-show, echo=FALSE, out.width="100%", fig.align='center', fig.cap="Posterior distribution of $\\theta$ from `BSET` without covariates. Dashed lines show estimated quantities: the Bayesian 95\\% credible interval upper bound and threshold $\\eta$ (<span style='color:#0072B2'>blue shades</span>), and the frequentist 95\\% confidence interval upper bound and threshold $\\varepsilon$ (<span style='color:#D55E00'>orange shades</span>)."----
knitr::include_graphics("BSET_quick_start_no_X_plot.png")

## ----quick-start-X, eval=FALSE------------------------------------------------
# result_X <- BSET::BSET(
#   data = df,
#   Y = "Y",
#   S = "S",
#   Z = "Z",
#   X = "X",
#   seed = 123,
#   plot = TRUE
# )
# result_X$theta_posterior_plot

## ----quick-start-X-show, echo=FALSE, out.width="100%", fig.align='center', fig.cap="Posterior distribution of $\\theta$ from `BSET`, adjusted for a baseline covariate. Dashed lines show estimated quantities: the Bayesian 95\\% credible interval upper bound and threshold $\\eta$ (<span style='color:#0072B2'>blue shades</span>), and the frequentist 95\\% confidence interval upper bound and threshold $\\varepsilon$ (<span style='color:#D55E00'>orange shades</span>)."----
knitr::include_graphics("BSET_quick_start_X_plot.png")

## ----set-seed-----------------------------------------------------------------
# Set the random seed for reproducibility
set.seed(1234)

# Sample size and treatment assignment probability used throughout
n <- 50
p <- 0.5

## ----generate-data-no-X-------------------------------------------------------
# Mean vector of potential outcomes for the primary outcome and the surrogate
mu_star <- c(6, 6, 2.5, 2.5)

# Covariance matrix of potential outcomes for the primary outcome and the surrogate
Sigma_star <- kronecker(
  diag(2),
  matrix(
    data = c(3, 3,
             3, 3.1),
    nrow = 2,
    ncol = 2)
)

# Generate the data
data_no_X <- BSET::DGP_no_X(
  n = n,
  p = p,
  mu_star = mu_star,
  Sigma_star = Sigma_star,
  model = "Gaussian"
)

## ----generate-data-X----------------------------------------------------------
# Binary covariate probability
q <- 0.5

# Mean vector for potential outcomes when X = 0 and when X = 1
mu_0 <- c(5, 5, 0, 0)
mu_1 <- c(5, -5, 0, -10)

# Covariance matrix for potential outcomes when X = 0 and when X = 1
Sigma_0 <- kronecker(
  diag(2),
  matrix(
    data = c(1, 1,
             1, 2),
    nrow = 2,
    ncol = 2)
  )

Sigma_1 <- kronecker(
  diag(2),
  matrix(
    data = c(1, 1,
             1, 2),
    nrow = 2,
    ncol = 2)
  )

# Generate the data
data_X <- BSET::DGP_X_binary(
  n = n,
  p = p,
  q = q,
  mu_0 = mu_0,
  mu_1 = mu_1,
  Sigma_0 = Sigma_0,
  Sigma_1 = Sigma_1
  )

## ----compute-delta-theta-no-X-------------------------------------------------
# Compute delta
delta_no_X <- BSET::compute_delta(MC_data = data_no_X)

# Compute theta
theta_no_X <- BSET::compute_theta(MC_data = data_no_X)

## ----compute-delta-theta-X----------------------------------------------------
# Compute delta
delta_X <- BSET::compute_delta(MC_data = data_X)

# Compute theta
theta_X <- BSET::compute_theta(MC_data = data_X)

## ----delta-theta-table, echo=FALSE--------------------------------------------
delta_theta_table <- data.frame(
  Setting = c("No X (Perfect Surrogate)", "Binary X"),
  "$\\widehat{\\delta}$" = c(delta_no_X$delta, delta_X$delta),
  "$\\widehat{\\theta}$" = c(theta_no_X$theta, theta_X$theta),
  check.names = FALSE 
)

knitr::kable(
  delta_theta_table, 
  caption = "Estimated values of $\\delta$ and $\\theta$ in the two settings.",
  escape = FALSE,
  booktabs = TRUE,
  align = "ccc",
  digits = 3
)

## ----compute-delta-theta-true-code, eval=FALSE--------------------------------
# # Compute Monte Carlo estimands (based on 1,000,000 samples) — this is slow
# estimands_Parast_et_al_2024 <- BSET::compute_estimands_Parast_et_al_2024(n_MC = 1000000)
# estimands_Carlotti_and_Parast_2026 <- BSET::compute_estimands_Carlotti_and_Parast_2026(n_MC = 1000000)

## ----compute-delta-theta-true-------------------------------------------------
# Load precomputed Monte Carlo estimands (based on 1,000,000 samples)
estimands_Parast_et_al_2024 <- BSET::estimands_Parast_et_al_2024
estimands_Carlotti_and_Parast_2026 <- BSET::estimands_Carlotti_and_Parast_2026

## ----estimands-Parast-et-al-2024, echo=FALSE----------------------------------
estimands_Parast_et_al_2024_table <- data.frame(
  Setting = c("Useless surrogate", "Perfect surrogate", "Imperfect surrogate", "Misspecified model"),
  "$\\delta$" = estimands_Parast_et_al_2024$delta_MC,
  "$\\theta$" = estimands_Parast_et_al_2024$theta_MC,
  check.names = FALSE
)

knitr::kable(
  estimands_Parast_et_al_2024_table,
  caption = "Monte Carlo estimates of $\\delta$ and $\\theta$ for the simulation settings considered in Parast et al. (2024).",
  booktabs = TRUE,
  align = "ccc",
  digits = 3
)

## ----estimands-Carlotti-and-Parast-2026, echo=FALSE---------------------------
estimands_Carlotti_and_Parast_2026_table <- data.frame(
  Setting = c("Perfect surrogate"),
  "$\\delta$" = estimands_Carlotti_and_Parast_2026$delta_MC,
  "$\\theta$" = estimands_Carlotti_and_Parast_2026$theta_MC,
  check.names = FALSE
)

knitr::kable(
  estimands_Carlotti_and_Parast_2026_table,
  caption = "Monte Carlo estimates of $\\delta$ and $\\theta$ for the simulation settings considered in Carlotti and Parast (2026).",
  booktabs = TRUE,
  align = "ccc",
  digits = 3
)

## ----compute-BF-distribution--------------------------------------------------
# Hypothesized value of V_S under the null
V_S_zero <- 0.5

# Prior parameters for the Beta distribution
a <- 1
b <- 1

# Alternative hypothesis for the Bayes factor
BF_alternative <- "greater"

# Compute the distribution of the Bayes factor
BF_distribution <- BSET::compute_BF_distribution(
  n = n,
  V_S_true = V_S_zero,
  V_S_zero = V_S_zero,
  a = a,
  b = b,
  BF_alternative = BF_alternative
)

## ----BF-distribution-plot, echo=FALSE, fig.width=8, fig.height=4.5, fig.align='center', fig.cap="Distribution of the Bayes factor under the null hypothesis $V_S = 0.5$."----
# Log transform the Bayes factor
BF_distribution$log_BF_values <- log(BF_distribution$BF_values)

ggplot(BF_distribution, aes(x = log_BF_values, y = BF_PMF)) +
  geom_col(fill = "steelblue", width = 0.07) +
  labs(
    x = expression(log(BF[n])),
    y = "PMF"
  ) +
  theme_minimal()

## ----compute-BF-alpha---------------------------------------------------------
# Type I error rate
alpha <- 0.05

# Compute the critical value of the Bayes factor corresponding to alpha
BF_alpha <- BF_distribution %>%
    dplyr::filter(BF_CDF >= 1 - alpha) %>%
    dplyr::slice(1) %>%
    dplyr::pull(BF_values)

## ----compute-V_S-star---------------------------------------------------------
# Type II error rate
beta <- 0.2

# Compute the value of v_S that satisfies the power constraint
V_S_star <- BSET::compute_V_S_star(
  n = n,
  alpha = alpha,
  beta = beta,
  V_S_zero = V_S_zero,
  a = a,
  b = b,
  BF_alternative = BF_alternative,
  root_tolerance = 1e-16
)$V_S_star

## ----compute-eta--------------------------------------------------------------
# Hypothesized value of the treatment effect on the primary outcome
v_Y <- estimands_Parast_et_al_2024$V_Y_MC[2]

# Compute the validation threshold eta
eta <- max(v_Y - V_S_star, 0)

## ----BSET-no-X-data-----------------------------------------------------------
# Prepare the data for BSET (no covariates)
BSET_no_X_data <- data.frame(
  Y = data_no_X$P_observed[, "Y"],
  S = data_no_X$P_observed[, "S"],
  Z = data_no_X$Z 
)

## ----BSET-no-X-posterior-sampling-parameters----------------------------------
# Posterior sampling parameters
n_chains <- 2
n_iter <- 2000
burn_in_ratio <- 0.25

## ----BSET-no-X-prior-parameters-----------------------------------------------
# Prior parameters
mu_0 <- rep(0, 4)
Sigma_0 <- diag(4)
s <- rep(1, 4)
tau <- 1

## ----BSET-no-X, eval=FALSE----------------------------------------------------
# # Run the BSET procedure without adjusting for covariates
# BSET_no_X_results <- BSET::BSET(
#   data = BSET_no_X_data,
#   Y = "Y",
#   S = "S",
#   Z = "Z",
#   delta_true = estimands_Parast_et_al_2024$delta_MC[2],
#   theta_true = estimands_Parast_et_al_2024$theta_MC[2],
#   seed = 123,
#   n_chains = n_chains,
#   n_iter = n_iter,
#   burn_in_ratio = burn_in_ratio,
#   a = a,
#   b = b,
#   alpha = alpha,
#   beta = beta,
#   V_S_zero = V_S_zero,
#   BF_alternative = BF_alternative,
#   root_tolerance = 1e-16,
#   mu_0 = mu_0,
#   Sigma_0 = Sigma_0,
#   s = s,
#   tau = tau,
#   plot = TRUE,
#   mute = TRUE,
#   parallel = TRUE
# )

## ----BSET-no-X-load, include=FALSE--------------------------------------------

## ----BSET-no-X-plot, echo=FALSE, out.width="100%", fig.align='center', fig.cap="Posterior distribution of $\\theta$ from the BSET procedure without adjusting for covariates. Dashed lines show estimated quantities: the Bayesian 95\\% credible interval upper bound and threshold $\\eta$ (<span style='color:#0072B2'>blue shades</span>), and the frequentist 95\\% confidence interval upper bound and threshold $\\varepsilon$ (<span style='color:#D55E00'>orange shades</span>). Solid lines show the true values of $\\theta$ (<span style='color:#551A8B'>purple</span>) and $\\delta$ (<span style='color:#CC79A7'>pink</span>)."----
knitr::include_graphics("BSET_no_X_plot.png")

## ----BSET-X-data--------------------------------------------------------------
# Prepare the data for BSET (with covariates)
BSET_X_data <- data.frame(
  Y = data_X$P_observed[, "Y"],
  S = data_X$P_observed[, "S"],
  Z = data_X$Z,
  X = data_X$X
)

## ----BSET-X-prior-parameters--------------------------------------------------
# Prior parameters for the regression coefficients of the potential outcomes on the covariates
d <- 2
mu_beta <- rep(0, d)
Sigma_beta <- 10*diag(d)

## ----BSET-X, eval=FALSE-------------------------------------------------------
# # Run the BSET procedure adjusting for covariates
# BSET_X_results <- BSET::BSET(
#   data = BSET_X_data,
#   Y = "Y",
#   S = "S",
#   Z = "Z",
#   X = "X",
#   delta_true = estimands_Carlotti_and_Parast_2026$delta_MC[1],
#   theta_true = estimands_Carlotti_and_Parast_2026$theta_MC[1],
#   seed = 123,
#   n_chains = n_chains,
#   n_iter = n_iter,
#   burn_in_ratio = burn_in_ratio,
#   a = a,
#   b = b,
#   alpha = alpha,
#   beta = beta,
#   V_S_zero = V_S_zero,
#   BF_alternative = BF_alternative,
#   root_tolerance = 1e-16,
#   mu_beta = mu_beta,
#   Sigma_beta = Sigma_beta,
#   s = s,
#   tau = tau,
#   plot = TRUE,
#   mute = TRUE,
#   parallel = TRUE
# )

## ----BSET-X-load, include=FALSE-----------------------------------------------

## ----BSET-X-plot, echo=FALSE, out.width="100%", fig.align='center', fig.cap="Posterior distribution of $\\theta$ from the BSET procedure adjusting for covariates. Dashed lines show estimated quantities: the Bayesian 95\\% credible interval upper bound and threshold $\\eta$ (<span style='color:#0072B2'>blue shades</span>), and the frequentist 95\\% confidence interval upper bound and threshold $\\varepsilon$ (<span style='color:#D55E00'>orange shades</span>). Solid lines show the true values of $\\theta$ (<span style='color:#551A8B'>purple</span>) and $\\delta$ (<span style='color:#CC79A7'>pink</span>)."----
knitr::include_graphics("BSET_X_plot.png")

