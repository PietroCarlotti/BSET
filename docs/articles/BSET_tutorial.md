# BSET Tutorial

## Introduction

`BSET` (Bayesian Surrogate Evaluation Test) is an `R` package for
assessing the validity of surrogate markers in clinical trials. It
provides hypothesis testing tools to evaluate whether a surrogate can
reliably estimate the causal effect of a treatment on a primary outcome.
The package implements the imputation-based Bayesian methodology of
Carlotti and Parast (2026), extending the frequentist rank-based
approach of Parast et al. (2024). BSET addresses key limitations of the
frequentist method, including the lack of causal interpretability and
the inability to adjust for covariates in the estimation process.

The package supports Bayesian testing both with and without baseline
covariates. Additionally, it includes comprehensive simulation suites to
replicate studies from both papers, enabling performance comparisons
between Bayesian and frequentist approaches across diverse clinical
scenarios.

## Installation

The package is not yet available on CRAN, but you can install the
development version of BSET from [GitHub](https://github.com/) with:

``` r
# install.packages("pak")
# pak::pak("PietroCarlotti/BSET")
```

To follow along with the examples in this tutorial, you will also need
to install the following packages:

``` r
library(ggplot2)
library(dplyr)
```

## Generating the data

We start by setting the random seed for reproducibility of the results.

``` r
# Set the random seed for reproducibility
set.seed(123)
```

The package includes two functions to generate data for the simulations:
`DGP` and `DGP_X`. The first one generates data from the simulation
settings of Parast et al. (2024), which do not include covariates; the
second one generates data which depends on a binary covariate $`X`$ and
is used for the simulations of Carlotti and Parast (2026).

For example, this how we can generate data from the second setting
defined in Parast et al. (2024) (i.e.: the setting where the surrogate
is perfect and covariates are not included):

``` r
# Sample size
n <- 6

# Treatment assignment probability
p <- 0.5

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
```

While this is how we can generate data from the first setting defined in
Carlotti and Parast (2026) (i.e.: the setting where the surrogate is
useful and covariates are included):

``` r
# Sample size
n <- 6

# Treatment assignment probability
p <- 0.5

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
```

## Computing the true values of $`\delta`$ and $`\theta`$

Given the data generated with the above functions, we can estimate
$`\delta`$ and $`\theta`$ using the `compute_delta` and `compute_theta`
functions, respectively.

For example, for the data generated without covariates:

``` r
# Compute delta
delta_no_X <- BSET::compute_delta(MC_data = data_no_X)

# Compute theta
theta_no_X <- BSET::compute_theta(MC_data = data_no_X)
```

And for the data generated with covariates:

``` r
# Compute delta
delta_X <- BSET::compute_delta(MC_data = data_X)

# Compute theta
theta_X <- BSET::compute_theta(MC_data = data_X)
```

The estimated values of $`\delta`$ and $`\theta`$ are shown in Table
@ref(tab:delta-theta-table).

|         Setting          | $`\widehat{\delta}`$ | $`\widehat{\theta}`$ |
|:------------------------:|:--------------------:|:--------------------:|
| No X (Perfect Surrogate) |        0.000         |        0.167         |
|         Binary X         |        0.875         |        0.000         |

Estimated values of $`\delta`$ and $`\theta`$ in the two settings.

To determine the values of $`\delta`$ and $`\theta`$ for all the
settings defined in the two papers, we use a large-scale Monte Carlo
simulation via the `compute_estimands_Parast_et_al_2024` function for
data without covariates, and the
`compute_estimands_Carlotti_and_Parast_2026` function for data including
covariates.

``` r
# Number of Monte Carlo samples to generate per setting
MC_samples <- 10#00000

# Compute the true values of delta and theta for all settings of Parast et al. (2024)
estimands_Parast_et_al_2024 <- BSET::compute_estimands_Parast_et_al_2024(MC_samples)

# Compute the true values of delta and theta for all settings of Carlotti and Parast (2026)
estimands_Carlotti_and_Parast_2026 <- BSET::compute_estimands_Carlotti_and_Parast_2026(MC_samples)
```

As shown in Table @ref(tab:estimands-Parast-et-al-2024), $`\delta`$ and
$`\theta`$ yield identical values (up to Monte Carlo error) across all
settings. Both estimands are close to $`0`$ for a perfect surrogate,
they are significantly higher for a useless surrogate, and take
intermediate values for imperfect surrogates or misspecified models.

|       Setting       | $`\delta`$ | $`\theta`$ |
|:-------------------:|:----------:|:----------:|
|  Useless surrogate  |   0.458    |    0.5     |
|  Perfect surrogate  |   0.042    |    0.1     |
| Imperfect surrogate |   0.083    |    0.1     |
| Misspecified model  |   0.238    |    0.2     |

Monte Carlo estimates of $`\delta`$ and $`\theta`$ for the simulation
settings considered in Parast et al. (2024).

Whereas, as shown in Table @ref(tab:estimands-Carlotti-and-Parast-2026),
when we include covariates in the data generating process, $`\delta`$
and $`\theta`$ yield different values. In particular, while $`\theta`$
is close to $`0`$ for a perfect surrogate, $`\delta`$ is significantly
higher.

|      Setting      | $`\delta`$ | $`\theta`$ |
|:-----------------:|:----------:|:----------:|
| Perfect surrogate |   0.286    |     0      |

Monte Carlo estimates of $`\delta`$ and $`\theta`$ for the simulation
settings considered in Carlotti and Parast (2026).

## Compute the validation threshold $`\eta`$

The validation threshold $`\eta`$ is the value of $`\theta`$ below which
the surrogate is considered valid. As explained in Carlotti and Parast
(2026), the computation of $`\eta`$ is based on the distribution of the
following Bayes factor:
``` math

  BF_{n} = \frac{1 - F_{\text{Beta} \left( a + n \hat{V}_{S}, \, b + n - n \hat{V}_{S} \right)} \left( \frac{1}{2} \right)}{1 - F_{\text{Beta} \left( a, \, b \right)} \left( \frac{1}{2} \right)}.
```
Given the true value of $`V_S`$, the distribution of $`BF_n`$ can be
computed with the function `compute_BF_distribution`. For example, the
distribution of the Bayes factor under the null hypothesis $`V_S = 0.5`$
can be computed as follows:

``` r
# Sample size
n <- 6

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
```

The distribution of the Bayes factor under the null hypothesis is shown
in Figure @ref(fig:BF-distribution-plot). Since the Bayes factor can
take values on a very wide range, we plot the distribution on the log
scale.

![Distribution of the Bayes factor under the null hypothesis \$V_S =
0.5\$.](BSET_tutorial_files/figure-html/BF-distribution-plot-1.png)

Distribution of the Bayes factor under the null hypothesis
$`V_S = 0.5`$.

Given the distribution of the Bayes factor under the null hypothesis, we
can compute the critical value $`BF_{n, \alpha}`$ corresponding to a
specified Type I error rate $`\alpha`$.

``` r
# Type I error rate
alpha <- 0.05

# Compute the critical value of the Bayes factor corresponding to alpha
BF_alpha <- BF_distribution %>%
    dplyr::filter(BF_CDF >= 1 - alpha) %>%
    dplyr::slice(1) %>%
    dplyr::pull(BF_values)
```

For example, for $`\alpha = 0.05`$ we have

``` math

  BF_{n, \alpha} = 2.857.
```

Once we have the value of $`BF_{n, \alpha}`$, we can compute the value
of $`v_S`$ that satisfies the following equation:

``` math

  P(\text{BF}_n \geq \text{BF}_{n, \alpha} \; | \; V_S = v_S) = 1 - \beta,
```

where $`1 - \beta`$ is the desired power of the test. A root-finding
algorithm is used to solve for the value of $`v_S`$ that satisfies the
above equation, which is implemented in the function `compute_V_S_star`.

For example, for a Type II error rate of $`\beta = 0.2`$, we can compute
the value of $`v_S`$ as follows:

``` r
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
```

Then, we have that

``` math

  v_S = 0.963.
```

Finally, the validation threshold $`\eta`$ can be computed as

``` math

  \eta = \max \{v_Y - v_S, 0\},
```

where $`v_Y`$ is the hypothesized value of the treatment effect on the
primary outcome (typically set equal to the estimate computed on the
available data).

For example, we can set $`v_Y`$ equal to the Monte Carlo estimate of
$`V_Y`$ for setting 1 of Parast et al. (2024) (i.e.: the setting where
the surrogate is perfect and covariates are not included).

``` r
# Hypothesized value of the treatment effect on the primary outcome
v_Y <- estimands_Parast_et_al_2024$V_Y_MC[2]

# Compute the validation threshold eta
eta <- max(v_Y - V_S_star, 0)
```

In this case, we have that

``` math

  \eta = 0.037.
```

## Run the BSET procedure

The package includes two functions to run the BSET procedure:
`BSET_no_X` and `BSET_X`. The first one runs the BSET procedure without
adjusting for covariates, while the second one runs the BSET procedure
adjusting for covariates.

### BSET with no covariates

As an example, we can run the BSET procedure without adjusting for
covariates on the data generated from the second setting of Parast et
al. (2024). First of all, we need to prepare the data in the format
required by the function `BSET_no_X`: a data frame with three columns,
where the first column contains the observed values of the primary
outcome $`Y`$, the second column contains the observed values of the
surrogate $`S`$, and the third column contains the treatment assignment
$`Z`$.

``` r
# Prepare the data for BSET_no_X
BSET_no_X_data <- data.frame(
  Y = data_no_X$P_observed[, "Y"],
  S = data_no_X$P_observed[, "S"],
  Z = data_no_X$Z 
)
```

We also need to define the posterior sampling parameters:

- `n_chains`: The number of MCMC chains to run.
- `n_iter`: The total number of MCMC iterations to run for each chain.
- `burn_in_ratio`: The proportion of MCMC iterations to discard as
  burn-in.

``` r
# Posterior sampling parameters
n_chains <- 2
n_iter <- 2000
burn_in_ratio <- 0.25
```

Then, we need to define the prior parameters for the distribution of the
potential outcomes in our Bayesian model. In particular, we need to
specify:

- `mu_0`: A numeric vector of length 4 containing the prior means for
  the potential outcomes.
- `Sigma_0`: A numeric positive-definite matrix of dimension
  $`4 \times 4`$ containing the prior covariance matrix for the
  potential outcomes.
- `s`: A numeric vector of length 4 containing the prior scale
  parameters for the potential outcomes.
- `tau`: A numeric scalar containing the prior scale parameter for the
  correlation matrix of the potential outcomes.

``` r
# Prior parameters
mu_0 <- rep(0, 4)
Sigma_0 <- diag(4)
s <- rep(1, 4)
tau <- 1
```

If the true values of $`\delta`$ and $`\theta`$ are known, we can also
specify them as input to the function `BSET_no_X` to check whether they
are included in the posterior credible intervals. In this case, we can
set $`\delta`$ and $`\theta`$ equal to the Monte Carlo estimates for
setting 2 of Parast et al. (2024).

We can also specify some additional parameters for the BSET procedure,
such as:

- `seed`: An integer scalar to set the random seed for reproducibility
  of the MCMC sampling process.
- `plot`: A logical scalar indicating whether to plot the posterior
  distribution of $`\theta`$.
- `mute`: A logical scalar indicating whether to suppress the output of
  the MCMC sampling process.
- `parallel`: A logical scalar indicating whether to run the MCMC
  sampling process in parallel or sequentially.

Finally, we can run the BSET procedure without adjusting for covariates
as follows:

``` r
# Run the BSET procedure without adjusting for covariates
BSET_no_X_results <- BSET::BSET_no_X(
  data = BSET_no_X_data,
  Y = "Y",
  S = "S",
  Z = "Z",
  delta_true = estimands_Parast_et_al_2024$delta_MC[2],
  theta_true = estimands_Parast_et_al_2024$theta_MC[2],
  seed = 123,
  n_chains = n_chains,
  n_iter = n_iter,
  burn_in_ratio = burn_in_ratio,
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
  parallel = TRUE
)
```

The posterior distribution of $`\theta`$ from the BSET procedure without
adjusting for covariates is shown in Figure @ref(fig:BSET-no-X-plot).

![Posterior distribution of \$\theta\$ from the BSET procedure without
adjusting for covariates. The \<span style='color:blue'\>blue\</span\>
vertical line indicates the upper bound of the 95% credible interval,
the \<span style='color:green'\>green\</span\> vertical line indicates
value of the validation threshold \$\eta\$, the \<span
style='color:orange'\>orange\</span\> vertical line indicates the true
value of \$\delta\$, and the \<span style='color:red'\>red\</span\>
vertical line indicates the true value of
\$\theta\$.](BSET_tutorial_files/figure-html/BSET-no-X-plot-1.png)

Posterior distribution of $`\theta`$ from the BSET procedure without
adjusting for covariates. The blue vertical line indicates the upper
bound of the 95% credible interval, the green vertical line indicates
value of the validation threshold $`\eta`$, the orange vertical line
indicates the true value of $`\delta`$, and the red vertical line
indicates the true value of $`\theta`$.

### BSET with covariates

As an example, we can run the BSET procedure adjusting for covariates on
the data generated from the first setting of Carlotti and Parast (2026).
First of all, we need to prepare the data in the format required by the
function `BSET_X`: a data frame with four columns, where the first
column contains the observed values of the primary outcome $`Y`$, the
second column contains the observed values of the surrogate $`S`$, the
third column contains the treatment assignment $`Z`$, and the fourth
column contains the values of the binary covariate $`X`$. We call $`d`$
the number of covariates, which in this case is equal to $`2`$
(including the intercept).

``` r
# Prepare the data for BSET_X
BSET_X_data <- data.frame(
  Y = data_X$P_observed[, "Y"],
  S = data_X$P_observed[, "S"],
  Z = data_X$Z,
  X = data_X$X
)
```

The posterior sampling parameters are defined in the same way as for the
BSET procedure without adjusting for covariates. Whereas the prior
parameters for the distribution of the potential outcomes in our
Bayesian model are defined as follows:

- `mu_beta`: A numeric vector of length $`d`$ containing the prior means
  for the regression coefficients of the potential outcomes on the
  covariates.
- `Sigma_beta`: A numeric positive-definite matrix of dimension
  $`d \times d`$ containing the prior covariance matrix for the
  regression coefficients of the potential outcomes on the covariates.

``` r
# Prior parameters for the regression coefficients of the potential outcomes on the covariates
d <- 2
mu_beta <- rep(0, d)
Sigma_beta <- 10*diag(d)
```

We can also specify the true values of $`\delta`$ and $`\theta`$ as
input to the function `BSET_X` to check whether if they are included in
the posterior credible intervals. In this case, we can set $`\delta`$
and $`\theta`$ equal to the Monte Carlo estimates of $`\delta`$ and
$`\theta`$ for setting 1 of Carlotti and Parast (2026).

Finally, we can run the BSET procedure adjusting for covariates as
follows:

``` r
# Run the BSET procedure adjusting for covariates
BSET_X_results <- BSET::BSET_X(
  data = BSET_X_data,
  Y = "Y",
  S = "S",
  Z = "Z",
  X = "X",
  delta_true = estimands_Carlotti_and_Parast_2026$delta_MC[1],
  theta_true = estimands_Carlotti_and_Parast_2026$theta_MC[1],
  seed = 123,
  n_chains = n_chains,
  n_iter = n_iter,
  burn_in_ratio = burn_in_ratio,
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
  parallel = TRUE
)
```

The posterior distribution of $`\theta`$ from the BSET procedure
adjusting for covariates is shown in Figure @ref(fig:BSET-X-plot).

![Posterior distribution of \$\theta\$ from the BSET procedure adjusting
for covariates. The \<span style='color:blue'\>blue\</span\> vertical
line indicates the upper bound of the 95% credible interval, the \<span
style='color:green'\>green\</span\> vertical line indicates value of the
validation threshold \$\eta\$, the \<span
style='color:orange'\>orange\</span\> vertical line indicates the true
value of \$\delta\$, and the \<span style='color:red'\>red\</span\>
vertical line indicates the true value of
\$\theta\$.](BSET_tutorial_files/figure-html/BSET-X-plot-1.png)

Posterior distribution of $`\theta`$ from the BSET procedure adjusting
for covariates. The blue vertical line indicates the upper bound of the
95% credible interval, the green vertical line indicates value of the
validation threshold $`\eta`$, the orange vertical line indicates the
true value of $`\delta`$, and the red vertical line indicates the true
value of $`\theta`$.
