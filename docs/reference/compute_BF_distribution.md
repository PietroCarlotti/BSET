# Compute the Distribution of the Bayes Factor from Carlotti and Parast (2026)

This function calculates the probability mass function and cumulative
distribution function of the Bayes factor defined in Carlotti and Parast
(2026) for the following hypothesis test: \$\$\begin{cases} H_0: V_S =
V_S^{0} \\ H_1: V_S \> V_S^{0} \end{cases}\$\$ where \\V_S\\ is the
surrogate's treatment effect on \\S\\ measured as the probability
\$\$V_S = P(S\_{1i} \> S\_{0i})\$\$ and \\V_S^{0}\\ is a hypothesized
value under the null hypothesis. These hypotheses can be tested by
fitting the following Beta-binomial model to the data: \$\$\hat{V}\_S
\mid V_S \sim \text{Binomial} (n, V_S)\$\$ \$\$V_S \sim \text{Beta} (a,
b),\$\$ where \\\hat{V}\_S\\ is the sample estimate of the surrogate's
treatment effect on \\S\\ computed as \$\$\hat{V}\_S = \frac{1}{n}
\sum\limits^{n}\_{i=1} I(S\_{1i} \> S\_{0i}).\$\$ The Bayes factor is
then computed as the ratio of the marginal likelihoods under the
alternatives: \$\$BF\_{n} = \frac{1 - F\_{\text{Beta} (a + n \hat{V}\_S,
b + n - n \hat{V}\_S)} (\frac{1}{2})}{1 - F\_{\text{Beta} (a, b)}
(\frac{1}{2})},\$\$ where \\F\_{\text{Beta}}\\ is the cumulative
distribution function of a Beta distribution, \\a\\ and \\b\\ are the
shape parameters of the Beta prior. Given the true value of \\V_S\\, the
distribution of the Bayes factor can be computed by evaluating \\BF_n\\
for all possible values of \\\hat{V}\_S\\ and their corresponding
probabilities under the Binomial distribution with parameters \\n\\ and
the true value of \\V_S\\.

## Usage

``` r
compute_BF_distribution(
  n,
  V_S_true,
  V_S_zero = 0.5,
  a = 1,
  b = 1,
  BF_alternative
)
```

## Arguments

- n:

  Integer. The sample size.

- V_S_true:

  Numeric. The true value of treatment effect on the surrogate.

- V_S_zero:

  Numeric. The hypothesized value of the surrogate's treatment effect
  under the null hypothesis (default is 0.5).

- a:

  Numeric. First shape parameter alpha for the Beta prior (default is
  1).

- b:

  Numeric. Second shape parameter beta for the Beta prior (default is
  1).

- BF_alternative:

  Character. The type of alternative hypothesis: either `"two_sided"` or
  `"greater"`.

## Value

A data frame containing:

- `BF_values`: The possible values of the Bayes Factor.

- `BF_PMF`: The probability mass function for the Bayes Factor.

- `BF_CDF`: The cumulative distribution function for the Bayes Factor.
