# poscosea

`poscosea` provides infinitesimal jackknife (IJK) approximations for
sensitivity analysis of posterior summaries. Given aligned posterior draws
of a quantity of interest and pointwise log-likelihood values, it computes
observation-level influence values for bootstrap and leave-k-out analyses.

For binomial N-mixture models, the package also provides site- and
visit-level likelihood components and a two-stage bootstrap approximation.
See the [N-mixture vignette](docs/articles/nmix-twostage.html) for a worked example.

## Installation

Install the development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("OhkuboYusaku/poscosea")
```

## Quick start

This example uses the model `y_i | theta ~ Normal(theta, 1)` with a flat prior
on `theta`. It draws from the resulting analytic posterior and computes the
pointwise log-likelihood at each corresponding draw.

```r
set.seed(123)
n_obs <- 6
y <- rnorm(n_obs, sd = 1)
signal_idx <- c(2, 5)
y[signal_idx] <- y[signal_idx] + 3

n_draws <- 2000
theta <- rnorm(n_draws, mean = mean(y), sd = 1 / sqrt(n_obs))
loglik <- matrix(0, n_draws, n_obs)
for (i in seq_len(n_obs)) {
  loglik[, i] <- dnorm(y[i], mean = theta, sd = 1, log = TRUE)
}

fit <- ijk(theta, loglik)
fit$estimate
head(fit$influence)
```

Generate bootstrap replicates and approximate leave-one-out estimates:

```r
boot <- ijk_bootstrap(fit$estimate, fit$influence, B = 500)
summary(boot)

loo <- ijk_lko(fit$estimate, fit$influence, k = 1)
plot(loo)
```

`k` must be an integer from 1 to the number of observations. Exhaustive
leave-k-out enumeration is limited to at most one million omission sets.

## N-mixture models

The package includes `nmix_draws()`, `nmix_loglik()`, `ijk_nmix()`, and
`ijk_bootstrap_2stage()` for binomial N-mixture models. For full setup and
examples, see `vignette("nmix-twostage", package = "poscosea")`.

## Documentation

- [Function reference](docs/reference/index.html)
- [Getting started](docs/articles/getting-started.html)
- [Two-stage N-mixture bootstrap](docs/articles/nmix-twostage.html)
