# Infinitesimal Jackknife Approximation for MCMC Output

Computes an infinitesimal jackknife (IJK) approximation of the posterior
mean from MCMC draws and corresponding log-likelihood values.

## Usage

``` r
ijk(theta, loglik)
```

## Arguments

- theta:

  Finite numeric vector of posterior draws of a scalar quantity, with at
  least two draws.

- loglik:

  Finite numeric matrix of log-likelihood draws. Rows correspond to
  posterior draws and columns correspond to observations; it must have
  at least one column and one row per value in `theta`.

## Value

An object of class `"ijk"` containing:

- estimate:

  Posterior mean estimate.

- influence:

  Observation-wise influence values.

- theta:

  Input posterior draws.

- loglik:

  Input log-likelihood matrix.

## Details

The function estimates observation-wise influence values that can be
used for bootstrap and leave-k-out cross-validation approximations
without refitting the original model.

## Examples

``` r
set.seed(123)
theta <- rnorm(1000)
loglik <- matrix(rnorm(1000 * 20), 1000, 20)

fit <- ijk(theta, loglik)
fit$estimate
#> [1] 0.01612787
```
