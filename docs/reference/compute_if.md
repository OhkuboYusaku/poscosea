# Compute Infinitesimal Jackknife Influence Values

Computes observation-wise influence values used in the infinitesimal
jackknife approximation.

## Usage

``` r
compute_if(theta, loglik)
```

## Arguments

- theta:

  Finite numeric vector of at least two posterior draws.

- loglik:

  Finite numeric matrix with one row per draw in `theta` and at least
  one observation column.

## Value

A numeric vector of influence values with one element per observation.
