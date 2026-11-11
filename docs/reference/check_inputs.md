# Validate Inputs

Internal utility for validating inputs supplied to the main ijk
function.

## Usage

``` r
check_inputs(theta, loglik)
```

## Arguments

- theta:

  Finite numeric vector with at least two values.

- loglik:

  Finite numeric matrix with at least one column and one row per value
  in `theta`.

## Value

Invisibly returns TRUE.
