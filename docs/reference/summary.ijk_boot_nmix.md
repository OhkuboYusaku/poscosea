# Summarize a Two-Stage IJK Bootstrap Object

Summarize a Two-Stage IJK Bootstrap Object

## Usage

``` r
# S3 method for class 'ijk_boot_nmix'
summary(object, level = 0.95, ...)

# S3 method for class 'ijk_boot_nmix'
print(x, digits = 3, ...)
```

## Arguments

- object, x:

  An object of class `"ijk_boot_nmix"`.

- level:

  Finite confidence level strictly between 0 and 1 for the percentile
  interval.

- ...:

  Currently unused.

- digits:

  Number of digits to display.

## Value

[`summary()`](https://rdrr.io/r/base/summary.html) returns a data frame
with one row per quantity: `estimate`, `bias`, `se` and the percentile
interval `lower`, `upper`.
[`print()`](https://rdrr.io/r/base/print.html) returns `x` invisibly.

## Examples

``` r
set.seed(1)
I <- 10; J <- 3; S <- 200
y <- matrix(rpois(I * J, 3), I, J)
lambda <- matrix(rgamma(S * I, 20, 4), S, I)
p <- matrix(runif(S * I * J, 0.4, 0.7), S, I * J)
fit <- ijk_nmix(rowSums(lambda), nmix_loglik(y, lambda, p))
boot <- ijk_bootstrap_2stage(fit, B = 200)
print(boot)
#> 
#> IJK Bootstrap Summary (N-mixture, type = "twostage")
#> 
#> Bootstrap replicates: 200 
#> 
#>   estimate   bias   se lower upper
#> 1     49.9 -0.121 1.19  47.3  52.2
```
