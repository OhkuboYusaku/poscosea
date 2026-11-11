# Summarize an IJK Bootstrap Object

Computes standard bootstrap summary statistics.

## Usage

``` r
# S3 method for class 'ijk_boot'
summary(object, ...)
```

## Arguments

- object:

  An object of class "ijk_boot".

- ...:

  Currently unused.

## Value

A list containing:

- estimate:

  Original estimate.

- bias:

  Bootstrap bias estimate.

- se:

  Bootstrap standard error.

- ci:

  95% percentile confidence interval.

## Examples

``` r
set.seed(123)
boot <- ijk_bootstrap(
  estimate = 0,
  influence = rnorm(20),
  B = 100
)

summary(boot)
#> $estimate
#> [1] 0
#> 
#> $bias
#> [1] -0.7199799
#> 
#> $se
#> [1] 4.195314
#> 
#> $ci
#>      2.5%     97.5% 
#> -8.140657  7.197027 
#> 
```
