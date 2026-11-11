# Print an IJK Bootstrap Object

Displays a concise summary of an object of class "ijk_boot".

## Usage

``` r
# S3 method for class 'ijk_boot'
print(x, digits = 3, ...)
```

## Arguments

- x:

  An object of class "ijk_boot".

- digits:

  Number of digits to display.

- ...:

  Currently unused.

## Value

The input object, invisibly.

## Examples

``` r
set.seed(123)
boot <- ijk_bootstrap(
  estimate = 0,
  influence = rnorm(20),
  B = 100
)

print(boot)
#> 
#> IJK Bootstrap Summary
#> 
#> Original estimate: 0 
#> 
#> Bootstrap replicates: 100 
#> 
#> Bias: -0.72 
#> Bootstrap SE: 4.195 
#> 95% percentile interval:
#> [-8.141, 7.197]
```
