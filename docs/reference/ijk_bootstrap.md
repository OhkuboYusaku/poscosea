# Bootstrap Using Infinitesimal Jackknife Approximation

Generates bootstrap replicates using infinitesimal jackknife (IJK)
influence values.

## Usage

``` r
ijk_bootstrap(estimate, influence, B = 1000)
```

## Arguments

- estimate:

  Finite numeric scalar giving the original estimate.

- influence:

  Non-empty finite numeric vector of influence values.

- B:

  Integer of at least 2 giving the number of bootstrap replicates.

## Value

An object of class "ijk_boot" containing:

- estimate:

  Original estimate.

- replicates:

  Bootstrap replicates.

- B:

  Number of bootstrap replicates.

## Details

The bootstrap distribution is obtained by resampling observation weights
from a multinomial distribution and applying the linear
influence-function approximation.

## Examples

``` r
set.seed(123)

influence <- rnorm(20)

boot <- ijk_bootstrap(
  estimate = 0,
  influence = influence,
  B = 200
)

boot
#> 
#> IJK Bootstrap Summary
#> 
#> Original estimate: 0 
#> 
#> Bootstrap replicates: 200 
#> 
#> Bias: 0.095 
#> Bootstrap SE: 4.203 
#> 95% percentile interval:
#> [-7.861, 7.632]
summary(boot)
#> $estimate
#> [1] 0
#> 
#> $bias
#> [1] 0.09549269
#> 
#> $se
#> [1] 4.203311
#> 
#> $ci
#>      2.5%     97.5% 
#> -7.860694  7.631705 
#> 
```
