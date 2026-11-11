# Leave-k-Out Approximation

Approximates leave-k-out cross-validation estimates using infinitesimal
jackknife influence values.

## Usage

``` r
ijk_lko(estimate, influence, k = 1)
```

## Arguments

- estimate:

  Finite numeric scalar giving the original estimate.

- influence:

  Non-empty finite numeric vector of influence values.

- k:

  Finite integer number of observations removed, satisfying
  `1 <= k <= length(influence)`. Exhaustive enumeration is limited to at
  most 1e6 omission sets.

## Value

An object of class `ijk_lko` containing:

- estimate:

  Original estimate.

- lko_estimates:

  Approximate estimates for all omitted sets.

- delta:

  Change from the original estimate for each omitted set.

- omitted:

  List of omitted observation indices for each estimate.

- k:

  Number of observations omitted in each set.

## Examples

``` r
influence <- rnorm(20)
ijk_lko(0, influence, k = 1)
#> $estimate
#> [1] 0
#> 
#> $lko_estimates
#>  [1]  0.828236578  0.074745458  0.454608335  0.898838790  1.386837727
#>  [6]  0.812561346 -0.510760550 -0.301850973  1.313677295  0.607611113
#> [11]  0.292080516 -0.648173999 -0.821998193  0.786546092 -0.157369194
#> [16]  0.009683662 -2.117555455  0.510126219 -1.721826447  0.138776467
#> 
#> $delta
#>  [1]  0.828236578  0.074745458  0.454608335  0.898838790  1.386837727
#>  [6]  0.812561346 -0.510760550 -0.301850973  1.313677295  0.607611113
#> [11]  0.292080516 -0.648173999 -0.821998193  0.786546092 -0.157369194
#> [16]  0.009683662 -2.117555455  0.510126219 -1.721826447  0.138776467
#> 
#> $omitted
#> $omitted[[1]]
#> [1] 1
#> 
#> $omitted[[2]]
#> [1] 2
#> 
#> $omitted[[3]]
#> [1] 3
#> 
#> $omitted[[4]]
#> [1] 4
#> 
#> $omitted[[5]]
#> [1] 5
#> 
#> $omitted[[6]]
#> [1] 6
#> 
#> $omitted[[7]]
#> [1] 7
#> 
#> $omitted[[8]]
#> [1] 8
#> 
#> $omitted[[9]]
#> [1] 9
#> 
#> $omitted[[10]]
#> [1] 10
#> 
#> $omitted[[11]]
#> [1] 11
#> 
#> $omitted[[12]]
#> [1] 12
#> 
#> $omitted[[13]]
#> [1] 13
#> 
#> $omitted[[14]]
#> [1] 14
#> 
#> $omitted[[15]]
#> [1] 15
#> 
#> $omitted[[16]]
#> [1] 16
#> 
#> $omitted[[17]]
#> [1] 17
#> 
#> $omitted[[18]]
#> [1] 18
#> 
#> $omitted[[19]]
#> [1] 19
#> 
#> $omitted[[20]]
#> [1] 20
#> 
#> 
#> $k
#> [1] 1
#> 
#> attr(,"class")
#> [1] "ijk_lko"
```
