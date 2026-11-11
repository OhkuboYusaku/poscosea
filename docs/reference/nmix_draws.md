# Extract N-Mixture Quantities from MCMC Draws

Picks the posterior draws of \\(\lambda_i\\), \\(p\_{ij}\\) and
(optionally) \\(N_i\\) out of a matrix of MCMC draws whose columns are
named as in JAGS/BUGS output (`lambda[1]`, `p[1,1]`, `N[1]`, ...).

## Usage

``` r
nmix_draws(draws, I, J, lambda = "lambda", p = "p", N = "N")
```

## Arguments

- draws:

  A non-empty numeric matrix or data frame of finite draws with unique,
  named columns, or a non-empty list of matrices/data frames with the
  same columns (chains are stacked).

- I:

  Positive integer giving the number of sites.

- J:

  Positive integer giving the number of visits.

- lambda, p, N:

  Non-empty node-name prefixes for \\(\lambda\\), \\(p\\) and \\(N\\).
  Set `N = NULL` to skip \\(N\\). If not all requested `N` columns are
  present, the returned `N` element is `NULL`.

## Value

A list with `lambda` (\\(S \times I\\)), `p` (\\(S \times IJ\\), columns
ordered as \\(i + (j-1)I\\)) and `N` (\\(S \times I\\) or `NULL` when
the columns are absent).

## Examples

``` r
I <- 2; J <- 2; S <- 3
draws <- matrix(runif(S * (I + I * J)), S, I + I * J)
colnames(draws) <- c(paste0("lambda[", 1:I, "]"),
                     paste0("p[", rep(1:I, J), ",", rep(1:J, each = I), "]"))
str(nmix_draws(draws, I, J))
#> List of 3
#>  $ lambda: num [1:3, 1:2] 0.9972 0.3039 0.1489 0.3285 0.0705 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr [1:2] "lambda[1]" "lambda[2]"
#>  $ p     : num [1:3, 1:4] 0.14 0.49 0.297 0.303 0.441 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr [1:4] "p[1,1]" "p[2,1]" "p[1,2]" "p[2,2]"
#>  $ N     : NULL
```
