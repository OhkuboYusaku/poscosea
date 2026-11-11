# Two-Stage Bootstrap Weights for Site/Visit Data

Draws bootstrap weights for hierarchical data with sites \\(i = 1,
\ldots, I\\) and visits \\(j = 1, \ldots, J\\). Sites are resampled with
replacement, \\(u \sim \mathrm{Mult}(I; 1/I, \ldots, 1/I)\\); then,
within each selected copy of site \\(i\\), its \\(J_i\\) observed visits
are resampled with replacement, so that the total weight of visit
\\(j\\) is \\(W\_{ij} \sim \mathrm{Mult}(u_i J_i; 1/J_i, \ldots,
1/J_i)\\) and the visit weight is \\(v\_{ij} = W\_{ij} / u_i\\). Visits
of sites with \\(u_i = 0\\) or \\(J_i = 1\\), and missing visits, get
\\(v\_{ij} = 1\\).

## Usage

``` r
nmix_weights(B, obs_mask)
```

## Arguments

- B:

  Positive integer giving the number of bootstrap replicates.

- obs_mask:

  Non-empty logical \\(I \times J\\) matrix of observed visits without
  missing values, e.g. `loglik$obs_mask` from
  [`nmix_loglik()`](https://ohkuboyusaku.github.io/poscosea/reference/nmix_loglik.md)
  or `!is.na(y)`.

## Value

A list with `site` (\\(B \times I\\)) and `visit` (\\(B \times IJ\\),
columns ordered as \\(i + (j-1)I\\)).

## Details

Supplying the same weights to several calls of
[`ijk_bootstrap_2stage()`](https://ohkuboyusaku.github.io/poscosea/reference/ijk_bootstrap_2stage.md)
gives bootstrap replicates that are coherent across quantities.

## Examples

``` r
set.seed(1)
w <- nmix_weights(5, matrix(TRUE, 4, 3))
rowSums(w$site)
#> [1] 4 4 4 4 4
```
