# Visit Weights Given Site Weights

Internal helper of
[`nmix_weights()`](https://ohkuboyusaku.github.io/poscosea/reference/nmix_weights.md).

## Usage

``` r
visit_weights(U, obs_mask)
```

## Arguments

- U:

  Non-empty matrix (\\(B \times I\\)) of non-negative integer site
  weights.

- obs_mask:

  Non-empty logical \\(I \times J\\) matrix of observed visits without
  missing values; its row count must match `ncol(U)`.

## Value

Numeric matrix (\\(B \times IJ\\)) of visit weights.
