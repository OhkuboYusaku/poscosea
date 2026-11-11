# Log-Likelihood Components of an N-Mixture Model

Computes, for every posterior draw, the site-level marginal
log-likelihood and the observation (visit)-level log-likelihood of a
binomial N-mixture model \$\$N_i \sim \mathrm{Poisson}(\lambda_i), \quad
y\_{ij} \mid N_i \sim \mathrm{Binomial}(N_i, p\_{ij}),\$\$ for sites
\\(i = 1, \ldots, I\\) and visits \\(j = 1, \ldots, J\\). These are the
inputs of
[`ijk_nmix()`](https://ohkuboyusaku.github.io/poscosea/reference/ijk_nmix.md)
and the two-stage bootstrap
[`ijk_bootstrap_2stage()`](https://ohkuboyusaku.github.io/poscosea/reference/ijk_bootstrap_2stage.md).

## Usage

``` r
nmix_loglik(y, lambda, p, N = NULL, Nmax_extra = 100)
```

## Arguments

- y:

  Numeric counts with \\(I\\) rows (sites) and \\(J\\) columns (visits);
  a numeric vector is treated as a one-column matrix. Counts must be
  finite non-negative integers, except that `NA` marks a missing visit.
  Every site must have at least one observed visit.

- lambda:

  Finite positive numeric matrix (\\(S \times I\\)) of posterior draws
  of \\(\lambda_i\\); at least one draw is required.

- p:

  Posterior draws of \\(p\_{ij}\\): either an array (\\(S \times I
  \times J\\)) or a matrix (\\(S \times IJ\\)) whose columns are ordered
  as \\(k = i + (j - 1) I\\) (the column-major order used by JAGS for
  `p[i,j]`). Values for observed visits must be finite and between 0 and
  1; values for missing visits are not used.

- N:

  Optional numeric matrix (\\(S \times I\\)) of finite non-negative
  integer posterior draws of \\(N_i\\); each draw must be at least the
  largest observed count at its site. Required for `obs_joint`.

- Nmax_extra:

  Non-negative integer; the summation grid for \\(N_i\\) extends
  `Nmax_extra` above \\(\max_j y\_{ij}\\).

## Value

An object of class `"nmix_loglik"` containing

- site:

  \\(S \times I\\) matrix of site-level marginal log-likelihoods.

- obs_rb:

  \\(S \times IJ\\) matrix (see Details). Missing visits are 0.

- obs_joint:

  \\(S \times IJ\\) matrix or `NULL`.

- condN:

  \\(S \times I\\) matrix of \\(E\[N_i \mid y_i, \theta\]\\).

- obs_mask:

  Logical \\(I \times J\\) matrix, `TRUE` for observed visits.

- tail_mass:

  Largest conditional probability at the upper end of the grid; a
  warning is issued when it exceeds `1e-8`.

## Details

The latent abundance \\(N_i\\) is summed out over the grid \\(N = \max_j
y\_{ij}, \ldots, \max_j y\_{ij} + \\) `Nmax_extra`. Missing visits (`NA`
in `y`) are dropped from the likelihood.

The observation-level log-likelihood is returned in two versions:

- `obs_rb`:

  \\(E\[\log \mathrm{Bin}(y\_{ij} \mid N, p\_{ij}) \mid y_i,
  \theta\]\\), the expectation over the conditional posterior of
  \\(N_i\\). Use it for quantities that depend on the parameters only
  (e.g. regression coefficients, \\(\lambda\\)). It has smaller Monte
  Carlo noise.

- `obs_joint`:

  \\(\log \mathrm{Bin}(y\_{ij} \mid N_i^{(s)}, p\_{ij}^{(s)})\\)
  evaluated at the posterior draws of \\(N_i\\). Use it for quantities
  that depend on the latent abundance (e.g. \\(\sum_i N_i\\)). Only
  available when `N` is supplied.

## See also

[`nmix_draws()`](https://ohkuboyusaku.github.io/poscosea/reference/nmix_draws.md),
[`ijk_nmix()`](https://ohkuboyusaku.github.io/poscosea/reference/ijk_nmix.md),
[`ijk_bootstrap_2stage()`](https://ohkuboyusaku.github.io/poscosea/reference/ijk_bootstrap_2stage.md)

## Examples

``` r
set.seed(1)
I <- 10; J <- 3; S <- 50
y <- matrix(rpois(I * J, 3), I, J)
lambda <- matrix(rgamma(S * I, 20, 4), S, I)
p <- matrix(runif(S * I * J, 0.4, 0.7), S, I * J)
ll <- nmix_loglik(y, lambda, p)
dim(ll$site)
#> [1] 50 10
```
