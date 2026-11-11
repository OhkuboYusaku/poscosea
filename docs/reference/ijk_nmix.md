# Infinitesimal Jackknife Influence for an N-Mixture Model

Computes site-level and observation (visit)-level influence values, i.e.
posterior covariances between the quantity of interest and the
log-likelihood components returned by
[`nmix_loglik()`](https://ohkuboyusaku.github.io/poscosea/reference/nmix_loglik.md).
These are used by
[`ijk_bootstrap_2stage()`](https://ohkuboyusaku.github.io/poscosea/reference/ijk_bootstrap_2stage.md)
for the two-stage (site and visit) bootstrap approximation.

## Usage

``` r
ijk_nmix(theta, loglik, obs = c("rb", "joint"))
```

## Arguments

- theta:

  Finite numeric vector (length \\(S\\), with \\(S \ge 2\\)) of
  posterior draws of the quantity of interest, or a finite numeric \\(S
  \times K\\) matrix for \\(K\\) quantities at once (e.g. a regression
  curve on a grid). The number of rows must match `loglik`.

- loglik:

  An object of class `"nmix_loglik"`.

- obs:

  Which observation-level log-likelihood to use: `"rb"` for quantities
  depending only on the parameters, `"joint"` for quantities that depend
  on the latent abundance \\(N\\).

## Value

An object of class `"ijk_nmix"` containing

- estimate:

  Posterior mean(s).

- influence_site:

  Site-level influence values (centred, as in
  [`ijk()`](https://ohkuboyusaku.github.io/poscosea/reference/ijk.md));
  a vector of length \\(I\\), or an \\(I \times K\\) matrix.

- influence_obs:

  Observation-level influence values; length \\(IJ\\) (ordered \\(i +
  (j-1)I\\)), 0 for missing visits.

- obs_mask:

  Observed-visit indicator.

- obs:

  The `obs` argument.

The site-level influence can be passed to
[`ijk_lko()`](https://ohkuboyusaku.github.io/poscosea/reference/ijk_lko.md)
for an approximate leave-one-site-out analysis.

## Details

The first-order approximation of the posterior mean under site weights
\\(u_i\\) and visit weights \\(v\_{ij}\\) is \$\$E\_{u,v}\[\theta\]
\approx E\[\theta\] + \sum_i (u_i - 1) c_i + \sum\_{i,j} (v\_{ij} - 1)
e\_{ij},\$\$ with \\(c_i = \mathrm{Cov}\[\theta, \ell_i\]\\) (site-level
marginal log-likelihood) and \\(e\_{ij} = \mathrm{Cov}\[\theta,
\ell\_{ij}\]\\) (observation-level log-likelihood).

## See also

[`nmix_loglik()`](https://ohkuboyusaku.github.io/poscosea/reference/nmix_loglik.md),
[`ijk_bootstrap_2stage()`](https://ohkuboyusaku.github.io/poscosea/reference/ijk_bootstrap_2stage.md)

## Examples

``` r
set.seed(1)
I <- 10; J <- 3; S <- 200
y <- matrix(rpois(I * J, 3), I, J)
lambda <- matrix(rgamma(S * I, 20, 4), S, I)
p <- matrix(runif(S * I * J, 0.4, 0.7), S, I * J)
ll <- nmix_loglik(y, lambda, p)
fit <- ijk_nmix(rowSums(lambda), ll)
fit$estimate
#> [1] 49.89695
```
