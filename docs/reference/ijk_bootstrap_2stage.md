# Two-Stage Bootstrap Approximation for an N-Mixture Model

Generates bootstrap replicates of a posterior mean from a single MCMC
run, using the first-order (infinitesimal jackknife) approximation with
site-level and visit-level weights: \$\$\theta^{(b)} \approx
E\[\theta\] + \sum_i (u\_{bi} - 1) c_i + \sum\_{i,j} (v\_{bij} - 1)
e\_{ij}.\$\$

## Usage

``` r
ijk_bootstrap_2stage(
  object,
  B = 1000,
  type = c("twostage", "site", "visit"),
  weights = NULL,
  latent_resid = NULL
)
```

## Arguments

- object:

  An object of class `"ijk_nmix"`.

- B:

  Positive integer giving the number of bootstrap replicates when
  `weights` is `NULL`; ignored when `weights` is supplied.

- type:

  One of `"twostage"`, `"site"`, `"visit"`.

- weights:

  Optional weights from
  [`nmix_weights()`](https://ohkuboyusaku.github.io/poscosea/reference/nmix_weights.md),
  to share the same resamples across quantities.

- latent_resid:

  Optional finite numeric vector (length \\(S\\)) or finite numeric
  matrix with \\(S\\) rows and \\(K\\) columns, containing latent
  prediction errors; one randomly chosen draw is added to each
  replicate.

## Value

An object of class `c("ijk_boot_nmix", "ijk_boot")` with `estimate`,
`replicates` (vector of length \\(B\\), or \\(B \times K\\) matrix), `B`
and `type`.

## Details

- `"twostage"`:

  Sites and, within them, visits are resampled (see
  [`nmix_weights()`](https://ohkuboyusaku.github.io/poscosea/reference/nmix_weights.md)).

- `"site"`:

  Only sites are resampled (\\(v = 1\\)); equivalent to
  [`ijk_bootstrap()`](https://ohkuboyusaku.github.io/poscosea/reference/ijk_bootstrap.md)
  with the site-level influence.

- `"visit"`:

  Sites are kept (\\(u = 1\\)) and only visits are resampled within each
  site.

For a quantity that involves the latent abundance itself, such as the
realised total \\(\sum_i N_i\\), the likelihood-based replicates only
capture the variation through the parameters. Pass `latent_resid` (see
[`nmix_latent_resid()`](https://ohkuboyusaku.github.io/poscosea/reference/nmix_latent_resid.md))
to add the prediction error of the latent variables.

## See also

[`nmix_loglik()`](https://ohkuboyusaku.github.io/poscosea/reference/nmix_loglik.md),
[`ijk_nmix()`](https://ohkuboyusaku.github.io/poscosea/reference/ijk_nmix.md),
[`nmix_weights()`](https://ohkuboyusaku.github.io/poscosea/reference/nmix_weights.md)

## Examples

``` r
set.seed(1)
I <- 10; J <- 3; S <- 200
y <- matrix(rpois(I * J, 3), I, J)
lambda <- matrix(rgamma(S * I, 20, 4), S, I)
p <- matrix(runif(S * I * J, 0.4, 0.7), S, I * J)
ll <- nmix_loglik(y, lambda, p)
fit <- ijk_nmix(rowSums(lambda), ll)
boot <- ijk_bootstrap_2stage(fit, B = 500)
summary(boot)
#>   estimate        bias       se  lower    upper
#> 1 49.89695 -0.06030114 1.182523 47.621 52.15243
```
