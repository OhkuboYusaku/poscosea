# Getting Started with {poscosea}

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`poscosea`](https://ohkuboyusaku.github.io/poscosea)`)`

## Overview

`poscosea` provides fast sensitivity analysis for posterior summaries
obtained from MCMC output. The package implements infinitesimal
jackknife (IJK) approximations based on

- posterior draws of a scalar quantity of interest,
- corresponding log-likelihood draws.

These influence values can then be used for

- nonparametric bootstrap approximations,
- leave-k-out cross-validation,
- sensitivity analysis without refitting the original Bayesian model.

## Installation

This package is available on GitHub.

You can install it using `remotes`:

\
[`install.packages`](https://rdrr.io/r/utils/install.packages.html)`(``"remotes"``)`\
`remotes``::`[`install_github`](https://rdrr.io/pkg/remotes/man/install_github.html)`(``"OhkuboYusaku/poscosea"``)`

Or using devtools:

\
[`install.packages`](https://rdrr.io/r/utils/install.packages.html)`(``"devtools"``)`\
`devtools``::`[`install_github`](https://rdrr.io/pkg/devtools/man/install-deprecated.html)`(``"OhkuboYusaku/poscosea"``)`

## Input Format

The main function requires

- a numeric vector of posterior draws (`theta`)
- a log-likelihood matrix (`loglik`)

Rows of `loglik` correspond to the posterior draws in `theta` and
columns correspond to observations.

\
[`set.seed`](https://rdrr.io/r/base/Random.html)`(``123``)`\
`n_obs`` ``<-`` ``6`\
`y`` ``<-`` `[`rnorm`](https://rdrr.io/r/stats/Normal.html)`(``n_obs``, sd ``=`` ``1``)`\
`# Introduce two influential observations.`\
`signal_idx`` ``<-`` `[`c`](https://rdrr.io/r/base/c.html)`(``2``, ``5``)`\
`y``[``signal_idx``]`` ``<-`` ``y``[``signal_idx``]`` ``+`` ``3`\
\
`# Under y_i | theta ~ N(theta, 1) and a flat prior on theta,`\
`# theta | y is Normal(mean(y), 1 / n_obs).`\
`n_draws`` ``<-`` ``2000`\
`theta`` ``<-`` `[`rnorm`](https://rdrr.io/r/stats/Normal.html)`(``n_draws``, mean ``=`` `[`mean`](https://rdrr.io/r/base/mean.html)`(``y``)``, sd ``=`` ``1`` ``/`` `[`sqrt`](https://rdrr.io/r/base/MathFun.html)`(``n_obs``)``)`\
\
`# Pointwise log-likelihood evaluated at the corresponding posterior draws.`\
`loglik`` ``<-`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(``0``, ``n_draws``, ``n_obs``)`\
`for`` ``(``i`` ``in`` `[`seq_len`](https://rdrr.io/r/base/seq.html)`(``n_obs``)``)`` ``{`\
`  ``loglik``[``, ``i``]`` ``<-`` `[`dnorm`](https://rdrr.io/r/stats/Normal.html)`(``y``[``i``]``, mean ``=`` ``theta``, sd ``=`` ``1``, log ``=`` ``TRUE``)`\
`}`

## Constructing an IJK Object

The [`ijk()`](https://ohkuboyusaku.github.io/poscosea/reference/ijk.md)
function computes the posterior mean and observation-wise influence
values.

\
`fit`` ``<-`` `[`ijk`](https://ohkuboyusaku.github.io/poscosea/reference/ijk.md)`(``theta``, ``loglik``)`

The estimated posterior mean is stored in

\
`fit``$``estimate`\
`#> [1] 1.458558`

The observation-wise influence values are stored in

\
[`head`](https://rdrr.io/r/utils/head.html)`(``fit``$``influence``)`\
`#> [1] -0.33446200  0.22035093  0.01858467 -0.22934285  0.28023625  0.04463300`

## Bootstrap Approximation

Bootstrap replicates can be generated from the influence values.

\
`boot`` ``<-`` `[`ijk_bootstrap`](https://ohkuboyusaku.github.io/poscosea/reference/ijk_bootstrap.md)`(`\
`  estimate ``=`` ``fit``$``estimate``,`\
`  influence ``=`` ``fit``$``influence``,`\
`  B ``=`` ``500`\
`)`

Summary statistics can be obtained from the bootstrap distribution.

\
[`print`](https://rdrr.io/r/base/print.html)`(``boot``)`\
`#> `\
`#> IJK Bootstrap Summary`\
`#> `\
`#> Original estimate: 1.459 `\
`#> `\
`#> Bootstrap replicates: 500 `\
`#> `\
`#> Bias: -0.003 `\
`#> Bootstrap SE: 0.521 `\
`#> 95% percentile interval:`\
`#> [0.455, 2.481]`

## Leave-One-Out Approximation

Leave-one-out approximations are obtained by setting `k = 1`.

\
`loo`` ``<-`` `[`ijk_lko`](https://ohkuboyusaku.github.io/poscosea/reference/ijk_lko.md)`(`\
`  estimate ``=`` ``fit``$``estimate``,`\
`  influence ``=`` ``fit``$``influence``,`\
`  k ``=`` ``1`\
`)`

The influence values can be visualized with
[`plot()`](https://rdrr.io/r/graphics/plot.default.html). For
leave-one-out (`k = 1`), each point shows the approximate change in the
estimate when one observation is removed. Large absolute values indicate
highly influential observations.

\
[`plot`](https://rdrr.io/r/graphics/plot.default.html)`(``loo``)`

![](getting-started_files/figure-html/unnamed-chunk-11-1.png)

## Leave-k-Out Approximation

More generally, `k` can be any integer from 1 to the number of
observations, provided the number of omission sets does not exceed 1e6.

\
`lko`` ``<-`` `[`ijk_lko`](https://ohkuboyusaku.github.io/poscosea/reference/ijk_lko.md)`(`\
`  estimate ``=`` ``fit``$``estimate``,`\
`  influence ``=`` ``fit``$``influence``,`\
`  k ``=`` ``2`\
`)`

### Plotting for k \> 1

The current
[`plot.ijk_lko()`](https://ohkuboyusaku.github.io/poscosea/reference/plot.ijk_lko.md)
method is available only for `k = 1`. For `k > 1`, the object contains
all leave-k-out estimates in `lko$lko_estimates`; inspect them directly
or summarize them with `summary(lko$lko_estimates)`.

## N-Mixture Models

For hierarchical count data analysed with a binomial N-mixture model,
see
[`vignette("nmix-twostage", package = "poscosea")`](https://ohkuboyusaku.github.io/poscosea/articles/nmix-twostage.md),
which describes the two-stage (site and visit) bootstrap approximation
implemented in
[`nmix_loglik()`](https://ohkuboyusaku.github.io/poscosea/reference/nmix_loglik.md),
[`ijk_nmix()`](https://ohkuboyusaku.github.io/poscosea/reference/ijk_nmix.md)
and
[`ijk_bootstrap_2stage()`](https://ohkuboyusaku.github.io/poscosea/reference/ijk_bootstrap_2stage.md).
