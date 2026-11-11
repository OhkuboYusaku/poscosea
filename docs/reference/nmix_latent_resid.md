# Latent Prediction Error of a Total Abundance

For the realised total \\(\sum\_{i \in A} N_i\\) over a set of sites
\\(A\\), returns \\(\sum\_{i \in A} N_i^{(s)} - \sum\_{i \in A} E\[N_i
\mid y_i, \theta^{(s)}\]\\) for each draw \\(s\\). Passing it to
[`ijk_bootstrap_2stage()`](https://ohkuboyusaku.github.io/poscosea/reference/ijk_bootstrap_2stage.md)
adds the uncertainty of the latent abundance that is not captured by the
likelihood weights.

## Usage

``` r
nmix_latent_resid(N, loglik, sites = NULL)
```

## Arguments

- N:

  Finite non-negative integer matrix (\\(S \times I\\)) of posterior
  draws of \\(N_i\\), with dimensions matching `loglik$condN`.

- loglik:

  An object of class `"nmix_loglik"`.

- sites:

  Optional unique integer site indices in `1:ncol(N)`; defaults to all
  sites.

## Value

Numeric vector of length \\(S\\).
