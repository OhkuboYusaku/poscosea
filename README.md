---
title: "Getting Started with {poscosea}"
author: "Ph.D. Ohkubo Yusaku "
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Getting Started with {poscosea}}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

```{r setup}
library(poscosea)
```

## Overview

`poscosea` provides fast sensitivity analysis for posterior summaries
obtained from MCMC output. The package implements infinitesimal jackknife (IJK) approximations based on

* posterior draws of a scalar quantity of interest,
* corresponding log-likelihood draws.

These influence values can then be used for

* nonparametric bootstrap approximations,
* leave-k-out cross-validation,
* sensitivity analysis without refitting the original Bayesian model.

## Installation

This package is available on GitHub.

You can install it using `remotes`:

```{r eval=FALSE}
install.packages("remotes")
remotes::install_github("OhkuboYusaku/poscosea")
```
Or using devtools:

```{r eval=FALSE}
install.packages("devtools")
devtools::install_github("OhkuboYusaku/poscosea")
```
## Input Format

The main function requires

* a numeric vector of posterior draws (`theta`)
* a log-likelihood matrix (`loglik`)

Rows of `loglik` correspond to MCMC draws and columns correspond to observations.

```{r}
n_obs <- 6
y <- rnorm(n_obs, sd = 1) # data
# --- Introducing a few influential observations ---
signal_idx <- c(2)
y[signal_idx]<- y[signal_idx] + 3
n_draws <- 30
theta <- rnorm(n_draws)
# log-likelihood matrix
loglik <- matrix(0, n_draws, n_obs)
for (i in seq_len(n_obs)){ 
  # likelihood: y_i ~ N(theta + beta_i, sigma^2)
  loglik[, i] <- dnorm(y[i], mean = theta, sd = 1, log = TRUE)
}
```

## Constructing an IJK Object

The `ijk()` function computes the posterior mean and observation-wise
influence values.

```{r}
fit <- ijk(theta, loglik)

fit
```

The estimated posterior mean is stored in

```{r}
fit$estimate
```

The observation-wise influence values are stored in

```{r}
head(fit$influence)
```

## Bootstrap Approximation

Bootstrap replicates can be generated from the influence values.

```{r}
boot <- ijk_bootstrap(
  estimate = fit$estimate,
  influence = fit$influence,
  B = 500
)

```

Summary statistic can be obtained from the bootstrap distribution.

```{r}
print(boot)
```

## Leave-One-Out Approximation

Leave-one-out approximations are obtained by setting `k = 1`.

```{r}
loo <- ijk_lko(
  estimate = fit$estimate,
  influence = fit$influence,
  k = 1
)
```

The influence structure can be visualized using plot(). For leave-one-out (k = 1), each point represents the change in the estimate when a single observation is removed. Large absolute values indicate highly influential observations.

```{r}
plot(loo)
```

# Leave-k-Out Approximation

More generally, any value of `k` can be used.

```{r}
lko <- ijk_lko(
  estimate = fit$estimate,
  influence = fit$influence,
  k = 2
)
```

Notes on plotting

The current implementation of plot.ijk_lko() is available only for k = 1. For k > 1, the object still contains all leave-k-out estimates in lko$lko_estimates, which can be summarized using summary() or manually inspected.
