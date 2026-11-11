#' Infinitesimal Jackknife Approximation for MCMC Output
#'
#' Computes an infinitesimal jackknife (IJK) approximation of the
#' posterior mean from MCMC draws and corresponding log-likelihood values.
#'
#' The function estimates observation-wise influence values that can be used
#' for bootstrap and leave-k-out cross-validation approximations without
#' refitting the original model.
#'
#' @param theta Finite numeric vector of posterior draws of a scalar
#'   quantity, with at least two draws.
#' @param loglik Finite numeric matrix of log-likelihood draws. Rows correspond
#'   to posterior draws and columns correspond to observations; it must have at
#'   least one column and one row per value in `theta`.
#'
#' @return
#' An object of class `"ijk"` containing:
#' \describe{
#'   \item{estimate}{Posterior mean estimate.}
#'   \item{influence}{Observation-wise influence values.}
#'   \item{theta}{Input posterior draws.}
#'   \item{loglik}{Input log-likelihood matrix.}
#' }
#'
#' @examples
#' set.seed(123)
#' theta <- rnorm(1000)
#' loglik <- matrix(rnorm(1000 * 20), 1000, 20)
#'
#' fit <- ijk(theta, loglik)
#' fit$estimate
#'
#' @export
ijk <- function(theta, loglik) {
  check_inputs(theta, loglik)
  estimate <- mean(theta)
  influence <- compute_if(theta, loglik)
  structure(
    list(
      estimate = estimate,
      influence = influence,
      theta = theta,
      loglik = loglik
    ),
    class = "ijk"
  )
}


#' Validate Inputs
#'
#' Internal utility for validating inputs supplied to
#' the main ijk function.
#'
#' @param theta Finite numeric vector with at least two values.
#' @param loglik Finite numeric matrix with at least one column and one row
#'   per value in `theta`.
#'
#' @return Invisibly returns TRUE.
#'
#' @keywords internal
check_inputs <- function(theta, loglik) {
  if (!is.numeric(theta) || !is.null(dim(theta)) || length(theta) < 2L)
    stop("'theta' must be a numeric vector with at least two draws.",
         call. = FALSE)
  if (any(!is.finite(theta)))
    stop("'theta' must contain only finite values.", call. = FALSE)
  if (!is.matrix(loglik) || !is.numeric(loglik))
    stop("'loglik' must be a numeric matrix.", call. = FALSE)
  if (ncol(loglik) < 1L)
    stop("'loglik' must contain at least one observation.", call. = FALSE)
  if (any(!is.finite(loglik)))
    stop("'loglik' must contain only finite values.", call. = FALSE)
  if (length(theta) != nrow(loglik))
    stop("length(theta) must equal nrow(loglik).", call. = FALSE)
  invisible(TRUE)
}


#' Compute Infinitesimal Jackknife Influence Values
#'
#' Computes observation-wise influence values used in the
#' infinitesimal jackknife approximation.
#'
#' @param theta Finite numeric vector of at least two posterior draws.
#' @param loglik Finite numeric matrix with one row per draw in `theta` and
#'   at least one observation column.
#'
#' @return
#' A numeric vector of influence values with one element per observation.
#'
#' @keywords internal
compute_if <- function(theta, loglik) {
  check_inputs(theta, loglik)
  if_val <- apply(loglik, 2, function(x) stats::cov(theta, x))

  if_val - mean(if_val)
}


#' Bootstrap Using Infinitesimal Jackknife Approximation
#'
#' Generates bootstrap replicates using infinitesimal jackknife
#' (IJK) influence values.
#'
#' The bootstrap distribution is obtained by resampling observation
#' weights from a multinomial distribution and applying the linear
#' influence-function approximation.
#'
#' @param estimate Finite numeric scalar giving the original estimate.
#' @param influence Non-empty finite numeric vector of influence values.
#' @param B Integer of at least 2 giving the number of bootstrap replicates.
#'
#' @return
#' An object of class "ijk_boot" containing:
#' \describe{
#'   \item{estimate}{Original estimate.}
#'   \item{replicates}{Bootstrap replicates.}
#'   \item{B}{Number of bootstrap replicates.}
#' }
#'
#' @examples
#' set.seed(123)
#'
#' influence <- rnorm(20)
#'
#' boot <- ijk_bootstrap(
#'   estimate = 0,
#'   influence = influence,
#'   B = 200
#' )
#'
#' boot
#' summary(boot)
#'
#' @export
ijk_bootstrap <- function(
    estimate,
    influence,
    B = 1000
) {

  if (!is.numeric(estimate) || !is.null(dim(estimate)) ||
      length(estimate) != 1L || !is.finite(estimate))
    stop("'estimate' must be a finite numeric scalar.", call. = FALSE)

  if (!is.numeric(influence) || !is.null(dim(influence)) ||
      length(influence) < 1L || any(!is.finite(influence)))
    stop("'influence' must be a non-empty finite numeric vector.",
         call. = FALSE)

  if (!is.numeric(B) || !is.null(dim(B)) || length(B) != 1L ||
      !is.finite(B) || B < 2 || B != floor(B))
    stop("'B' must be an integer of at least 2.", call. = FALSE)

  R <- length(influence)

  reps <- replicate(B, {

    w <- stats::rmultinom(
      1,
      size = R,
      prob = rep(1 / R, R)
    )[, 1]

    estimate + sum((w - 1) * influence)
  })

  structure(
    list(
      estimate = estimate,
      replicates = reps,
      B = B
    ),
    class = "ijk_boot"
  )
}


#' Print an IJK Bootstrap Object
#'
#' Displays a concise summary of an object of class "ijk_boot".
#'
#' @param x An object of class "ijk_boot".
#' @param digits Number of digits to display.
#' @param ... Currently unused.
#'
#' @return
#' The input object, invisibly.
#'
#' @examples
#' set.seed(123)
#' boot <- ijk_bootstrap(
#'   estimate = 0,
#'   influence = rnorm(20),
#'   B = 100
#' )
#'
#' print(boot)
#'
#' @export
print.ijk_boot <- function(
    x,
    digits = 3,
    ...
) {
  reps <- x$replicates
  bias <- mean(reps) - x$estimate
  se <- stats::sd(reps)
  ci <- stats::quantile(
    reps,
    probs = c(0.025, 0.975)
  )
  cat("\n")
  cat("IJK Bootstrap Summary\n\n")
  cat(
    "Original estimate:",
    round(x$estimate, digits),
    "\n\n"
  )
  cat(
    "Bootstrap replicates:",
    length(reps),
    "\n\n"
  )
  cat(
    "Bias:",
    round(bias, digits),
    "\n"
  )
  cat(
    "Bootstrap SE:",
    round(se, digits),
    "\n"
  )
  cat(
    "95% percentile interval:\n"
  )
  cat(
    "[",
    round(ci[1], digits),
    ", ",
    round(ci[2], digits),
    "]\n",
    sep = ""
  )
  invisible(x)
}


#' Summarize an IJK Bootstrap Object
#'
#' Computes standard bootstrap summary statistics.
#'
#' @param object An object of class "ijk_boot".
#' @param ... Currently unused.
#'
#' @return
#' A list containing:
#' \describe{
#'   \item{estimate}{Original estimate.}
#'   \item{bias}{Bootstrap bias estimate.}
#'   \item{se}{Bootstrap standard error.}
#'   \item{ci}{95% percentile confidence interval.}
#' }
#'
#' @examples
#' set.seed(123)
#' boot <- ijk_bootstrap(
#'   estimate = 0,
#'   influence = rnorm(20),
#'   B = 100
#' )
#'
#' summary(boot)
#'
#' @export
summary.ijk_boot <- function(
    object,
    ...
) {
  reps <- object$replicates
  list(
    estimate = object$estimate,
    bias = mean(reps) - object$estimate,
    se = stats::sd(reps),
    ci = stats::quantile(
      reps,
      probs = c(0.025, 0.975)
    )
  )
}

#' Leave-k-Out Approximation
#'
#' Approximates leave-k-out cross-validation estimates
#' using infinitesimal jackknife influence values.
#'
#' @param estimate Finite numeric scalar giving the original estimate.
#' @param influence Non-empty finite numeric vector of influence values.
#' @param k Finite integer number of observations removed, satisfying
#'   `1 <= k <= length(influence)`. Exhaustive enumeration is limited to at most
#'   1e6 omission sets.
#'
#' @return
#' An object of class `ijk_lko` containing:
#' \describe{
#'   \item{estimate}{Original estimate.}
#'   \item{lko_estimates}{Approximate estimates for all omitted sets.}
#'   \item{delta}{Change from the original estimate for each omitted set.}
#'   \item{omitted}{List of omitted observation indices for each estimate.}
#'   \item{k}{Number of observations omitted in each set.}
#' }
#'
#' @examples
#' influence <- rnorm(20)
#' ijk_lko(0, influence, k = 1)
#'
#' @export
ijk_lko <- function(
    estimate,
    influence,
    k = 1
) {
  if (!is.numeric(estimate) || !is.null(dim(estimate)) ||
      length(estimate) != 1L || !is.finite(estimate))
    stop("'estimate' must be a finite numeric scalar.", call. = FALSE)
  if (!is.numeric(influence) || !is.null(dim(influence)) ||
      length(influence) < 1L || any(!is.finite(influence)))
    stop("'influence' must be a non-empty finite numeric vector.",
         call. = FALSE)
  R <- length(influence)

  if (!is.numeric(k) || !is.null(dim(k)) || length(k) != 1L ||
      !is.finite(k) || k != floor(k))
    stop("'k' must be a finite integer.", call. = FALSE)
  if (k < 1 || k > R)
    stop(
      "'k' must satisfy 1 <= k <= length(influence).",
      call. = FALSE
    )

  ncomb <- choose(R, k)

  if (ncomb > 1e6) {
    stop(
      "choose(R, k) = ",
      format(ncomb, scientific = TRUE),
      " is too large for exhaustive leave-k-out.",
      call. = FALSE
    )
  }

  omit_sets <- utils::combn(
    R,
    k,
    simplify = FALSE
  )

  est <- vapply(
    omit_sets,
    function(idx) {
      estimate - sum(influence[idx])
    },
    numeric(1)
  )

  structure(
    list(
      estimate = estimate,
      lko_estimates = est,
      delta = est - estimate,
      omitted = omit_sets,
      k = k
    ),
    class = "ijk_lko"
  )
}

#' Plot Leave-k-out Influence Diagnostics
#'
#' Produces an influence diagnostic plot for an object of class
#' \code{"ijk_lko"}.
#'
#' The vertical axis shows the change in the posterior estimate
#' after removing observations:
#'
#' \deqn{
#' \Delta_i = \hat{\theta}_{(-i)} - \hat{\theta}.
#' }
#'
#' For leave-one-out (\code{k = 1}), the horizontal axis corresponds
#' to observation IDs. Optionally, the most influential observations
#' can be labelled.
#'
#' @param x An object of class \code{"ijk_lko"}.
#' @param label Integer from 0 to the number of observations, giving the
#'   number of most influential observations to label. Default is \code{5}.
#' @param pch Plotting character.
#' @param ... Additional graphical parameters passed to
#'   \code{\link[graphics]{plot}}.
#'
#' @return
#' Produces a diagnostic plot and returns \code{x} invisibly.
#'
#' @examples
#' set.seed(123)
#'
#' influence <- rnorm(20)
#'
#' # Leave-one-out (k = 1)
#' lko1 <- ijk_lko(
#'   estimate = 0,
#'   influence = influence,
#'   k = 1
#' )
#'
#' plot(lko1)
#'
#' @export
#' @method plot ijk_lko
plot.ijk_lko <- function(
    x,
    label = 5,
    pch = 19,
    ...
) {
  if (!inherits(x, "ijk_lko") || !is.numeric(x$k) || length(x$k) != 1L ||
      !is.finite(x$k) || !is.numeric(x$delta) || length(x$delta) < 1L ||
      any(!is.finite(x$delta)))
    stop("'x' must be a valid object returned by ijk_lko().", call. = FALSE)
  if (x$k != 1) {
    stop(
      "plot.ijk_lko() is currently implemented only for k = 1.",
      call. = FALSE
    )
  }
  delta <- x$delta
  n <- length(delta)
  if (!is.numeric(label) || !is.null(dim(label)) || length(label) != 1L ||
      !is.finite(label) || label < 0 || label > n || label != floor(label))
    stop("'label' must be an integer between 0 and length(x$delta).",
         call. = FALSE)

  graphics::plot(
    seq_len(n),
    delta,
    pch = pch,
    xlab = "Observation ID",
    ylab = expression(
      Delta[i]
    ),
    main = "Leave-one-out influence diagnostics",
    ...
  )

  graphics::abline(
    h = 0,
    lty = 2
  )

  top <- order(
    abs(delta),
    decreasing = TRUE
  )[seq_len(min(label, n))]
  pos_vec <- ifelse(delta[top] > 0, 1, 3)
  graphics::text(
    x = top,
    y = delta[top],
    labels = paste0("ID=", top),
    pos = pos_vec
  )

  invisible(x)
}
