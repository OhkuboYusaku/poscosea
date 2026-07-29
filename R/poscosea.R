#' Infinitesimal Jackknife Approximation for MCMC Output
#'
#' Computes an infinitesimal jackknife (IJK) approximation of the
#' posterior mean from MCMC draws and corresponding log-likelihood values.
#'
#' The function estimates observation-wise influence values that can be used
#' for bootstrap and leave-k-out cross-validation approximations without
#' refitting the original model.
#'
#' @param theta Numeric vector of posterior draws of a scalar quantity.
#' @param loglik Numeric matrix of log-likelihood draws.
#' Rows correspond to posterior draws and columns correspond to observations.
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
#' @param theta Numeric vector.
#' @param loglik Numeric matrix.
#'
#' @return Invisibly returns TRUE.
#'
#' @keywords internal
check_inputs <- function(theta, loglik) {
  if (!is.numeric(theta))
    stop("'theta' must be numeric.", call. = FALSE)
  if (anyNA(theta))
    stop("'theta' contains missing values.", call. = FALSE)
  if (!is.matrix(loglik))
    stop("'loglik' must be a matrix.", call. = FALSE)
  if (!is.numeric(loglik))
    stop("'loglik' must be numeric.", call. = FALSE)
  if (anyNA(loglik))
    stop("'loglik' contains missing values.", call. = FALSE)
  if (length(theta) != nrow(loglik))
    stop(
      "length(theta) must equal nrow(loglik).",
      call. = FALSE
    )
  invisible(TRUE)
}


#' Compute Infinitesimal Jackknife Influence Values
#'
#' Computes observation-wise influence values used in the
#' infinitesimal jackknife approximation.
#'
#' @param theta Numeric vector of posterior draws.
#' @param loglik Numeric matrix of log-likelihood draws.
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
#' @param estimate Numeric scalar giving the original estimate.
#' @param influence Numeric vector of influence values.
#' @param B Positive integer giving the number of bootstrap replicates.
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

  if (!is.numeric(estimate) || length(estimate) != 1)
    stop(
      "'estimate' must be a numeric scalar.",
      call. = FALSE
    )

  if (!is.numeric(influence))
    stop(
      "'influence' must be numeric.",
      call. = FALSE
    )

  if (!is.numeric(B) || length(B) != 1 || B < 1)
    stop(
      "'B' must be a positive integer.",
      call. = FALSE
    )

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
#' @param ... Additional arguments passed to methods.
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
#' @param ... Additional arguments passed to methods.
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
#' @param estimate Scalar estimate.
#' @param influence Numeric vector of influence values.
#' @param k Number of observations removed.
#'
#' @return
#' Numeric vector of leave-k-out estimates.
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
  if (!is.numeric(estimate) || length(estimate) != 1)
    stop("'estimate' must be a numeric scalar.",
         call. = FALSE)
  if (!is.numeric(influence))
    stop("'influence' must be numeric.",
         call. = FALSE)
  R <- length(influence)

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
#' @param label Integer giving the number of most influential
#'   observations to label. Default is \code{5}.
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
#' summary(lko1)
#' plot(lko1)
#'
#' # Leave-two-out (k = 2)
#' lko2 <- ijk_lko(
#'   estimate = 0,
#'   influence = influence,
#'   k = 2
#' )
#'
#' summary(lko2)
#'
#' @export
#' @method plot ijk_lko
plot.ijk_lko <- function(
    x,
    label = 5,
    pch = 19,
    ...
) {
  if (x$k != 1) {
    stop(
      "plot.ijk_lko() is currently implemented only for k = 1.",
      call. = FALSE
    )
  }
  delta <- x$delta
  n <- length(delta)

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
