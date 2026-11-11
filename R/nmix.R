#' Log-Likelihood Components of an N-Mixture Model
#'
#' Computes, for every posterior draw, the site-level marginal
#' log-likelihood and the observation (visit)-level log-likelihood of a
#' binomial N-mixture model
#' \deqn{N_i \sim \mathrm{Poisson}(\lambda_i), \quad
#'       y_{ij} \mid N_i \sim \mathrm{Binomial}(N_i, p_{ij}),}
#' for sites \eqn{i = 1, \ldots, I} and visits \eqn{j = 1, \ldots, J}.
#' These are the inputs of [ijk_nmix()] and the two-stage bootstrap
#' [ijk_bootstrap_2stage()].
#'
#' The latent abundance \eqn{N_i} is summed out over the grid
#' \eqn{N = \max_j y_{ij}, \ldots, \max_j y_{ij} + } `Nmax_extra`.
#' Missing visits (`NA` in `y`) are dropped from the likelihood.
#'
#' The observation-level log-likelihood is returned in two versions:
#' \describe{
#'   \item{`obs_rb`}{\eqn{E[\log \mathrm{Bin}(y_{ij} \mid N, p_{ij}) \mid y_i, \theta]},
#'     the expectation over the conditional posterior of \eqn{N_i}.
#'     Use it for quantities that depend on the parameters only
#'     (e.g. regression coefficients, \eqn{\lambda}).
#'     It has smaller Monte Carlo noise.}
#'   \item{`obs_joint`}{\eqn{\log \mathrm{Bin}(y_{ij} \mid N_i^{(s)}, p_{ij}^{(s)})}
#'     evaluated at the posterior draws of \eqn{N_i}. Use it for quantities
#'     that depend on the latent abundance (e.g. \eqn{\sum_i N_i}).
#'     Only available when `N` is supplied.}
#' }
#'
#' @param y Numeric counts with \eqn{I} rows (sites) and \eqn{J} columns
#'   (visits); a numeric vector is treated as a one-column matrix. Counts must
#'   be finite non-negative integers, except that `NA` marks a missing visit.
#'   Every site must have at least one observed visit.
#' @param lambda Finite positive numeric matrix (\eqn{S \times I}) of
#'   posterior draws of \eqn{\lambda_i}; at least one draw is required.
#' @param p Posterior draws of \eqn{p_{ij}}: either an array
#'   (\eqn{S \times I \times J}) or a matrix (\eqn{S \times IJ}) whose
#'   columns are ordered as \eqn{k = i + (j - 1) I} (the column-major order
#'   used by JAGS for `p[i,j]`). Values for observed visits must be finite and
#'   between 0 and 1; values for missing visits are not used.
#' @param N Optional numeric matrix (\eqn{S \times I}) of finite
#'   non-negative integer posterior draws of \eqn{N_i}; each draw must be at
#'   least the largest observed count at its site. Required for `obs_joint`.
#' @param Nmax_extra Non-negative integer; the summation grid for
#'   \eqn{N_i} extends `Nmax_extra` above \eqn{\max_j y_{ij}}.
#'
#' @return An object of class `"nmix_loglik"` containing
#' \describe{
#'   \item{site}{\eqn{S \times I} matrix of site-level marginal log-likelihoods.}
#'   \item{obs_rb}{\eqn{S \times IJ} matrix (see Details). Missing visits are 0.}
#'   \item{obs_joint}{\eqn{S \times IJ} matrix or `NULL`.}
#'   \item{condN}{\eqn{S \times I} matrix of \eqn{E[N_i \mid y_i, \theta]}.}
#'   \item{obs_mask}{Logical \eqn{I \times J} matrix, `TRUE` for observed visits.}
#'   \item{tail_mass}{Largest conditional probability at the upper end of the
#'     grid; a warning is issued when it exceeds `1e-8`.}
#' }
#'
#' @seealso [nmix_draws()], [ijk_nmix()], [ijk_bootstrap_2stage()]
#'
#' @examples
#' set.seed(1)
#' I <- 10; J <- 3; S <- 50
#' y <- matrix(rpois(I * J, 3), I, J)
#' lambda <- matrix(rgamma(S * I, 20, 4), S, I)
#' p <- matrix(runif(S * I * J, 0.4, 0.7), S, I * J)
#' ll <- nmix_loglik(y, lambda, p)
#' dim(ll$site)
#'
#' @export
nmix_loglik <- function(y, lambda, p, N = NULL, Nmax_extra = 100) {
  y <- as.matrix(y)
  if (!is.numeric(y) || length(dim(y)) != 2L ||
      nrow(y) < 1L || ncol(y) < 1L)
    stop("'y' must be a non-empty numeric matrix.", call. = FALSE)
  if (any(is.nan(y)) || any(!is.na(y) & !is.finite(y)))
    stop("'y' may contain NA for missing visits, but no NaN or Inf.",
         call. = FALSE)
  y_obs <- y[!is.na(y)]
  if (any(y_obs < 0) || any(y_obs != floor(y_obs)))
    stop("'y' must contain non-negative integers or NA.", call. = FALSE)
  I <- nrow(y)
  J <- ncol(y)
  obs_mask <- !is.na(y)
  if (any(rowSums(obs_mask) == 0))
    stop("Every site must have at least one observed visit.", call. = FALSE)

  lambda <- as.matrix(lambda)
  if (!is.numeric(lambda) || nrow(lambda) < 1L || ncol(lambda) != I)
    stop("'lambda' must be a numeric S x I matrix with at least one draw.",
         call. = FALSE)
  S <- nrow(lambda)
  if (any(!is.finite(lambda)) || any(lambda <= 0))
    stop("'lambda' must contain finite positive values.", call. = FALSE)

  if (length(dim(p)) == 3L) {
    if (!is.numeric(p) || !all(dim(p) == c(S, I, J)))
      stop("A 3-d 'p' must be numeric with dimension c(S, I, J).",
           call. = FALSE)
    p <- matrix(p, S, I * J)
  } else if (!is.matrix(p)) {
    stop("'p' must be an S x (I*J) matrix or an S x I x J array.",
         call. = FALSE)
  }
  if (!is.numeric(p) || !all(dim(p) == c(S, I * J)))
    stop("'p' must be a numeric S x (I*J) matrix or S x I x J array.",
         call. = FALSE)
  p_obs <- p[, as.vector(obs_mask), drop = FALSE]
  if (any(!is.finite(p_obs)) || any(p_obs < 0 | p_obs > 1))
    stop("'p' must be finite and between 0 and 1 for observed visits.",
         call. = FALSE)

  if (!is.null(N)) {
    N <- as.matrix(N)
    if (!is.numeric(N) || !all(dim(N) == c(S, I)))
      stop("'N' must be a numeric S x I matrix.", call. = FALSE)
    if (any(!is.finite(N)) || any(N < 0) || any(N != floor(N)))
      stop("'N' must contain finite non-negative integers.", call. = FALSE)
    for (i in seq_len(I)) {
      if (any(N[, i] < max(y[i, obs_mask[i, ]])))
        stop("'N' draws must be at least the observed count at each site.",
             call. = FALSE)
    }
  }
  if (!is.numeric(Nmax_extra) || !is.null(dim(Nmax_extra)) ||
      length(Nmax_extra) != 1L || !is.finite(Nmax_extra) ||
      Nmax_extra < 0 || Nmax_extra != floor(Nmax_extra))
    stop("'Nmax_extra' must be a non-negative integer.", call. = FALSE)

  site <- condN <- matrix(NA_real_, S, I)
  obs_rb <- matrix(0, S, I * J)
  obs_joint <- if (is.null(N)) NULL else matrix(0, S, I * J)
  tail_mass <- 0
  n_zero <- 0

  for (i in seq_len(I)) {
    obs <- which(obs_mask[i, ])
    ymax <- max(y[i, obs])
    Ns <- seq.int(ymax, ymax + Nmax_extra)
    Nrow <- matrix(Ns, S, length(Ns), byrow = TRUE)
    lt <- stats::dpois(Nrow, lambda[, i], log = TRUE)
    lb <- vector("list", J)
    for (j in obs) {
      k <- i + (j - 1L) * I
      lb[[j]] <- stats::dbinom(y[i, j], Nrow, p[, k], log = TRUE)
      lt <- lt + lb[[j]]
      if (!is.null(N))
        obs_joint[, k] <- stats::dbinom(y[i, j], N[, i], p[, k], log = TRUE)
    }
    m <- apply(lt, 1, max)
    zero <- !is.finite(m)            # draws with zero likelihood on the grid
    m[zero] <- 0
    w <- exp(lt - m)
    sw <- rowSums(w)
    site[, i] <- m + log(sw)         # -Inf for those draws
    w <- w / sw
    condN[, i] <- drop(w %*% Ns)
    for (j in obs) obs_rb[, i + (j - 1L) * I] <- rowSums(w * lb[[j]])
    tail_mass <- max(tail_mass, w[, ncol(w)], na.rm = TRUE)
    n_zero <- n_zero + sum(zero)
  }

  if (n_zero > 0)
    warning(n_zero, " site x draw combination(s) have zero likelihood ",
            "(site log-likelihood -Inf).", call. = FALSE)
  if (tail_mass > 1e-8)
    warning("Conditional probability at the upper end of the N grid is ",
            format(tail_mass, digits = 3), "; consider increasing 'Nmax_extra'.",
            call. = FALSE)
  if (!is.null(obs_joint) && any(!is.finite(obs_joint)))
    warning("Non-finite values in 'obs_joint' (N draws below the observed counts?).",
            call. = FALSE)

  structure(
    list(site = site, obs_rb = obs_rb, obs_joint = obs_joint, condN = condN,
         obs_mask = obs_mask, tail_mass = tail_mass),
    class = "nmix_loglik"
  )
}


#' Extract N-Mixture Quantities from MCMC Draws
#'
#' Picks the posterior draws of \eqn{\lambda_i}, \eqn{p_{ij}} and
#' (optionally) \eqn{N_i} out of a matrix of MCMC draws whose columns are
#' named as in JAGS/BUGS output (`lambda[1]`, `p[1,1]`, `N[1]`, ...).
#'
#' @param draws A non-empty numeric matrix or data frame of finite draws
#'   with unique, named columns, or a non-empty list of matrices/data frames
#'   with the same columns (chains are stacked).
#' @param I Positive integer giving the number of sites.
#' @param J Positive integer giving the number of visits.
#' @param lambda,p,N Non-empty node-name prefixes for \eqn{\lambda}, \eqn{p}
#'   and \eqn{N}. Set `N = NULL` to skip \eqn{N}. If not all requested `N`
#'   columns are present, the returned `N` element is `NULL`.
#'
#' @return A list with `lambda` (\eqn{S \times I}), `p`
#'   (\eqn{S \times IJ}, columns ordered as \eqn{i + (j-1)I}) and `N`
#'   (\eqn{S \times I} or `NULL` when the columns are absent).
#'
#' @examples
#' I <- 2; J <- 2; S <- 3
#' draws <- matrix(runif(S * (I + I * J)), S, I + I * J)
#' colnames(draws) <- c(paste0("lambda[", 1:I, "]"),
#'                      paste0("p[", rep(1:I, J), ",", rep(1:J, each = I), "]"))
#' str(nmix_draws(draws, I, J))
#'
#' @export
nmix_draws <- function(draws, I, J, lambda = "lambda", p = "p", N = "N") {
  valid_index <- function(x) {
    is.numeric(x) && is.null(dim(x)) && length(x) == 1L &&
      is.finite(x) && x >= 1 && x == floor(x)
  }
  if (!valid_index(I) || !valid_index(J))
    stop("'I' and 'J' must be positive integers.", call. = FALSE)
  if (!is.character(lambda) || length(lambda) != 1L || is.na(lambda) ||
      !nzchar(lambda) || !is.character(p) || length(p) != 1L || is.na(p) ||
      !nzchar(p))
    stop("'lambda' and 'p' must be non-empty node-name prefixes.",
         call. = FALSE)
  if (!is.null(N) && (!is.character(N) || length(N) != 1L || is.na(N) ||
                      !nzchar(N)))
    stop("'N' must be NULL or a non-empty node-name prefix.", call. = FALSE)
  if (is.list(draws) && !is.data.frame(draws)) {
    if (length(draws) == 0L || any(!vapply(draws, function(x) {
      (is.matrix(x) || is.data.frame(x)) && is.numeric(as.matrix(x))
    }, logical(1))))
      stop("'draws' must be a matrix/data frame or a non-empty list of them.",
           call. = FALSE)
    draws <- do.call(rbind, lapply(draws, as.matrix))
  } else if (!is.matrix(draws) && !is.data.frame(draws)) {
    stop("'draws' must be a matrix, data frame, or list of chains.",
         call. = FALSE)
  }
  draws <- as.matrix(draws)
  if (!is.numeric(draws) || nrow(draws) < 1L || ncol(draws) < 1L ||
      any(!is.finite(draws)))
    stop("'draws' must be a non-empty numeric matrix of finite values.",
         call. = FALSE)
  cn <- colnames(draws)
  if (is.null(cn) || anyNA(cn) || any(!nzchar(cn)) || anyDuplicated(cn))
    stop("'draws' must have unique, non-empty column names.", call. = FALSE)

  lam_cols <- paste0(lambda, "[", seq_len(I), "]")
  p_cols <- paste0(p, "[", rep(seq_len(I), J), ",", rep(seq_len(J), each = I), "]")
  miss <- setdiff(c(lam_cols, p_cols), cn)
  if (length(miss) > 0)
    stop("Columns not found in 'draws': ",
         paste(utils::head(miss, 5), collapse = ", "),
         if (length(miss) > 5) ", ...", call. = FALSE)

  N_out <- NULL
  if (!is.null(N)) {
    N_cols <- paste0(N, "[", seq_len(I), "]")
    if (all(N_cols %in% cn)) N_out <- draws[, N_cols, drop = FALSE]
  }
  list(lambda = draws[, lam_cols, drop = FALSE],
       p = draws[, p_cols, drop = FALSE],
       N = N_out)
}


#' Infinitesimal Jackknife Influence for an N-Mixture Model
#'
#' Computes site-level and observation (visit)-level influence values,
#' i.e. posterior covariances between the quantity of interest and the
#' log-likelihood components returned by [nmix_loglik()]. These are used by
#' [ijk_bootstrap_2stage()] for the two-stage (site and visit) bootstrap
#' approximation.
#'
#' The first-order approximation of the posterior mean under site weights
#' \eqn{u_i} and visit weights \eqn{v_{ij}} is
#' \deqn{E_{u,v}[\theta] \approx E[\theta] + \sum_i (u_i - 1) c_i +
#'       \sum_{i,j} (v_{ij} - 1) e_{ij},}
#' with \eqn{c_i = \mathrm{Cov}[\theta, \ell_i]} (site-level marginal
#' log-likelihood) and \eqn{e_{ij} = \mathrm{Cov}[\theta, \ell_{ij}]}
#' (observation-level log-likelihood).
#'
#' @param theta Finite numeric vector (length \eqn{S}, with \eqn{S \ge 2})
#'   of posterior draws of the quantity of interest, or a finite numeric
#'   \eqn{S \times K} matrix for \eqn{K} quantities at once (e.g. a
#'   regression curve on a grid). The number of rows must match `loglik`.
#' @param loglik An object of class `"nmix_loglik"`.
#' @param obs Which observation-level log-likelihood to use: `"rb"` for
#'   quantities depending only on the parameters, `"joint"` for quantities
#'   that depend on the latent abundance \eqn{N}.
#'
#' @return An object of class `"ijk_nmix"` containing
#' \describe{
#'   \item{estimate}{Posterior mean(s).}
#'   \item{influence_site}{Site-level influence values (centred, as in
#'     [ijk()]); a vector of length \eqn{I}, or an \eqn{I \times K} matrix.}
#'   \item{influence_obs}{Observation-level influence values; length
#'     \eqn{IJ} (ordered \eqn{i + (j-1)I}), 0 for missing visits.}
#'   \item{obs_mask}{Observed-visit indicator.}
#'   \item{obs}{The `obs` argument.}
#' }
#'
#' The site-level influence can be passed to [ijk_lko()] for an
#' approximate leave-one-site-out analysis.
#'
#' @seealso [nmix_loglik()], [ijk_bootstrap_2stage()]
#'
#' @examples
#' set.seed(1)
#' I <- 10; J <- 3; S <- 200
#' y <- matrix(rpois(I * J, 3), I, J)
#' lambda <- matrix(rgamma(S * I, 20, 4), S, I)
#' p <- matrix(runif(S * I * J, 0.4, 0.7), S, I * J)
#' ll <- nmix_loglik(y, lambda, p)
#' fit <- ijk_nmix(rowSums(lambda), ll)
#' fit$estimate
#'
#' @export
ijk_nmix <- function(theta, loglik, obs = c("rb", "joint")) {
  if (!inherits(loglik, "nmix_loglik"))
    stop("'loglik' must be created by nmix_loglik().", call. = FALSE)
  obs <- match.arg(obs)
  theta_m <- as.matrix(theta)
  if (!is.numeric(theta_m) || nrow(theta_m) < 2L || ncol(theta_m) < 1L)
    stop("'theta' must be a numeric vector or matrix with at least two draws.",
         call. = FALSE)
  if (any(!is.finite(theta_m)))
    stop("'theta' must contain only finite values.", call. = FALSE)
  if (!is.matrix(loglik$site) || !is.numeric(loglik$site) ||
      ncol(loglik$site) < 1L)
    stop("'loglik' must contain a non-empty numeric 'site' matrix.",
         call. = FALSE)
  S <- nrow(loglik$site)
  if (nrow(theta_m) != S)
    stop("The number of draws in 'theta' must equal nrow(loglik$site).", call. = FALSE)

  if (!is.matrix(loglik$obs_mask) || !is.logical(loglik$obs_mask) ||
      nrow(loglik$obs_mask) < 1L || ncol(loglik$obs_mask) < 1L ||
      anyNA(loglik$obs_mask) || nrow(loglik$obs_mask) != ncol(loglik$site))
    stop("'loglik' contains an invalid observed-visit mask.", call. = FALSE)
  L_obs <- if (obs == "joint") loglik$obs_joint else loglik$obs_rb
  if (is.null(L_obs))
    stop("obs = 'joint' requires 'N' in nmix_loglik().", call. = FALSE)
  if (!is.matrix(L_obs) || !is.numeric(L_obs) ||
      !all(dim(L_obs) == c(S, length(loglik$obs_mask))))
    stop("'loglik' contains an observation log-likelihood matrix with invalid dimensions.",
         call. = FALSE)
  if (any(!is.finite(loglik$site)) || any(!is.finite(L_obs)))
    stop("The log-likelihood contains non-finite values; check the draws ",
         "and 'Nmax_extra'.", call. = FALSE)

  est <- colMeans(theta_m)
  tc <- sweep(theta_m, 2, est)
  # Cov[theta, L] = crossprod(L, theta - mean) / (S - 1)  (no need to centre L)
  c_site <- crossprod(loglik$site, tc) / (S - 1)
  c_site <- sweep(c_site, 2, colMeans(c_site))
  e_obs <- crossprod(L_obs, tc) / (S - 1)
  e_obs[!as.vector(loglik$obs_mask), ] <- 0

  one <- ncol(theta_m) == 1
  structure(
    list(estimate = if (one) unname(est) else est,
         influence_site = if (one) c_site[, 1] else c_site,
         influence_obs = if (one) e_obs[, 1] else e_obs,
         obs_mask = loglik$obs_mask,
         obs = obs),
    class = "ijk_nmix"
  )
}


#' Two-Stage Bootstrap Weights for Site/Visit Data
#'
#' Draws bootstrap weights for hierarchical data with sites
#' \eqn{i = 1, \ldots, I} and visits \eqn{j = 1, \ldots, J}.
#' Sites are resampled with replacement,
#' \eqn{u \sim \mathrm{Mult}(I; 1/I, \ldots, 1/I)}; then, within each
#' selected copy of site \eqn{i}, its \eqn{J_i} observed visits are
#' resampled with replacement, so that the total weight of visit \eqn{j} is
#' \eqn{W_{ij} \sim \mathrm{Mult}(u_i J_i; 1/J_i, \ldots, 1/J_i)} and the
#' visit weight is \eqn{v_{ij} = W_{ij} / u_i}. Visits of sites with
#' \eqn{u_i = 0} or \eqn{J_i = 1}, and missing visits, get \eqn{v_{ij} = 1}.
#'
#' Supplying the same weights to several calls of [ijk_bootstrap_2stage()]
#' gives bootstrap replicates that are coherent across quantities.
#'
#' @param B Positive integer giving the number of bootstrap replicates.
#' @param obs_mask Non-empty logical \eqn{I \times J} matrix of observed
#'   visits without missing values, e.g. `loglik$obs_mask` from [nmix_loglik()]
#'   or `!is.na(y)`.
#'
#' @return A list with `site` (\eqn{B \times I}) and `visit`
#'   (\eqn{B \times IJ}, columns ordered as \eqn{i + (j-1)I}).
#'
#' @examples
#' set.seed(1)
#' w <- nmix_weights(5, matrix(TRUE, 4, 3))
#' rowSums(w$site)
#'
#' @export
nmix_weights <- function(B, obs_mask) {
  if (!is.numeric(B) || !is.null(dim(B)) || length(B) != 1L ||
      !is.finite(B) || B < 1 || B != floor(B))
    stop("'B' must be a positive integer.", call. = FALSE)
  if (!is.matrix(obs_mask) || !is.logical(obs_mask) ||
      nrow(obs_mask) < 1L || ncol(obs_mask) < 1L || anyNA(obs_mask))
    stop("'obs_mask' must be a non-empty logical matrix without missing values.",
         call. = FALSE)
  I <- nrow(obs_mask)
  U <- t(stats::rmultinom(B, size = I, prob = rep(1 / I, I)))
  list(site = U, visit = visit_weights(U, obs_mask))
}


#' Visit Weights Given Site Weights
#'
#' Internal helper of [nmix_weights()].
#'
#' @param U Non-empty matrix (\eqn{B \times I}) of non-negative integer
#'   site weights.
#' @param obs_mask Non-empty logical \eqn{I \times J} matrix of observed
#'   visits without missing values; its row count must match `ncol(U)`.
#'
#' @return Numeric matrix (\eqn{B \times IJ}) of visit weights.
#'
#' @keywords internal
visit_weights <- function(U, obs_mask) {
  if (!is.matrix(U) || !is.numeric(U) || nrow(U) < 1L || ncol(U) < 1L ||
      any(!is.finite(U)) || any(U < 0) || any(U != floor(U)))
    stop("'U' must be a non-empty matrix of non-negative integer weights.",
         call. = FALSE)
  if (!is.matrix(obs_mask) || !is.logical(obs_mask) ||
      any(dim(obs_mask) < 1L) || anyNA(obs_mask))
    stop("'obs_mask' must be a non-empty logical matrix without missing values.",
         call. = FALSE)
  B <- nrow(U)
  I <- nrow(obs_mask)
  J <- ncol(obs_mask)
  if (ncol(U) != I)
    stop("ncol(U) must equal nrow(obs_mask).", call. = FALSE)
  V <- matrix(1, B, I * J)
  J_obs <- rowSums(obs_mask)
  for (jj in setdiff(sort(unique(J_obs)), c(0, 1))) {
    sites <- which(J_obs == jj)
    Ui <- U[, sites, drop = FALSE]
    remaining <- Ui * jj
    vis <- matrix(t(vapply(sites, function(i) which(obs_mask[i, ]), integer(jj))),
                  ncol = jj)
    for (k in seq_len(jj)) {
      Wk <- if (k < jj) {
        matrix(stats::rbinom(length(remaining), size = remaining,
                             prob = 1 / (jj - k + 1)), nrow = B)
      } else {
        remaining
      }
      remaining <- remaining - Wk
      cols <- sites + (vis[, k] - 1L) * I
      V[, cols] <- ifelse(Ui > 0, Wk / pmax(Ui, 1), 1)
    }
  }
  V
}


#' Two-Stage Bootstrap Approximation for an N-Mixture Model
#'
#' Generates bootstrap replicates of a posterior mean from a single MCMC
#' run, using the first-order (infinitesimal jackknife) approximation with
#' site-level and visit-level weights:
#' \deqn{\theta^{(b)} \approx E[\theta] + \sum_i (u_{bi} - 1) c_i +
#'       \sum_{i,j} (v_{bij} - 1) e_{ij}.}
#'
#' \describe{
#'   \item{`"twostage"`}{Sites and, within them, visits are resampled
#'     (see [nmix_weights()]).}
#'   \item{`"site"`}{Only sites are resampled (\eqn{v = 1}); equivalent to
#'     [ijk_bootstrap()] with the site-level influence.}
#'   \item{`"visit"`}{Sites are kept (\eqn{u = 1}) and only visits are
#'     resampled within each site.}
#' }
#'
#' For a quantity that involves the latent abundance itself, such as the
#' realised total \eqn{\sum_i N_i}, the likelihood-based replicates only
#' capture the variation through the parameters. Pass `latent_resid`
#' (see [nmix_latent_resid()]) to add the prediction error of the latent
#' variables.
#'
#' @param object An object of class `"ijk_nmix"`.
#' @param B Positive integer giving the number of bootstrap replicates when
#'   `weights` is `NULL`; ignored when `weights` is supplied.
#' @param type One of `"twostage"`, `"site"`, `"visit"`.
#' @param weights Optional weights from [nmix_weights()], to share the same
#'   resamples across quantities.
#' @param latent_resid Optional finite numeric vector (length \eqn{S}) or
#'   finite numeric matrix with \eqn{S} rows and \eqn{K} columns, containing
#'   latent prediction errors; one randomly chosen draw is added to each
#'   replicate.
#'
#' @return An object of class `c("ijk_boot_nmix", "ijk_boot")` with
#'   `estimate`, `replicates` (vector of length \eqn{B}, or
#'   \eqn{B \times K} matrix), `B` and `type`.
#'
#' @seealso [nmix_loglik()], [ijk_nmix()], [nmix_weights()]
#'
#' @examples
#' set.seed(1)
#' I <- 10; J <- 3; S <- 200
#' y <- matrix(rpois(I * J, 3), I, J)
#' lambda <- matrix(rgamma(S * I, 20, 4), S, I)
#' p <- matrix(runif(S * I * J, 0.4, 0.7), S, I * J)
#' ll <- nmix_loglik(y, lambda, p)
#' fit <- ijk_nmix(rowSums(lambda), ll)
#' boot <- ijk_bootstrap_2stage(fit, B = 500)
#' summary(boot)
#'
#' @export
ijk_bootstrap_2stage <- function(object, B = 1000,
                                 type = c("twostage", "site", "visit"),
                                 weights = NULL, latent_resid = NULL) {
  if (!inherits(object, "ijk_nmix"))
    stop("'object' must be created by ijk_nmix().", call. = FALSE)
  type <- match.arg(type)
  obs_mask <- object$obs_mask
  if (!is.matrix(obs_mask) || !is.logical(obs_mask) ||
      nrow(obs_mask) < 1L || ncol(obs_mask) < 1L || anyNA(obs_mask))
    stop("'object' has an invalid observed-visit mask.", call. = FALSE)
  I <- nrow(obs_mask)
  if (is.null(weights)) {
    weights <- nmix_weights(B, obs_mask)
  } else {
    if (!is.list(weights) || !is.matrix(weights$site) ||
        !is.matrix(weights$visit) || !is.numeric(weights$site) ||
        !is.numeric(weights$visit))
      stop("'weights' must be a valid list from nmix_weights().",
           call. = FALSE)
    if (nrow(weights$site) < 1L || nrow(weights$visit) != nrow(weights$site) ||
        ncol(weights$site) != I || ncol(weights$visit) != length(obs_mask))
      stop("'weights' do not match the dimensions of the data.",
           call. = FALSE)
    if (any(!is.finite(weights$site)) || any(!is.finite(weights$visit)) ||
        any(weights$site < 0) || any(weights$site != floor(weights$site)) ||
        any(weights$visit < 0) || any(!is.finite(weights$visit)))
      stop("'weights' must contain finite non-negative values.",
           call. = FALSE)
    if (any(rowSums(weights$site) != I))
      stop("Each row of 'weights$site' must sum to the number of sites.",
           call. = FALSE)
    B <- nrow(weights$site)
  }

  cs <- as.matrix(object$influence_site)
  eo <- as.matrix(object$influence_obs)
  if (!is.numeric(cs) || !is.numeric(eo) ||
      !is.numeric(object$estimate) || !is.null(dim(object$estimate)) ||
      any(!is.finite(cs)) || any(!is.finite(eo)) ||
      nrow(cs) != I || ncol(cs) < 1L || nrow(eo) != length(obs_mask) ||
      ncol(eo) != ncol(cs) || length(object$estimate) != ncol(cs) ||
      any(!is.finite(object$estimate)))
    stop("'object' contains inconsistent or non-finite influence values.",
         call. = FALSE)
  d <- switch(
    type,
    site = (weights$site - 1) %*% cs,
    visit = (visit_weights(matrix(1L, B, I), obs_mask) - 1) %*% eo,
    twostage = (weights$site - 1) %*% cs + (weights$visit - 1) %*% eo
  )
  reps <- sweep(d, 2, object$estimate, "+")

  if (!is.null(latent_resid)) {
    lr <- as.matrix(latent_resid)
    if (!is.numeric(lr) || nrow(lr) < 1L || ncol(lr) != ncol(reps) ||
        any(!is.finite(lr)))
      stop("'latent_resid' must be a finite numeric vector or matrix with one column per quantity.",
           call. = FALSE)
    reps <- reps + lr[sample.int(nrow(lr), B, replace = TRUE), , drop = FALSE]
  }
  colnames(reps) <- names(object$estimate)
  if (ncol(reps) == 1) reps <- reps[, 1]

  structure(
    list(estimate = object$estimate, replicates = reps, B = B, type = type),
    class = c("ijk_boot_nmix", "ijk_boot")
  )
}


#' Latent Prediction Error of a Total Abundance
#'
#' For the realised total \eqn{\sum_{i \in A} N_i} over a set of sites
#' \eqn{A}, returns \eqn{\sum_{i \in A} N_i^{(s)} -
#' \sum_{i \in A} E[N_i \mid y_i, \theta^{(s)}]} for each draw \eqn{s}.
#' Passing it to [ijk_bootstrap_2stage()] adds the uncertainty of the
#' latent abundance that is not captured by the likelihood weights.
#'
#' @param N Finite non-negative integer matrix (\eqn{S \times I}) of
#'   posterior draws of \eqn{N_i}, with dimensions matching `loglik$condN`.
#' @param loglik An object of class `"nmix_loglik"`.
#' @param sites Optional unique integer site indices in `1:ncol(N)`; defaults
#'   to all sites.
#'
#' @return Numeric vector of length \eqn{S}.
#'
#' @export
nmix_latent_resid <- function(N, loglik, sites = NULL) {
  if (!inherits(loglik, "nmix_loglik"))
    stop("'loglik' must be created by nmix_loglik().", call. = FALSE)
  N <- as.matrix(N)
  if (!is.numeric(N) || !all(dim(N) == dim(loglik$condN)) ||
      any(!is.finite(N)) || any(N < 0) || any(N != floor(N)))
    stop("'N' must be a finite non-negative integer matrix matching loglik$condN.",
         call. = FALSE)
  if (is.null(sites)) sites <- seq_len(ncol(N))
  if (!is.numeric(sites) || !is.null(dim(sites)) || length(sites) < 1L ||
      any(!is.finite(sites)) || any(sites != floor(sites)) ||
      any(sites < 1 | sites > ncol(N)) || anyDuplicated(sites))
    stop("'sites' must be unique valid site indices.", call. = FALSE)
  if (any(!is.finite(loglik$condN[, sites, drop = FALSE])))
    stop("'loglik$condN' must be finite for the selected sites.",
         call. = FALSE)
  rowSums(N[, sites, drop = FALSE]) -
    rowSums(loglik$condN[, sites, drop = FALSE])
}


#' Summarize a Two-Stage IJK Bootstrap Object
#'
#' @param object,x An object of class `"ijk_boot_nmix"`.
#' @param level Finite confidence level strictly between 0 and 1 for the
#'   percentile interval.
#' @param digits Number of digits to display.
#' @param ... Currently unused.
#'
#' @return `summary()` returns a data frame with one row per quantity:
#'   `estimate`, `bias`, `se` and the percentile interval `lower`, `upper`.
#'   `print()` returns `x` invisibly.
#'
#' @examples
#' set.seed(1)
#' I <- 10; J <- 3; S <- 200
#' y <- matrix(rpois(I * J, 3), I, J)
#' lambda <- matrix(rgamma(S * I, 20, 4), S, I)
#' p <- matrix(runif(S * I * J, 0.4, 0.7), S, I * J)
#' fit <- ijk_nmix(rowSums(lambda), nmix_loglik(y, lambda, p))
#' boot <- ijk_bootstrap_2stage(fit, B = 200)
#' print(boot)
#'
#' @export
summary.ijk_boot_nmix <- function(object, level = 0.95, ...) {
  if (!is.numeric(level) || !is.null(dim(level)) || length(level) != 1L ||
      !is.finite(level) || level <= 0 || level >= 1)
    stop("'level' must be a finite number strictly between 0 and 1.",
         call. = FALSE)
  reps <- as.matrix(object$replicates)
  if (!is.numeric(reps) || nrow(reps) < 2L || ncol(reps) < 1L ||
      any(!is.finite(reps)))
    stop("'object' must contain at least two finite bootstrap replicates.",
         call. = FALSE)
  a <- (1 - level) / 2
  data.frame(
    estimate = object$estimate,
    bias = colMeans(reps) - object$estimate,
    se = apply(reps, 2, stats::sd),
    lower = apply(reps, 2, stats::quantile, probs = a, names = FALSE),
    upper = apply(reps, 2, stats::quantile, probs = 1 - a, names = FALSE),
    row.names = if (is.null(colnames(reps))) NULL else colnames(reps)
  )
}

#' @rdname summary.ijk_boot_nmix
#' @export
print.ijk_boot_nmix <- function(x, digits = 3, ...) {
  cat("\nIJK Bootstrap Summary (N-mixture, type = \"", x$type, "\")\n\n", sep = "")
  cat("Bootstrap replicates:", x$B, "\n\n")
  print(summary(x), digits = digits)
  invisible(x)
}
