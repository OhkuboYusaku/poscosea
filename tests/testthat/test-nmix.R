## Small synthetic "posterior draws" for an N-mixture model -------------------
sim_nmix_draws <- function(I = 12, J = 3, S = 60, seed = 1) {
  set.seed(seed)
  lam_true <- exp(1 + 0.3 * rnorm(I))
  N_true <- rpois(I, lam_true)
  p_true <- matrix(stats::plogis(rnorm(I * J, 0, 0.5)), I, J)
  y <- matrix(rbinom(I * J, rep(N_true, J), p_true), I, J)
  y[2, 3] <- NA            # one missing visit
  y[5, 2:3] <- NA          # a site with a single observed visit
  lambda <- exp(matrix(rep(log(lam_true), each = S), S, I) +
                  matrix(rnorm(S * I, 0, 0.1), S, I))
  p <- stats::plogis(matrix(rep(stats::qlogis(as.vector(p_true)), each = S), S, I * J) +
                       matrix(rnorm(S * I * J, 0, 0.1), S, I * J))
  ymax <- apply(y, 1, max, na.rm = TRUE)
  N <- matrix(rep(ymax, each = S), S, I) + matrix(rpois(S * I, 1), S, I)
  list(y = y, lambda = lambda, p = p, N = N, I = I, J = J, S = S)
}

## Brute-force reference for one draw and one site
brute_site <- function(y_i, lam, p_i, Nmax_extra = 300) {
  obs <- !is.na(y_i)
  Ns <- seq.int(max(y_i[obs]), max(y_i[obs]) + Nmax_extra)
  lb <- sapply(which(obs), function(j) stats::dbinom(y_i[j], Ns, p_i[j], log = TRUE))
  lt <- stats::dpois(Ns, lam, log = TRUE) + rowSums(matrix(lb, length(Ns)))
  w <- exp(lt - max(lt)); w <- w / sum(w)
  list(site = max(lt) + log(sum(exp(lt - max(lt)))),
       condN = sum(w * Ns),
       obs_rb = colSums(w * matrix(lb, length(Ns))))
}

d <- sim_nmix_draws()
ll <- nmix_loglik(d$y, d$lambda, d$p, N = d$N)

test_that("nmix_loglik returns objects of the right shape", {
  expect_s3_class(ll, "nmix_loglik")
  expect_equal(dim(ll$site), c(d$S, d$I))
  expect_equal(dim(ll$obs_rb), c(d$S, d$I * d$J))
  expect_equal(dim(ll$obs_joint), c(d$S, d$I * d$J))
  expect_equal(dim(ll$condN), c(d$S, d$I))
  expect_identical(ll$obs_mask, !is.na(d$y))
  expect_null(nmix_loglik(d$y, d$lambda, d$p)$obs_joint)
})

test_that("site and observation log-likelihoods match a brute-force computation", {
  for (s in c(1, 17, d$S)) for (i in c(1, 2, 5, d$I)) {
    p_i <- d$p[s, i + (seq_len(d$J) - 1) * d$I]
    ref <- brute_site(d$y[i, ], d$lambda[s, i], p_i)
    expect_equal(ll$site[s, i], ref$site, tolerance = 1e-10)
    expect_equal(ll$condN[s, i], ref$condN, tolerance = 1e-10)
    k_obs <- i + (which(!is.na(d$y[i, ])) - 1) * d$I
    expect_equal(ll$obs_rb[s, k_obs], ref$obs_rb, tolerance = 1e-10)
  }
})

test_that("missing visits contribute zero and obs_joint uses the N draws", {
  miss <- which(is.na(as.vector(d$y)))
  expect_true(all(ll$obs_rb[, miss] == 0))
  expect_true(all(ll$obs_joint[, miss] == 0))
  k <- 3 + (1 - 1) * d$I
  expect_equal(ll$obs_joint[, k],
               stats::dbinom(d$y[3, 1], d$N[, 3], d$p[, k], log = TRUE))
})

test_that("a 3-d array for p gives the same result", {
  p_arr <- array(d$p, c(d$S, d$I, d$J))
  ll2 <- nmix_loglik(d$y, d$lambda, p_arr)
  expect_equal(ll2$site, ll$site)
  expect_equal(ll2$obs_rb, ll$obs_rb)
})

test_that("nmix_loglik validates its inputs", {
  expect_error(nmix_loglik(d$y, d$lambda[, -1], d$p))
  expect_error(nmix_loglik(d$y, d$lambda, d$p[, -1]))
  y_bad <- d$y; y_bad[1, ] <- NA
  expect_error(nmix_loglik(y_bad, d$lambda, d$p))
  expect_warning(nmix_loglik(d$y, d$lambda * 20, d$p, Nmax_extra = 5))
})

test_that("draws with zero likelihood are flagged and rejected by ijk_nmix", {
  p0 <- d$p
  p0[1, ] <- 0                                   # p = 0 cannot produce y > 0
  expect_warning(ll0 <- nmix_loglik(d$y, d$lambda, p0), "zero likelihood")
  expect_true(any(ll0$site[1, ] == -Inf))
  expect_true(all(is.finite(ll0$site[-1, ])))
  expect_error(ijk_nmix(rowSums(d$lambda), ll0), "non-finite")
})

test_that("nmix_draws extracts columns in the right order", {
  I <- d$I; J <- d$J
  draws <- cbind(d$lambda, d$p, d$N)
  colnames(draws) <- c(paste0("lam[", 1:I, "]"),
                       paste0("p[", rep(1:I, J), ",", rep(1:J, each = I), "]"),
                       paste0("N[", 1:I, "]"))
  draws <- draws[, sample(ncol(draws))]          # shuffle the columns
  ex <- nmix_draws(draws, I, J, lambda = "lam")
  expect_equal(unname(ex$lambda), unname(d$lambda))
  expect_equal(unname(ex$p), unname(d$p))
  expect_equal(unname(ex$N), unname(d$N))
  ex2 <- nmix_draws(list(draws[1:30, ], draws[31:60, ]), I, J, lambda = "lam")
  expect_equal(ex2$p, ex$p)
  expect_error(nmix_draws(draws, I, J))           # "lambda[...]" is absent
})

test_that("two-stage weights resample visits within the selected sites", {
  set.seed(2)
  mask <- ll$obs_mask
  w <- nmix_weights(400, mask)
  I <- nrow(mask); J <- ncol(mask)
  expect_equal(dim(w$site), c(400, I))
  expect_equal(dim(w$visit), c(400, I * J))
  expect_true(all(rowSums(w$site) == I))
  J_obs <- rowSums(mask)
  for (i in seq_len(I)) {
    cols <- i + (which(mask[i, ]) - 1) * I
    sel <- w$site[, i] > 0
    # observed visits of a selected site have weights summing to J_i
    expect_equal(rowSums(w$visit[sel, cols, drop = FALSE]), rep(J_obs[i], sum(sel)))
    # sites that were not selected, or with a single visit, keep weight 1
    expect_true(all(w$visit[!sel, cols] == 1))
    if (J_obs[i] == 1) expect_true(all(w$visit[, cols] == 1))
  }
  miss <- which(!as.vector(mask))
  expect_true(all(w$visit[, miss] == 1))
})

test_that("ijk_nmix computes covariances consistent with ijk()", {
  theta <- rowSums(d$lambda)
  fit <- ijk_nmix(theta, ll)
  expect_s3_class(fit, "ijk_nmix")
  expect_equal(fit$estimate, mean(theta))
  expect_equal(fit$influence_site, compute_if(theta, ll$site))
  expect_equal(mean(fit$influence_site), 0, tolerance = 1e-12)
  expect_length(fit$influence_obs, d$I * d$J)
  k <- 4
  expect_equal(fit$influence_obs[k], stats::cov(theta, ll$obs_rb[, k]))
  expect_true(all(fit$influence_obs[!as.vector(ll$obs_mask)] == 0))

  fit_j <- ijk_nmix(rowSums(d$N), ll, obs = "joint")
  expect_equal(fit_j$influence_obs[k], stats::cov(rowSums(d$N), ll$obs_joint[, k]))
  expect_error(ijk_nmix(theta, nmix_loglik(d$y, d$lambda, d$p), obs = "joint"))
  expect_error(ijk_nmix(theta[-1], ll))
  expect_error(ijk_nmix(theta, list(site = ll$site)))
})

test_that("matrix theta gives one column per quantity", {
  th <- cbind(a = rowSums(d$lambda), b = d$lambda[, 1])
  fit <- ijk_nmix(th, ll)
  expect_equal(dim(fit$influence_site), c(d$I, 2))
  expect_equal(dim(fit$influence_obs), c(d$I * d$J, 2))
  expect_equal(fit$influence_site[, "b"], ijk_nmix(th[, "b"], ll)$influence_site)
  boot <- ijk_bootstrap_2stage(fit, B = 50)
  expect_equal(dim(boot$replicates), c(50, 2))
  expect_equal(rownames(summary(boot)), c("a", "b"))
})

test_that("two-stage bootstrap replicates are built from site and visit parts", {
  theta <- rowSums(d$lambda)
  fit <- ijk_nmix(theta, ll)
  set.seed(3)
  w <- nmix_weights(200, ll$obs_mask)
  b_two <- ijk_bootstrap_2stage(fit, type = "twostage", weights = w)
  b_site <- ijk_bootstrap_2stage(fit, type = "site", weights = w)
  expect_s3_class(b_two, "ijk_boot_nmix")
  expect_s3_class(b_two, "ijk_boot")
  expect_length(b_two$replicates, 200)
  expect_equal(b_site$replicates,
               fit$estimate + drop((w$site - 1) %*% fit$influence_site))
  expect_equal(b_two$replicates - b_site$replicates,
               drop((w$visit - 1) %*% fit$influence_obs))
  b_vis <- ijk_bootstrap_2stage(fit, B = 30, type = "visit")
  expect_length(b_vis$replicates, 30)
})

test_that("site bootstrap variance matches its analytic value", {
  theta <- rowSums(d$lambda)
  fit <- ijk_nmix(theta, ll)
  set.seed(4)
  b <- ijk_bootstrap_2stage(fit, B = 20000, type = "site")
  c_i <- fit$influence_site
  v_exact <- sum(c_i^2) - sum(c_i)^2 / length(c_i)
  expect_equal(stats::var(b$replicates), v_exact, tolerance = 0.05)
  b2 <- ijk_bootstrap_2stage(fit, B = 20000, type = "twostage")
  expect_gt(stats::var(b2$replicates), stats::var(b$replicates))
})

test_that("zero influence gives constant replicates and latent_resid is added", {
  fit <- ijk_nmix(rep(5, d$S), ll)
  b <- ijk_bootstrap_2stage(fit, B = 40)
  expect_equal(b$replicates, rep(5, 40))
  lr <- nmix_latent_resid(d$N, ll)
  expect_length(lr, d$S)
  expect_equal(lr, rowSums(d$N) - rowSums(ll$condN))
  set.seed(5)
  b_lat <- ijk_bootstrap_2stage(fit, B = 40, latent_resid = lr)
  # each replicate is the estimate plus one of the latent residual draws
  expect_true(all(vapply(b_lat$replicates - 5, function(r) any(abs(r - lr) < 1e-10), logical(1))))
  expect_error(ijk_bootstrap_2stage(fit, B = 10, latent_resid = cbind(lr, lr)))
})

test_that("print and summary methods work", {
  fit <- ijk_nmix(rowSums(d$lambda), ll)
  b <- ijk_bootstrap_2stage(fit, B = 100)
  s <- summary(b)
  expect_s3_class(s, "data.frame")
  expect_named(s, c("estimate", "bias", "se", "lower", "upper"))
  expect_output(print(b), "type = \"twostage\"")
  expect_invisible(print(b))
})

test_that("site influence can be used for leave-one-site-out", {
  fit <- ijk_nmix(rowSums(d$lambda), ll)
  loo <- ijk_lko(fit$estimate, fit$influence_site, k = 1)
  expect_length(loo$delta, d$I)
})
