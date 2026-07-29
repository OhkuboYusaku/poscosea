library(testthat)

test_that("ijk object is created correctly", {
  set.seed(123)
  theta <- rnorm(100)
  loglik <- matrix(
    rnorm(100 * 10),
    nrow = 100,
    ncol = 10
  )
  fit <- ijk(theta, loglik)
  expect_s3_class(fit, "ijk")
  expect_true(is.numeric(fit$estimate))
  expect_length(fit$estimate, 1)
  expect_true(is.numeric(fit$influence))
  expect_length(fit$influence, 10)
})

test_that("compute_if returns one value per observation", {
  set.seed(123)
  theta <- rnorm(100)
  loglik <- matrix(
    rnorm(100 * 15),
    nrow = 100,
    ncol = 15
  )
  infl <- compute_if(theta, loglik)
  expect_type(infl, "double")
  expect_length(infl, 15)
})


test_that("compute_if is centered", {
  set.seed(123)
  theta <- rnorm(100)
  loglik <- matrix(
    rnorm(100 * 20),
    nrow = 100,
    ncol = 20
  )
  infl <- compute_if(theta, loglik)
  expect_equal(
    mean(infl),
    0,
    tolerance = 1e-12
  )
})

test_that("input dimension mismatch is detected", {
  theta <- rnorm(100)
  loglik <- matrix(
    rnorm(50 * 10),
    nrow = 50,
    ncol = 10
  )
  expect_error(
    ijk(theta, loglik)
  )
})

test_that("non-numeric theta is rejected", {
  expect_error(
    ijk(
      letters[1:10],
      matrix(rnorm(100), 10, 10)
    )
  )
})

test_that("bootstrap object is created correctly", {

  set.seed(123)

  out <- ijk_bootstrap(
    estimate = 0,
    influence = rnorm(20),
    B = 50
  )
  expect_s3_class(out, "ijk_boot")
  expect_true(is.numeric(out$estimate))
  expect_equal(out$estimate, 0)
  expect_true(is.numeric(out$replicates))
  expect_length(out$replicates, 50)
  expect_equal(out$B, 50)

})

test_that("bootstrap with zero influence is constant", {
  out <- ijk_bootstrap(
    estimate = 5,
    influence = rep(0, 20),
    B = 100
  )
  expect_equal(
    out$replicates,
    rep(5, 100)
  )

})

test_that("leave-k-out returns all combinations", {

  out <- ijk_lko(
    estimate = 0,
    influence = rnorm(20),
    k = 2
  )

  expect_s3_class(
    out,
    "ijk_lko"
  )

  expect_length(
    out$lko_estimates,
    choose(20, 2)
  )

  expect_length(
    out$omitted,
    choose(20, 2)
  )

})

test_that("leave-k-out with zero influence is constant", {

  out <- ijk_lko(
    estimate = 5,
    influence = rep(0, 20),
    k = 3
  )

  expect_equal(
    out$lko_estimates,
    rep(
      5,
      choose(20, 3)
    )
  )

  expect_equal(
    out$delta,
    rep(
      0,
      choose(20, 3)
    )
  )
})

test_that("invalid k is rejected", {

  expect_error(
    ijk_lko(
      estimate = 0,
      influence = rnorm(10),
      k = 11
    )
  )

})

test_that("leave-k-out stores metadata correctly", {

  out <- ijk_lko(
    estimate = 3,
    influence = rnorm(10),
    k = 2
  )

  expect_equal(
    out$estimate,
    3
  )

  expect_equal(
    out$k,
    2
  )

})

test_that("plot method runs", {

  out <- ijk_lko(
    estimate = 0,
    influence = rnorm(20),
    k = 1
  )

  expect_invisible(
    plot(out)
  )

})
