test_that("ijk object is created", {
  theta <- rnorm(100)
  loglik <- matrix(
    rnorm(100 * 10),
    100,
    10
  )
  fit <- ijk(theta, loglik)
  expect_s3_class(fit, "ijk")
})
