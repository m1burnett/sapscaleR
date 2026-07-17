if (requireNamespace("testthat", quietly = TRUE)) {
  library(testthat)
  library(sapscaleR)
  test_check("sapscaleR")
} else {
  message("testthat is not installed; skipping tests.")
}
