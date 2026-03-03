# Helpers for expensive statistical validation tests.

isFullValidationEnabled <- function() {
  flag <- trimws(Sys.getenv("MEDUSA_FULL_VALIDATION", unset = ""))
  tolower(flag) %in% c("1", "true", "yes")
}


skip_if_not_full_validation <- function() {
  if (!isFullValidationEnabled()) {
    testthat::skip("Set MEDUSA_FULL_VALIDATION=true to run slow validation tests.")
  }

  invisible(TRUE)
}
