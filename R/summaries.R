#' Calculate the 0.05 quantile
#'
#' @param ... Arguments passed to the `quantile()` function.
#' @return Sample quantile.
#' @export
q05 <- function(...) {
  quantile(..., probs = 0.05)
}

#' Calculate the 0.95 quantile
#'
#' @param ... Arguments passed to the `quantile()` function.
#' @return Sample quantile.
#' @export
q95 <- function(...) {
  quantile(..., probs = 0.95)
}

#' Calculate the 0.25 quantile
#'
#' @param ... Arguments passed to the `quantile()` function.
#' @return Sample quantile.
#' @export
q25 <- function(...) {
  quantile(..., probs = 0.25)
}

#' Calculate the 0.50 quantile
#'
#' @param ... Arguments passed to the `quantile()` function.
#' @return Sample quantile.
#' @export
q50 <- function(...) {
  quantile(..., probs = 0.50)
}

#' Calculate the 0.75 quantile
#'
#' @param ... Arguments passed to the `quantile()` function.
#' @return Sample quantile.
#' @export
q75 <- function(...) {
  quantile(..., probs = 0.75)
}
