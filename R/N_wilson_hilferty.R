#' Wilson Hilferty transformation
#'
#' Returns the Wilson-Hilferty transformation of Multivariate Normal random variables.
#'
#' @param Z matrix of data following a Multivariate Normal distribution (with p columns).
#' @param m the mean vector (with dimension p x 1).
#' @param S the covariance matrix (with dimension p x p).
#'
#' @return a vector with the transformed distance.
#'
N_wilson_hilferty <- function(Z,m,S){
  d <- stats::mahalanobis(Z, m, S)
  p <- nrow(S)
  d <- d / p
  z_w.h <- (d ^ (1 / 3) - (1 - 2 / (9 * p))) / sqrt(2 / (9 * p))
  z_w.h
}
