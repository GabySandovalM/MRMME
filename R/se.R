#' Compute standard errors for the parameter estimates
#'
#' @param theta a vector with the parameter estimates.
#' @param X covariates matrix. It has dimension n x p.
#' @param Y response variables matrix. It has dimension n x q.
#' @param type the type of estimator of the Fisher information matrix.
#'
#' @return a list with the covariance matrix of the parameters estimates and
#' a vector with standard errors.
#'
se <- function(theta, X, Y, type = "emp"){
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  Z <- cbind(X, Y)
  n <- dim(Y)[1]
  q <- dim(Y)[2]
  p <- dim(X)[2]
  r <- p + q
  d <- (p + 1) * (2 + p + 2 * q) / 2
  type <- type

  fim <- fim(theta, X, Y, type)
  D <-  Rfast::spdinv(fim)
  D <- 0.5 * (D + t(D))
  se <- sqrt(diag(D))
  colnames(D) <- nms(p,q)
  rownames(D) <- nms(p,q)
  names(se) <- nms(p,q,TRUE)
  res <- list("covmat" = D,
              "se" = se)
  res
}
