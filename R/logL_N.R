#' Title
#'
#' @param theta
#' @param X
#' @param Y
#'
#' @return
#'
#' @examples
logL <- function(theta, X, Y) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  Z <- cbind(X,Y)
  n <- dim(Y)[1]
  q <- dim(Y)[2]
  p <- dim(X)[2]
  r <- p + q
  d <- (p + 1) * (2 + p + 2 * q) / 2

  # theta --> (a,B,phi,mu,Sigma)
  theta <- as.vector(theta)
  a <- matrix(theta[1:q], ncol = 1, nrow = q)
  B <- matrix(theta[(q + 1):(q + p * q)], nrow = q, ncol = p)
  phi <- theta[q + p * q + 1]
  mu_x <- theta[(q + p * q + 2):(q + p * q + 1 + p)]
  Sigma_x <- ks::invvech(theta[(q + p * q + 2 + p):d])
  Sigma_x <- 0.5 * (Sigma_x + t(Sigma_x))

  # eta and psi
  Lambda <- rbind(diag(1, p, p), B)
  alpha <- c(rep(0,p),a)
  eta <- c(alpha + Lambda %*% mu_x)
  psi <- Lambda %*% tcrossprod(Sigma_x, Lambda) + diag(phi, r, r)
  psi <- 0.5 * (psi + t(psi))

  # logL MRMME
  cte <- - r * 0.5 * log(2 * pi)
  det2 <- det(phi*diag(1,r,r)) * det(diag(1,p,p) + (1/phi) * crossprod(Lambda) %*% Sigma_x)
  log_det <- log(det2)
  delta_i <- stats::mahalanobis(Z,eta,psi)
  logL_i <- cte - log_det/2 - 1/2 * delta_i
  logL_MRMME <- sum(logL_i)
  logL_MRMME
}
