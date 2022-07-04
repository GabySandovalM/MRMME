#' Fitting MRMME via Chan and Mak estimators
#'
#' @param X covariates matrix. It has dimension n x p.
#' @param Y response variables matrix. It has dimension n x q.
#'
fit_ChanMak <- function(X, Y) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  Z <- cbind(X, Y)
  n <- dim(Y)[1]
  q <- dim(Y)[2]
  p <- dim(X)[2]
  r <- p + q
  d <- (p + 1) * (2 + p + 2 * q) / 2

  # S, D and R
  z_bar <- colMeans(Z)
  Z_cen <- t(t(Z) - z_bar)
  S <- 1/n * crossprod(Z_cen)
  D <- diag(eigen(S)$values)
  D_max <- diag(eigen(S)$values[1:p],p,p)
  D_min <- diag(eigen(S)$values[(p+1):r],q,q)
  R <- eigen(S)$vectors
  R_11 <- R[1:p,1:p]
  R_21 <- R[(p+1):r,1:p]

  # mu
  mu_hat <- colMeans(X)

  # phi
  phi_hat <- c("phi" = sum(diag(D_min))/q)

  # B
  B_hat <- R_21 %*% solve(R_11)
  row.names(B_hat) <- colnames(Y)
  colnames(B_hat) <- colnames(X)

  # a
  a_hat <- colMeans(Y) - B_hat %*% mu_hat
  row.names(a_hat) <- colnames(Y)
  colnames(a_hat) <- "intercept"

  # Sigma
  Lambda <- rbind(diag(1, p, p), B_hat)
  C_hat <- crossprod(Lambda)
  C_inv <- solve(C_hat)
  Sigma_hat <- C_inv %*% t(Lambda) %*% S %*% Lambda %*% C_inv - phi_hat*C_inv
  Sigma_hat <- 0.5 * (Sigma_hat + t(Sigma_hat))

  theta <- c(a_hat,B_hat,phi_hat,mu_hat, matrixcalc::vech(Sigma_hat))
  names(theta) <- nms(p,q,TRUE)

  AIC <- 2 * d - 2 * logL(theta, X = X, Y = Y)
  BIC <- d * log(n) - 2 * logL(theta, X = X, Y = Y)
  logLz <- logL(theta, X = X, Y = Y)

  # Results
  res <- list(
    "a" = a_hat,
    "B" = B_hat,
    "phi" = phi_hat,
    "mu_x" = mu_hat,
    "Sigma_x" = Sigma_hat,
    "coef" = cbind(a_hat,B_hat),
    "theta" = theta,
    "AIC" = AIC,
    "BIC" = BIC,
    "logL" = logLz,
    #"d" = d,
    "X" = X,
    "Y" = Y,
    "iter" = "-"
  )
  res
}
