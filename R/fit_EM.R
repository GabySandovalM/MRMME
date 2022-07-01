#' Fitting MRMME via EM algorithm.
#'
#' @param X covariates matrix. It has dimension n x p.
#' @param Y response variables matrix. It has dimension n x q.
#' @param crit convergence criterion. Default is 1e-10.
#'
#'
fit_EM <- function(X, Y, crit = 1e-10) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  Z <- cbind(X, Y)
  n <- dim(Y)[1]
  q <- dim(Y)[2]
  p <- dim(X)[2]
  r <- p + q
  d <- (p + 1) * (2 + p + 2 * q) / 2


# Initial values ----------------------------------------------------------
  mod_lm <- stats::lm(Y ~ X)
  a <- matrix(stats::coef(mod_lm)[1,],ncol=1) # lm
  B <- t(stats::coef(mod_lm)[-1,]) # lm
  phi <- 1
  mu_x <- colMeans(X) # observed
  Sigma_x <- stats::cov(X) # observed

  # a <- rep(0,p)
  # B <- matrix(0, nrow = q, ncol = p)
  # phi <- 1
  # mu_x <- rep(0,p)
  # Sigma_x <- diag(1, p)

  #log-likelihood
  theta_est <- c(a, B, phi, mu_x, matrixcalc::vech(Sigma_x))
  logL_k <- logL(theta_est, X, Y)

  # Iterations
  s <- 0
  dif <- 1
  logs <- NULL

  while (dif > crit) {
    s <- s + 1
    Lambda <- rbind(diag(1, p, p), B)
    alpha <- c(rep(0, p), a)
# E-step ------------------------------------------------------------------
    # psi:  Z_i variance
    psi <- Lambda %*% Sigma_x %*% t(Lambda) + diag(phi, r, r)
    psi_inv <- Rfast::spdinv(psi)

    # m y M (x_i|Z_i,theta) ~ N(m_i,M)
    aux <- t(Z) - c(alpha + Lambda %*% mu_x)
    m <- t(mu_x + tcrossprod(Sigma_x, Lambda) %*% psi_inv %*% aux)
    M <-
      Sigma_x - tcrossprod(Sigma_x, Lambda) %*% psi_inv %*% Lambda %*% Sigma_x

# M-step ------------------------------------------------------------------

    y_bar <- colMeans(Y)
    m_bar <- colMeans(m)
    # S_ym
    S_ym <- 1 / n * tcrossprod(t(Y) - y_bar, t(m) - m_bar)
    # S_mm
    S_mm <- 1 / n * tcrossprod(t(m) - m_bar)

    # theta_hat:
    # Sigma
    Sigma_x <- S_mm + M
    Sigma_x <- 0.5 * (Sigma_x + t(Sigma_x))
    # mu
    mu_x <- m_bar
    # B
    B <- S_ym %*% Rfast::spdinv(Sigma_x)
    # a
    a <- y_bar - B %*% m_bar
    row.names(a) <- colnames(Y)
    colnames(a) <- "intercept"
    # phi_hat
    #  S(theta_hat)
    #  use ||A||^2 = tr(A^t A)
    S_theta_v <- NULL
    for (i in 1:n) {
      S_theta_v[i] <- sum(diag(crossprod(X[i,] - m[i,]))) +
        sum(diag(crossprod(Y[i,] - a - B %*% m[i,]))) +
        sum(diag(M + t(B) %*% B %*% M))
    }
    S_theta <- sum(S_theta_v)
    phi <- c("phi" = S_theta / (n * r))

    #log-likelihood
    theta_est <- c(a, B, phi, mu_x, matrixcalc::vech(Sigma_x))
    logL_k1 <- logL(theta_est, X, Y)
    logs[s] <- c(logL_k1)
    dif  <- abs((logL_k1 / logL_k) - 1)
    logL_k <- logL_k1
  }
  names(theta_est) <- nam2(p, q)

# Results -----------------------------------------------------------------

  AIC <- 2 * d - 2 * logL(theta_est, X = X, Y = Y)
  BIC <- d * log(n) - 2 * logL(theta_est, X = X, Y = Y)
  logLz <- logL(theta_est, X = X, Y = Y)

  res <- list(
    "a" = a,
    "B" = B,
    "phi" = phi,
    "mu_x" = mu_x,
    "Sigma_x" = Sigma_x,
    "coef" = cbind(a, B),
    "theta" = theta_est,
    #"d" = d,
    "iter" = s,
    # "logs" = logs,
    # "dif" = dif,
    "AIC" = AIC,
    "BIC" = BIC,
    "logL" = logLz,
    "X" = X,
    "Y" = Y
  )
  return(res)
}
