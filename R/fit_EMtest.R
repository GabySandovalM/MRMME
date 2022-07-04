#' Fitting restricted MRMME via EM algorithm
#'
#' @param X1 a matrix of covariates to be tested with dimension n x p1.
#' It is a subset of X.
#' @param X2 a matrix of covariates with dimension n x p2. It is a subset of X
#' that does not include the p1 covariates in `X1`.
#' @param Y response variables matrix. It has dimension n x q.
#' @param crit convergence criterion. Default is 1e-10.
#'
fit_EMtest <- function(X1, X2, Y, crit = 1e-10) {
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  Y <- as.matrix(Y)
  n <- dim(Y)[1]
  q <- dim(Y)[2]
  p1 <- dim(X1)[2]
  p2 <- dim(X2)[2]
  p <- p1 + p2
  r <- p + q
  d <- (p + 1) * (2 + p + 2 * q) / 2
  X <- cbind(X1, X2)
  Z <- cbind(X, Y)

  # Initial Values ----
  mod_lm <- stats::lm(Y ~ X2)
  a <- matrix(stats::coef(mod_lm)[1,],ncol=1) # lm
  B2 <- t(stats::coef(mod_lm)[-1,]) # lm
  phi <- 1
  mu_x <- colMeans(X) # observed
  Sigma_x <- stats::cov(X) # observed

  # a <- rep(0, p)
  # B2 <- matrix(0, nrow = q, ncol = p2)
  # phi <- 1
  # mu_x <- rep(0, p)
  # Sigma_x <- diag(1, p)

  #log-likelihood
  B <- matrix(c(rep(0, (q * p1)), B2), nrow = q, ncol = p)
  theta_est <- c(a, B, phi, mu_x, matrixcalc::vech(Sigma_x))
  logL_k <- logL(theta_est, X, Y)

  # Iterations
  s <- 0
  dif <- 1

  while (dif > crit) {
    s <- s + 1

    # E - step ----
    B2 <- matrix(B[, -c(1:p1)], nrow = q)

    Lambda2 <- rbind(diag(1, p2, p2), B2)
    LambdaR1 <- cbind(diag(1, p1, p1), matrix(0, nrow = p1, ncol = p2))
    LambdaR2 <- cbind(matrix(0, nrow = (q + p2), ncol = p1), Lambda2)
    LambdaR <- rbind(LambdaR1, LambdaR2)

    alpha <- c(rep(0, p), a)

    # psi:  Z_i variance
    psiR <- LambdaR %*% Sigma_x %*% t(LambdaR) + diag(phi, r, r)
    psiR_inv <- Rfast::spdinv(psiR)


    # m y M (x_i|Z_i,theta) ~ N(m_i,M)
    aux <- t(Z) - c(alpha + LambdaR %*% mu_x)
    m <- t(mu_x + tcrossprod(Sigma_x, LambdaR) %*% psiR_inv %*% aux)
    M <-
      Sigma_x - tcrossprod(Sigma_x, LambdaR) %*% psiR_inv %*% LambdaR %*% Sigma_x
    m1 <- m[, 1:p1]
    m2 <- matrix(m[, -c(1:p1)], ncol = p2)
    M11 <- M[1:p1, 1:p1]
    M21 <- M[-c(1:p1), 1:p1]
    M12 <- M[1:p1, -c(1:p1)]
    M22 <- M[-c(1:p1), -c(1:p1)]

    # M - step ----
    y_bar <- colMeans(Y)
    m_bar <- colMeans(m)
    m1_bar <- m_bar[1:p1]
    m2_bar <- m_bar[-c(1:p1)]

    # S_ym
    S_ym2 <- (1 / n) * tcrossprod(t(Y) - y_bar, t(m2) - m2_bar)
    # S_mm
    S_mm <- (1 / n) * tcrossprod(t(m) - m_bar)
    S_m2m2 <- (1 / n) * tcrossprod(t(m2) - m2_bar)

    # theta_hat:
    # Sigma
    Sigma_x <- S_mm + M
    Sigma_x <- 0.5 * (Sigma_x + t(Sigma_x))
    # mu
    mu_x <- m_bar
    # B2
    B2 <- S_ym2 %*% Rfast::spdinv(S_m2m2 + M22)
    # a
    a <- y_bar - B2 %*% m2_bar
    row.names(a) <- colnames(Y)
    colnames(a) <- "intercept"
    # phi_hat
    #  S(theta_hat)
    #  use ||A||^2 = tr(A^t A)
    S_theta_v <- NULL
    for (i in 1:n) {
      S_theta_v[i] <- sum(diag(crossprod(X[i, ] - m[i, ]))) +
        sum(diag(crossprod(Y[i, ] - a - B2 %*% m2[i, ]))) +
        sum(diag(M)) + sum(diag(t(B2) %*% B2 %*% M22))
    }
    S_theta <- sum(S_theta_v)
    phi <- c("phi" = S_theta / (n * r))

    B <- matrix(c(rep(0, (q * p1)), B2), nrow = q, ncol = p)

    #log-likelihood
    theta_est <- c(a, B, phi, mu_x, matrixcalc::vech(Sigma_x))
    logL_k1 <- logL(theta_est, X, Y)
    dif <- abs((logL_k1 / logL_k) - 1)
    logL_k <- logL_k1
  }
  names(theta_est) <- nms(p, q)
  AIC <- 2 * d - 2 * logL(theta_est, X = X, Y = Y)
  BIC <- d * log(n) - 2 * logL(theta_est, X = X, Y = Y)


  res <- list(
    "a" = a,
    "B" = B,
    "phi" = phi,
    "mu_x" = mu_x,
    "Sigma_x" = Sigma_x,
    "coef" = cbind(a, B),
    "theta" = theta_est,
    "d" = d,
    "AIC" = AIC,
    "BIC" = BIC,
    "logL" = logL_k
  )
  res
}

