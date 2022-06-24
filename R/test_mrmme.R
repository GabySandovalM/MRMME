#' Title
#'
#' @param X
#' @param Y
#' @param vars.pos
#' @param FIM
#'
#' @return
#'
#' @examples
test_mrmme <- function(X, Y, vars.pos = 1, FIM = "emp", crit = 1e-10) {
  Y <- as.matrix(Y)
  X <- as.matrix(X)
  X1 <- as.matrix(X[, c(vars.pos)])
  X2 <- as.matrix(X[, -c(vars.pos)])
  X <- cbind(X1, X2)

  mod_nr <- fit_EM(X, Y, crit = crit)
  mod_r <- fit_EMtest(X1, X2, Y, crit = crit)

  n <- dim(Y)[1]
  q <- dim(Y)[2]
  p <- dim(X)[2]
  p1 <- dim(X1)[2]
  p2 <- dim(X2)[2]
  r <- p + q
  d <- (p + 1) * (2 + p + 2 * q) / 2
  p_s <- length(matrixcalc::vech(mod_r$Sigma))

  theta_tilde <- mod_r$theta
  theta_hat <- mod_nr$theta
  S_tilde <- matrix(score_matrix_N(theta_tilde, X, Y, total = TRUE), ncol = 1)

  # A matrix and g vector
  k <- p1 * q
  A <- diag(c(rep(0, q),
              rep(1, p1 * q),
              rep(0, p2 * q),
              0,
              rep(0, p),
              rep(0, p_s)))
  colnames(A) <- nam(p, q)
  rownames(A) <- nam(p, q)
  A <- A[colSums(A) != 0, ]
  g <- matrix(rep(0, k), ncol = 1)

  # Information
  FIM_tilde <- fim_N(theta_tilde, X, Y, type = FIM)
  FIM_inv_tilde <- Rfast::spdinv(FIM_tilde)

  FIM_hat <- fim_N(theta_hat, X, Y, type = FIM)
  FIM_inv_hat <- Rfast::spdinv(FIM_hat)

  # LR
  LR <- 2 * (logL(theta_hat, X, Y) - logL(theta_tilde, X, Y))

  # WD
  g_aux <- (A %*% theta_hat) - g
  WD <-  t(g_aux) %*% solve(A %*% FIM_inv_hat %*% t(A)) %*% g_aux

  # SC
  SC <- t(S_tilde) %*% FIM_inv_tilde %*% S_tilde

  # GR
  GR <- t(S_tilde) %*% (theta_hat - theta_tilde)
  GR

  statistics <- round(c(
    "LR" = LR,
    "WD" = WD,
    "SC" = SC,
    "GR" = GR
  ), 5)

  d_nr <- mod_nr$d
  d_r <- mod_nr$d - (p1 * q)
  df <- d_nr - d_r

  p.value <- stats::pchisq(statistics, df, lower.tail = FALSE)

  AIC_nr <- mod_nr$AIC
  AIC_r <- mod_r$AIC

  BIC_nr <- d_nr * log(n) - 2 * mod_nr$logL
  BIC_r <- d_r * log(n) - 2 * mod_r$logL

  logL_nr <- mod_nr$logL
  logL_r <- mod_r$logL

  results <- list(
    "stat" = statistics,
    "p-value" = p.value,
    "df" = df,
    "tilde" = theta_tilde,
    "hat" = theta_hat,
    "d_r" = d_r,
    "d_nr" = d_nr,
    "AIC_r" = AIC_r,
    "AIC_nr" = AIC_nr,
    "logL_r" = logL_r,
    "logL_nr" = logL_nr,
    "BIC_r" = BIC_r,
    "BIC_nr" = BIC_nr
  )
  results
}
