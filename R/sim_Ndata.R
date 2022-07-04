#' Simulate data from a MRMME
#'
#' @param n number of observations.
#' @param p number of covariates.
#' @param q number of response variables.
#' @param a vector of intercepts.
#' @param B a matrix fo regression coefficients.
#' @param phi a numeric value representing the measurement error variance.
#' @param mu_x mean vector of the unobserved covariates.
#' @param Sigma_x covariance matrix of the unobserved covariates.
#' @param testB if `TRUE` the function return \emph{Y} simulates under \eqn{H_0}.
#' The default value is `FALSE``.`
#' @param vars.pos a numerical value indicating the position (column number in X)
#' of the covariate to be tested under \eqn{H_0}. By default the function the
#' assume the covariate tho be tested is the first.
#'
#' @return a list with simulated data from a MRMME and the parameters used to simulate.
#'
sim_Ndata <- function(n, p, q, a, B, phi, mu_x, Sigma_x, testB = FALSE, vars.pos = c(1)) {
  r <- p + q

  # True x values, from a N(mu_x, Sigma_x)
  x <- matrix(, ncol = p, nrow = n)
  for (i in 1:n) {
    x[i,] <- mvtnorm::rmvnorm(1, mu_x, Sigma_x, method = "chol")
  }

  # Simulate errors
  Sigma_e <- diag(phi, r, r)

  # Errors vector, they are independent and come from a multivariate normal
  e <- matrix(, ncol = r, nrow = n)
  for (i in 1:n) {
    e[i,] <- mvtnorm::rmvnorm(1, rep(0, r), Sigma_e, method = "chol")
  }
  e1 <- as.matrix(e[, 1:p])
  e2 <- as.matrix(e[, -(1:p)])

  # Simulation of Y
  Y <- matrix(, ncol = q, nrow = n)
  for (i in 1:n) {
    Y[i,] <- a + B %*% x[i,] + e2[i,]
  }

  # Simulation of Y under H0
  Y_h0 <- matrix(, ncol = q, nrow = n)
  B_h0 <- B
  B_h0[,vars.pos] <- 0

  for (i in 1:n) {
    Y_h0[i,] <- a + B_h0 %*% x[i,] + e2[i,]
  }

  # X observed
  X <- x + e1

  # Simulation of Z
  alpha <- matrix(c(rep(0, p), a), ncol = 1)
  Lambda <- rbind(diag(1, p), B)
  Z <- matrix(, ncol = r, nrow = n)
  for (i in 1:n) {
    Z[i,] <- alpha + Lambda %*% x[i,] + e[i,]
  }

  # eta
  eta <- alpha + Lambda %*% matrix(mu_x, ncol = 1)

  # psi
  psi <- Lambda %*% Sigma_x %*% t(Lambda) + diag(phi, r, r)

  theta_sim <- c(a, B, phi, mu_x, matrixcalc::vech(Sigma_x))
  names(theta_sim) <- nms(p, q)

  # Resultados
  data = list(
    "X" = X,
    "Y" = Y,
    "Z" = Z,
    "a" = a,
    "B" = B,
    "phi" = phi,
    "mu.x" = mu_x,
    "Sigma.x" = Sigma_x,
    "eta" = eta,
    "psi" = psi,
    "theta" = theta_sim

  )
  if (testB == FALSE) {
    return(data)
  } else{
    data_h0 = list(
      "X" = X,
      "Y" = Y,
      "Y.h0" = Y_h0,
      "Z" = Z,
      "a" = a,
      "B" = B,
      "phi" = phi,
      "mu.x" = mu_x,
      "Sigma.x" = Sigma_x,
      "eta" = eta,
      "psi" = psi,
      "theta" = theta_sim
    )
    return(data_h0)
  }
}

