#' function to generate the names of parameters

nam <- function(p,q){
  name_theta <- NULL
  for (i in 1:q) {
    name_a <- paste0("a_", i)
    name_theta <- c(name_theta, name_a)
  }

  for (k in 1:p) {
    for (l in 1:q) {
      name_b <- paste0("b_", l, k)
      name_theta <- c(name_theta, name_b)
    }
  }

  name_theta <- c(name_theta, "phi")

  for (i in 1:p) {
    name_mu <- paste0("mu_", i)
    name_theta <- c(name_theta, name_mu)
  }

  for (k in 1:p) {
    for (l in 1:p) {
      if (l >= k) {
        name_s <- paste0("sigma_", l, k)
        name_theta <- c(name_theta, name_s)
      }
    }
  }
  name_theta
}


nam2 <- function(p,q){
  name_theta <- NULL
  for (i in 1:q) {
    name_a <- paste0("a[", i,"]")
    name_theta <- c(name_theta, name_a)
  }

  for (k in 1:p) {
    for (l in 1:q) {
      name_b <- paste0("b[", l, k,"]")
      name_theta <- c(name_theta, name_b)
    }
  }

  name_theta <- c(name_theta, "phi")

  for (i in 1:p) {
    name_mu <- paste0("mu[", i, "]")
    name_theta <- c(name_theta, name_mu)
  }

  for (k in 1:p) {
    for (l in 1:p) {
      if (l >= k) {
        name_s <- paste0("sigma[", l, k, "]")
        name_theta <- c(name_theta, name_s)
      }
    }
  }
  name_theta
}

#' Returns the Wilson-Hilferty transformation of Multivariate Normal random variables
#' Z matrix of data following a Multivariate Normal distribution (with p columns).
#' m the mean vector (with dimension p x 1).
#' S the covariance matrix (with dimension p x p).

N_wilson_hilferty <- function(Z,m,S){
  d <- stats::mahalanobis(Z, m, S)
  p <- nrow(S)
  d <- d / p
  z_w.h <- (d ^ (1 / 3) - (1 - 2 / (9 * p))) / sqrt(2 / (9 * p))
  z_w.h
}


#' generate the envelopes in the qqplot of wh-transformation, it generates
#' envelopes to X, Y and Z.

envlp <- function(X, Y, theta, nsim = 50, conf = 0.95){
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  n <- dim(Y)[1]
  q <- dim(Y)[2]
  p <- dim(X)[2]
  r <- p + q
  d <- (p + 1) * (2 + p + 2 * q) / 2

  # theta --> (a,B,phi,mu,Sigma)
  theta <- c(theta)
  a <- matrix(theta[1:q], ncol = 1, nrow = q)
  B <- matrix(theta[(q + 1):(q + p * q)], nrow = q, ncol = p)
  phi <- theta[q + p * q + 1]
  mu_x <- theta[(q + p * q + 2):(q + p * q + 1 + p)]
  Sigma_x <- ks::invvech(theta[(q + p * q + 2 + p):d])
  Sigma_x <- 0.5 * (Sigma_x + t(Sigma_x))

  # sims
  simX <- matrix(0, nrow = n, ncol = nsim)
  simY <- matrix(0, nrow = n, ncol = nsim)
  simZ <- matrix(0, nrow = n, ncol = nsim)
  for (i in 1:nsim) {
    data_sim <- sim_Ndata(
      n = n,
      p = p,
      q = q,
      a = a,
      B = B,
      phi = phi,
      mu_x = mu_x,
      Sigma_x = Sigma_x)

    fit <- fit_ChanMak(X = data_sim$X, Y = data_sim$Y)

    # eta, psi, mu_Y, Sigma_Y
    Lambda <- rbind(diag(1, p, p), fit$B)
    alpha <- c(rep(0, p), fit$a)
    eta_fit <- c(alpha + Lambda %*% fit$mu_x)
    psi_fit <-
      Lambda %*% tcrossprod(fit$Sigma_x, Lambda) + diag(fit$phi, r, r)
    psi_fit <- 0.5 * (psi_fit + t(psi_fit))
    mu_Y_fit <- fit$a + fit$B %*% matrix(fit$mu_x, ncol = 1)
    Sigma_Y_fit <-
      fit$B %*% fit$Sigma_x %*% t(fit$B) + fit$phi * diag(1, q, q)
    Sigma_Y_fit <- 0.5 * (Sigma_Y_fit + t(Sigma_Y_fit))

    Xwh <- N_wilson_hilferty(data_sim$X, fit$mu_x, fit$Sigma_x)
    Ywh <- N_wilson_hilferty(data_sim$Y, mu_Y_fit, Sigma_Y_fit)
    Zwh <- N_wilson_hilferty(cbind(data_sim$X, data_sim$Y), eta_fit, psi_fit)

    simX [, i] <- sort(Xwh)
    simY [, i] <- sort(Ywh)
    simZ [, i] <- sort(Zwh)
  }

  #bandas para X
  limsX <- matrix(, nrow = n, ncol = 2)
  for (i in 1:n) {
    limsX[i,] <- stats::quantile(simX[i,], probs = c((1 - conf) / 2, 1 - (1 - conf) / 2))
  }
  #bandas para Y
  limsY <- matrix(, nrow = n, ncol = 2)
  for (i in 1:n) {
    limsY[i,] <- stats::quantile(simY[i,], probs = c((1 - conf) / 2, 1 - (1 - conf) / 2))
  }
  #bandas para Z
  limsZ <- matrix(, nrow = n, ncol = 2)
  for (i in 1:n) {
    limsZ[i,] <- stats::quantile(simZ[i,], probs = c((1 - conf) / 2, 1 - (1 - conf) / 2))
  }
  res <- list("bandX" = limsX,
              "bandY" = limsY,
              "bandZ" = limsZ)
  return(res)
}

# GENERATE SIGNIFICANCE INDEX

defast = function(x){
  if(x>0.1){ast = " "}else
  {
    if(x>0.05){ast = "."}else
    {
      if(x>0.01){ast = "*"}else
      {
        if(x>0.001){ast = "**"}else
        {
          {
            ast = "***"
          }
        }
      }
    }
  }
  return(ast)
}
