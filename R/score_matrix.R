#' Computes the score vectors in MRMME.
#'
#' Given the MLE estimates of theta and the observed data, this function
#' computes the individual score vector for each observation.
#'
#' @param theta a vector with the parameter estimates.
#' @param X covariates matrix. It has dimension n x p.
#' @param Y response variables matrix. It has dimension n x q.
#' @param total if `FALSE` (the default) returns the score matrix of dimension
#' n x d, where n is the number of observations and d is the total number of parameters.
#' It means the i\emph{th} row corresponds to the i\emph{th} score vector.
#' If `TRUE` returns a d-dimensional score vector (the sum of the score matrix columns).
#'
#'
#' @return if `total = TRUE` it returns a matrix. If `total = FALSE` it returns a vector that corresponds
#' to the total score. Given the MLE, each element in the total score is approximately equal to zero.
#'
score_matrix <- function(theta, X, Y, total = FALSE) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  Z <- cbind(X, Y)
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

  # eta and psi
  Lambda <- rbind(diag(1, p, p), B)
  alpha <- c(rep(0, p), a)
  eta <- c(alpha + Lambda %*% mu_x)
  psi <- Lambda %*% tcrossprod(Sigma_x, Lambda) + diag(phi, r, r)
  psi <- 0.5 * (psi + t(psi))
  psi_inv <- Rfast::spdinv(psi)

  # Score function
  score <- function(psi_j, partial_delta_j) {
    score =  -0.5 * (matrixcalc::matrix.trace(psi_inv %*% psi_j) + partial_delta_j)
    score
  }

  # jth - partial derivative
  partial_delta_j <- function(psi_j_inv, eta_j) {
    pdj <- stats::mahalanobis(Z, eta, psi_j_inv, inverted = TRUE) -
      2 * (t(t(Z) - eta) %*% psi_inv %*% eta_j)
    pdj
  }

  # score a----
  score_a <- matrix(, nrow = n, ncol = q)
  index <- NULL
  for (i in 1:q) {
    psi_j_a <- matrix(0, nrow = r, ncol = r)
    psi_j_inv_a <- -psi_inv %*% psi_j_a %*% psi_inv
    eta_j_a <- matrix(c(rep(0, p), diag(q)[, i]), ncol = 1)
    partial_delta_j_a <- partial_delta_j(psi_j_inv_a, eta_j_a)
    score_a[, i] <- score(psi_j_a, partial_delta_j_a)
  }

  # score vec(B)----
  # E_mat: matrix with 1 in k,l position
  E_mat <- list()
  for (l in 1:p) {
    E_mat[[l]] <- list()
    for (k in 1:q) {
      E_aux =  matrix(0, nrow = q, ncol = p)
      E_aux[k, l] =  1
      E_mat[[l]][[k]] <- E_aux
    }
  }
  E_mat <- unlist(E_mat, recursive = FALSE)

  eta_j_b <- list()
  for (k in 1:(p * q)) {
    aux <- E_mat[[k]] %*% mu_x
    eta_j_b[[k]] <- rbind(matrix(rep(0, p), ncol = 1), aux)
  }


  psi_j_b <- list()
  for (k in 1:(p * q)) {
    psi_j_b[[k]] <- list()
    row1 = cbind(matrix(0, nrow = p, ncol = p), tcrossprod(Sigma_x, E_mat[[k]]))
    row2 = cbind(E_mat[[k]] %*% Sigma_x,
                 B %*% tcrossprod(Sigma_x, E_mat[[k]]) +
                   E_mat[[k]] %*% Sigma_x %*% t(B))
    psi_j_b[[k]] <- rbind(row1, row2)
  }

  psi_j_inv_b <- list()
  for (k in 1:(p * q)) {
    psi_j_inv_b[[k]] <- -psi_inv %*% psi_j_b[[k]] %*% psi_inv
  }

  partial_delta_j_b <- list()
  for (k in 1:(p * q)) {
    partial_delta_j_b[[k]] <-
      partial_delta_j(psi_j_inv_b[[k]], eta_j_b[[k]])
  }

  score_b <- matrix(, nrow = n, ncol = (q * p))
  for (k in 1:(q * p)) {
    score_b[, k] <- score(psi_j_b[[k]], partial_delta_j_b[[k]])
  }

  # score phi----
  eta_j_phi <- matrix(rep(0, r), ncol = 1)
  psi_j_phi <- diag(r)
  psi_j_inv_phi <- -psi_inv %*% psi_j_phi %*% psi_inv
  partial_delta_j_phi <- partial_delta_j(psi_j_inv_phi, eta_j_phi)
  score_phi <- score(psi_j_phi, partial_delta_j_phi)


  # score mu----
  score_mu <- matrix(, nrow = n, ncol = p)
  for (i in 1:p) {
    eta_j_mu <- matrix(c(diag(p)[, i], c(B %*% diag(p)[, i])), ncol = 1)
    psi_j_mu <- matrix(0, nrow = r, ncol = r)
    psi_j_inv_mu <- -psi_inv %*% psi_j_mu %*% psi_inv
    partial_delta_j_mu <- partial_delta_j(psi_j_inv_mu, eta_j_mu)
    score_mu[, i] <- score(psi_j_mu, partial_delta_j_mu)
  }

  # score vech(Sigma)----
  E_mat_Sigma <- list()
  for (l in 1:p) {
    E_mat_Sigma[[l]] <- list()
    for (k in 1:p) {
      if (k >= l) {
        E_aux =  matrix(0, nrow = p, ncol = p)
        E_aux[k, l] =  1
        E_aux[l, k] =  1
        E_mat_Sigma[[l]][[k]] <- E_aux
      }
    }
  }

  eta_j_Sigma <- list()
  for (l in 1:p) {
    eta_j_Sigma[[l]] <- list()
    for (k in 1:p) {
      if (k >= l) {
        eta_j_Sigma[[l]][[k]] <- matrix(rep(0, r), ncol = 1)
      }
    }
  }

  psi_j_Sigma <- list()
  for (l in 1:p) {
    psi_j_Sigma[[l]] <- list()
    for (k in 1:p) {
      if (k >= l) {
        psi_j_Sigma[[l]][[k]] <-
          Lambda %*% E_mat_Sigma[[l]][[k]] %*% t(Lambda)
      }
    }
  }

  psi_j_inv_Sigma <- list()
  for (l in 1:p) {
    psi_j_inv_Sigma[[l]] <- list()
    for (k in 1:p) {
      if (k >= l) {
        psi_j_inv_Sigma[[l]][[k]] <-
          -psi_inv %*% psi_j_Sigma[[l]][[k]] %*% psi_inv
      }
    }
  }

  partial_delta_j_Sigma <- list()
  for (l in 1:p) {
    partial_delta_j_Sigma[[l]] <- list()
    for (k in 1:p) {
      if (k >= l) {
        partial_delta_j_Sigma[[l]][[k]] <-
          partial_delta_j(psi_j_inv_Sigma[[l]][[k]], eta_j_Sigma[[l]][[k]])
      }
    }
  }

  score_Sigma <- list()
  for (l in 1:p) {
    score_Sigma[[l]] <- list()
    for (k in 1:p) {
      if (k >= l) {
        score_Sigma[[l]][[k]] <-
          score(psi_j_Sigma[[l]][[k]], partial_delta_j_Sigma[[l]][[k]])
      }
    }
  }
  score_Sigma <-
    matrix(unlist(score_Sigma),
           ncol = p * (p + 1) / 2,
           byrow = FALSE)

  # score matrix----
  score_theta <- cbind(score_a, score_b, score_phi, score_mu, score_Sigma)
  colnames(score_theta) <- nms(p, q)

  if (total == TRUE) {
    score_t <- colSums(score_theta)
    return(score_t)
  } else{
    return(score_theta)
  }
}
