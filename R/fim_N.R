#' Title
#'
#' @param theta
#' @param X
#' @param Y
#' @param type
#'
#' @return
#'
#' @examples
fim_N <- function(theta, X, Y, type = "emp") {
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  n <- dim(Y)[1]
  q <- dim(Y)[2]
  p <- dim(X)[2]
  r <- p + q
  d <- (p + 1) * (2 + p + 2 * q) / 2

  # expected ----
  if (type == "exp") {
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
    alpha <- c(rep(0, p), a)
    eta <- c(alpha + Lambda %*% mu_x)
    psi <- Lambda %*% tcrossprod(Sigma_x, Lambda) + diag(phi, r, r)
    psi <- 0.5 * (psi + t(psi))
    psi_inv <-  Rfast::spdinv(psi)

    # aux's
    V <-
      rbind(matrix(0, nrow = p, ncol = q), diag(nrow = q, ncol = q))
    A <- diag(nrow = r, ncol = r)

    # block a
    cov_a <- n * crossprod(V, psi_inv) %*% V

    # block B
    E_mat <- list()
    for (j in 1:(q * p)) {
      E_mat[[j]] <- list()
      for (l in 1:p) {
        E_mat[[j]][[l]] <- list()
        for (k in 1:q) {
          E_aux =  matrix(0, nrow = q, ncol = p)
          E_aux[k, l] =  1
          E_mat[[j]][[l]][[k]] <-
            rbind(diag(0, nrow = p, ncol = p), E_aux)
        }
      }
    }

    E_mat_aux <- unlist(E_mat[[1]], recursive = FALSE)

    cov_B <- matrix(0, nrow = q * p, ncol = q * p)
    for (i in 1:(q * p)) {
      izq <- E_mat_aux[[i]]
      fim_aux <- NULL
      for (k in 1:p) {
        for (l in 1:q) {
          der <- E_mat[[i]][[k]][[l]]###
          p1 <- crossprod(izq, psi_inv) %*% der
          s1 <- crossprod(mu_x, p1) %*% mu_x
          p2 <- Lambda %*% tcrossprod(Sigma_x, izq)
          p3 <- izq %*% tcrossprod(Sigma_x, Lambda)
          p4 <- Lambda %*% tcrossprod(Sigma_x, der)
          s2 <- psi_inv %*% (p2 + p3) %*% psi_inv %*% p4
          aux <- n * (s1 + sum(diag(s2)))
          fim_aux <- c(fim_aux, aux)
        }
      }
      cov_B[, i] <- fim_aux
    }

    # block phi
    cov_phi <- (n / 2) * sum(diag(psi_inv %*% A %*% psi_inv %*% A))

    # block mu
    cov_mu <- n * crossprod(Lambda, psi_inv) %*% Lambda

    # block Sigma
    E_mat_s <- list()
    for (j in 1:(p * p)) {
      E_mat_s[[j]] <- list()
      for (l in 1:p) {
        E_mat_s[[j]][[l]] <- list()
        for (k in 1:p) {
          if (k >= l) {
            E_aux =  matrix(0, nrow = p, ncol = p)
            E_aux[k, l] =  1
            E_aux[l, k] =  1
            E_mat_s[[j]][[l]][[k]] <- E_aux
          }
        }
      }
    }

    E_mat_aux_s <- unlist(E_mat_s[[1]], recursive = FALSE)
    E_mat_s_aux <- E_mat_aux_s[lengths(E_mat_aux_s) != 0]

    p_aux <- length(matrixcalc::vech(Sigma_x))
    cov_S <- matrix(0, nrow = p_aux, ncol = p_aux)
    for (i in 1:p_aux) {
      izq <- E_mat_s_aux[[i]]
      fim_aux <- NULL
      for (k in 1:p) {
        for (l in 1:p) {
          if (k >= l) {
            der <- E_mat_s[[i]][[l]][[k]]
            p1 <- crossprod(Lambda, psi_inv) %*% Lambda
            s1 <- izq %*% p1 %*% der %*% p1
            aux <- (n / 2) * sum(diag(s1))
            fim_aux <- c(fim_aux, aux)
          }
        }
      }
      cov_S[i, ] <- fim_aux
    }

    # block a with b
    cov_a_B <- matrix(0, nrow = q, ncol = (q * p))
    for (i in 1:(q * p)) {
      cov_a_B[, i] <-
        n * crossprod(V, psi_inv) %*% E_mat_aux[[i]] %*% mu_x
    }

    # block a with phi
    cov_a_phi <- matrix(0, nrow = q, ncol = 1)

    # block a with mu
    cov_a_mu <- n * crossprod(V, psi_inv) %*% Lambda

    # block a with Sigma
    cov_a_S <- matrix(0, nrow = q, ncol = p_aux)

    # block B with phi
    cov_B_phi <- matrix(0, nrow = (q * p), ncol = 1)
    for (i in 1:(q * p)) {
      p1 <- Sigma_x %*% crossprod(E_mat_aux[[i]], psi_inv)
      s1 <- A %*% psi_inv %*% Lambda %*% p1
      cov_B_phi[i, ] <- n * sum(diag(s1))
    }

    # block B with mu
    cov_B_mu <- matrix(0, nrow = (q * p), ncol = p)
    for (i in 1:(q * p)) {
      p1 <- crossprod(Lambda, psi_inv)
      p2 <- E_mat_aux[[i]] %*% mu_x
      cov_B_mu[i, ] <- n * p1 %*% p2
    }

    # block B with Sigma
    cov_B_S <- matrix(0, nrow = (q * p), ncol = p_aux)
    for (i in 1:(q * p)) {
      cov_aux <- NULL
      for (j in 1:p_aux) {
        p1 <- E_mat_s_aux[[j]] %*% crossprod(Lambda, psi_inv)
        p2 <-
          Lambda %*% Sigma_x %*% crossprod(E_mat_aux[[i]], psi_inv) %*% Lambda
        aux <- n * sum(diag(p1 %*% p2))
        cov_aux <- c(cov_aux, aux)
      }
      cov_B_S[i,] <- cov_aux
    }

    # block phi with mu
    cov_phi_mu <- matrix(0, nrow = 1, ncol = p)

    # block phi with Sigma
    cov_phi_S <- NULL
    for (j in 1:p_aux) {
      p1 <- psi_inv %*% Lambda
      p2 <- E_mat_s_aux[[j]] %*% crossprod(Lambda, psi_inv) %*% A
      cov_phi_S[j] <- (n / 2) * sum(diag(p1 %*% p2))
    }
    cov_phi_S <- matrix(cov_phi_S, nrow = 1)

    # block mu with Sigma
    cov_mu_S <- matrix(0, nrow = p, ncol = p_aux)

    # join together by cols
    a <- cbind(cov_a, cov_a_B, cov_a_phi, cov_a_mu, cov_a_S)
    B <- cbind(t(cov_a_B), cov_B, cov_B_phi, cov_B_mu, cov_B_S)
    phi <-
      c(t(cov_a_phi), t(cov_B_phi), cov_phi, cov_phi_mu, cov_phi_S)
    mu <-
      cbind(t(cov_a_mu), t(cov_B_mu), t(cov_phi_mu), cov_mu, cov_mu_S)
    S <-
      cbind(t(cov_a_S), t(cov_B_S), t(cov_phi_S), t(cov_mu_S), cov_S)

    # join together by rows
    fim_CM <- rbind(a, B, phi, mu, S)
    fim_CM <- 0.5 * (fim_CM + t(fim_CM)) # this is the FIM total
    colnames(fim_CM) <- nam(p, q)
    rownames(fim_CM) <- nam(p, q)
    fim_CM # this is the FIM total

  } else {
    # observed ----
    if (type == "obs") {
      fim_obs <- -numDeriv::jacobian(
        func = score_matrix_N,
        x = theta,
        X = X,
        Y = Y,
        total = TRUE
      )
      colnames(fim_obs) <- nam(p, q)
      rownames(fim_obs) <- nam(p, q)
      fim_obs <-  0.5 * (fim_obs + t(fim_obs))
      fim_obs # This is the FIM total

    } else {
      # empirical ----
      if (type == "emp") {
        score_mat <- score_matrix_N(theta, X, Y)
        fim_emp <- crossprod(score_mat)
        fim_emp <- 0.5 * (fim_emp + t(fim_emp))
        colnames(fim_emp) <- nam(p, q)
        rownames(fim_emp) <- nam(p, q)
        fim_emp # This is the FIM total
      }
    }
  }
}
