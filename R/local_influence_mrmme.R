#' Local influence diagnostics in MRMME
#'
#' @param theta a vector with the parameter estimates.
#' @param X covariates matrix. It has dimension n x p.
#' @param Y response variables matrix. It has dimension n x q.
#'
#' @return
#'
influence_mrmme <- function(theta, X, Y) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  Z <- cbind(X, Y)
  n <- dim(Y)[1]
  q <- dim(Y)[2]
  p <- dim(X)[2]
  r <- p + q
  d <- (p + 1) * (2 + p + 2 * q) / 2
  p_s <- p * (p + 1) / 2

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
  psi_inv <- Rfast::spdinv(psi)

  # m_i y M
  aux_m <- t(Z) - c(alpha + Lambda %*% mu_x)
  m <- t(mu_x + tcrossprod(Sigma_x, Lambda) %*% psi_inv %*% aux_m)
  M <- Sigma_x - tcrossprod(Sigma_x, Lambda) %*% psi_inv %*% Lambda %*% Sigma_x
  m_bar <- colMeans(m)
  I_q <- diag(1, q, q)
  D_p <- matrixcalc::duplication.matrix(p)
  S_mm <- 1 / n * tcrossprod(t(m) - m_bar)
  Sigma_x_inv <- Rfast::spdinv(Sigma_x)

  # Qpp ---------------------------------------------------------------------
  # second derivative of Q
  Q_a <- -(n / phi) * diag(1, q, q)
  Q_aB <- -(n / phi) * kronecker(t(m_bar), I_q)
  aux_B1 <-  kronecker(S_mm + tcrossprod(m_bar), I_q)
  aux_B2 <- kronecker(M, I_q)
  Q_B <- -(n / phi) * (aux_B1 + aux_B2)
  Q_phi <- matrix(-(n * r) / (2 * phi ^ 2), nrow = 1, ncol = 1)
  Q_mu <- -n * Sigma_x_inv
  aux_Sig1 <- kronecker(Sigma_x_inv %*% M %*% Sigma_x_inv, Sigma_x_inv)
  aux_Sig2 <- kronecker(Sigma_x_inv, Sigma_x_inv)
  Q_Sigma <- n * t(D_p) %*% (aux_Sig1 - (1 / 2 * aux_Sig2)) %*% D_p
  Qpp <- pracma::blkdiag(Q_a, Q_B, Q_phi, Q_mu, Q_Sigma)
  Qpp[1:q, (q + 1):(q * (1 + p))] <- Q_aB

  # Esquema 1 ---------------------------------------------------------------
  Delta_1 <- matrix(0, nrow = d, ncol = r)
  for (j in 1:r) {
    # e_j, c_qj
    e_j <- diag(r)[j, ]
    c_qj <- cbind(matrix(0, nrow = q, ncol = p), diag(1, q, q)) %*% e_j
    Delta_j <- matrix(0, nrow = d, ncol = 1)
    for (i in 1:n) {
      aux_delta <- c(t(Z[i, ] - alpha - Lambda %*% m[i, ]) %*% e_j)
      kron_prod <- kronecker(m[i, ], c_qj)
      Q_ia_omega  <- 1 / phi * aux_delta * c_qj
      Q_iB_omega  <- 1 / phi * aux_delta * kron_prod
      Q_iphi_omega  <- 1 / (2 * phi ^ 2) * aux_delta ^ 2
      Q_imu_omega  <- matrix(0, nrow = p, ncol = 1)
      Q_iSigma_omega  <- matrix(0, nrow = p_s, ncol = 1)
      aux_delta_j <-
        rbind(Q_ia_omega,
              Q_iB_omega,
              Q_iphi_omega,
              Q_imu_omega,
              Q_iSigma_omega)
      Delta_j <- Delta_j + aux_delta_j
    }
    Delta_1[, j] <- Delta_j
  }
  rownames(Delta_1) <- nms(p, q)


  # Esquema 2 ---------------------------------------------------------------

  Delta_2 <- matrix(0, nrow = d, ncol = n)
  for (i in 1:n) {
    aux_Y <- t(Y[i, ] - a - B %*% m[i, ])
    Q_ia_omega  <- 1 / phi * aux_Y
    Q_iB_omega  <- 1 / phi * (kronecker(m[i,], aux_Y) - c(B %*% M))
    # use ||A||^2 = tr(A^t A)
    aux_Z <- Z[i, ] - alpha - Lambda %*% m[i, ]
    S_i <-
      sum(diag(crossprod(aux_Z))) + sum(diag(M + t(B) %*% B %*% M))
    Q_iphi_omega  <- 1 / (2 * phi ^ 2) * S_i
    Q_imu_omega  <- matrix(0, nrow = p, ncol = 1)
    Q_iSigma_omega  <- matrix(0, nrow = p_s, ncol = 1)
    Delta_i <-
      c(Q_ia_omega,
        Q_iB_omega,
        Q_iphi_omega,
        Q_imu_omega,
        Q_iSigma_omega)
    Delta_2[, i] <- Delta_i
  }
  rownames(Delta_2) <- nms(p, q)

  # B_inf1 -------------------------------------------------------------------
  F_inf_1 <- t(Delta_1) %*% Rfast::spdinv(-Qpp) %*% Delta_1
  F_inf_1 <- 2 * abs(diag(F_inf_1)) / matrixcalc::matrix.trace(2 * F_inf_1)
  cut_point_1 <- mean(F_inf_1) + 2 * stats::sd(F_inf_1)

  # B_inf1 -------------------------------------------------------------------
  F_inf_2 <- t(Delta_2) %*% Rfast::spdinv(-Qpp) %*% Delta_2
  F_inf_2 <- 2 * abs(diag(F_inf_2)) / matrixcalc::matrix.trace(2 * F_inf_2)
  cut_point_2 <- mean(F_inf_2) + 2 * stats::sd(F_inf_2)

  inf_1 <- which(F_inf_1 > cut_point_1)
  inf_2 <- which(F_inf_2 > cut_point_2)

  # Plot scheme 1----
  ind1 <- 1:length(F_inf_1) #index
  p1 <-
    ggplot(data.frame(F_inf_1, ind1), aes(label = ind1)) +
    #scale_y_continuous(limits = c(0,1)) +
    geom_hline(aes(yintercept = cut_point_1), col = "red") +
    geom_point(aes(x = ind1, y = F_inf_1)) +
    labs(x = "", y = "F values", title = "Scheme I") +
    geom_text(aes(
      x = ind1,
      y = F_inf_1,
      label = ifelse(F_inf_1 >= cut_point_1, ind1, "")
    ),
    hjust = 0,
    vjust = 1,
    size = 2) +
    theme_bw()


  # Plot scheme 2----
  ind2 <- 1:length(F_inf_2) #index
  p2 <-
    ggplot(data.frame(F_inf_2, ind2), aes(label = ind2)) +
    #scale_y_continuous(limits = c(0,1)) +
    geom_hline(aes(yintercept = cut_point_2), col = "red") +
    geom_point(aes(x = ind2, y = F_inf_2)) +
    labs(x = "", y = "F values", title = "Scheme II") +
    geom_text(
      aes(
        x = ind2,
        y = F_inf_2,
        label = ifelse(F_inf_2 > cut_point_2, ind2, "")
      ),
      hjust = 0,
      vjust = 1,
      size = 2
    ) +
    theme_bw() +
    theme(panel.grid.minor.x = element_blank())

  p2t <- (p1 | p2) +
    patchwork::plot_annotation(title = "Influence diagnostics: Local influence")

  res = list(
    "F1" = F_inf_1,
    "F2" = F_inf_2,
    "cut_point_SCH1" = cut_point_1,
    "cut_point_SCH2" = cut_point_2,
    "influentials_SCH1" =  inf_1,
    "influentials_SCH2" =  inf_2,
   #"plot1" = p1,
   #"plot2" = p2,
    "plot" = p2t
  )

  return(res)

}

