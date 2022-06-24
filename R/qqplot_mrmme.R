#' Qqplot of the trasnformed distances
#'
#' @param X
#' @param Y
#' @param theta
#' @param nsim by default
#' @param conf by default
#' @param plot.all
#'
#' @examples
qqplot_Nmrmme <-
  function(X,
           Y,
           theta,
           nsim = 50,
           conf = 0.95,
           plot.all = FALSE) {
    X <- as.matrix(X)
    Y <- as.matrix(Y)
    Z <- cbind(X, Y)
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
    alpha <- c(rep(0, p), a)
    eta <- c(alpha + Lambda %*% mu_x)
    psi <- Lambda %*% tcrossprod(Sigma_x, Lambda) + diag(phi, r, r)
    psi <- 0.5 * (psi + t(psi))
    mu_Y <- a + B %*% matrix(mu_x, ncol = 1)
    Sigma_Y <- B %*% Sigma_x %*% t(B) + phi * diag(1, q, q)
    Sigma_Y <- 0.5 * (Sigma_Y + t(Sigma_Y))


    wh <- data.frame(
      N_wilson_hilferty(X, mu_x, Sigma_x),
      N_wilson_hilferty(Y, mu_Y, Sigma_Y),
      N_wilson_hilferty(Z, eta, psi)
    )
    colnames(wh) <- c("X", "Y", "Z")
    bands <- envlp(X, Y,
                   theta,
                   nsim = nsim,
                   conf = conf)
    bandX <- bands$bandX
    bandY <- bands$bandY
    bandZ <- bands$bandZ

    #par(mfrow=c(1,3))
    #Grafico de Z
    plots <- list()
    plots[["Z"]]  <- ggplot2::ggplot(data.frame(cbind(wh$Z, bandZ))) +
      ggplot2::stat_qq(ggplot2::aes(sample = X1), size = 0.8) +
      ggplot2::geom_qq(ggplot2::aes(sample = X2), col = "red", geom = "path") +
      ggplot2:: geom_qq(ggplot2::aes(sample = X3), col =  "red", geom = "path") +
      ggplot2::geom_abline(slope = 1, intercept = 0, col = "gray40", size=0.3) +
      ggplot2::theme_bw() +
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank()) +
      ggplot2::labs(title = "Z")

    plots[["X"]] <- ggplot2::ggplot(data.frame(cbind(wh$X, bandX))) +
      ggplot2::stat_qq(ggplot2::aes(sample = X1), size = 0.8) +
      ggplot2::geom_qq(ggplot2::aes(sample = X2), col = "red", geom = "path") +
      ggplot2::geom_qq(ggplot2::aes(sample = X3), col =  "red", geom = "path") +
      ggplot2::theme_bw() +
      ggplot2::labs(title = "X")

    plots[["Y"]] <- ggplot2::ggplot(data.frame(cbind(wh$Y, bandY))) +
      ggplot2::stat_qq(ggplot2::aes(sample = X1), size = 0.8) +
      ggplot2::geom_qq(ggplot2::aes(sample = X2), col = "red", geom = "path") +
      ggplot2::geom_qq(ggplot2::aes(sample = X3), col =  "red", geom = "path") +
      ggplot2::theme_bw() +
      ggplot2::labs(title = "Y")

    if (plot.all) {
      library(patchwork)
      (plots$Z | (plots$X / plots$Y)) +
        patchwork::plot_annotation(title = "QQ-plot of the transformed distances")
    } else {
      plots$Z # +
      # patchwork::plot_annotation(title = "QQ-plot of the transformed distances")
    }
  }
