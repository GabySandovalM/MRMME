#' QQ-plot of the transformed distances.
#'
#' @param X covariates matrix. It has dimension n x p.
#' @param Y response variables matrix. It has dimension n x q.
#' @param theta a vector with the parameter estimates.
#' @param nsim A numeric value that indicate the number of simulation used to compute the envelopes.
#' By default is 50.
#' @param conf A numeric value between 0 and 1, that indicate the confidence level used to compute the envelopes.
#' By default is 0.95
#' @param plot.all if `TRUE` displays the qq-plot of Z, Y and X. If `FALSE` (the default) displays the qq-plot of Z.
#'
#' @import ggplot2
#' @import patchwork
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
      X = N_wilson_hilferty(X, mu_x, Sigma_x),
      Y = N_wilson_hilferty(Y, mu_Y, Sigma_Y),
      Z = N_wilson_hilferty(Z, eta, psi)
    )
    bands <- envlp(X, Y,
                   theta,
                   nsim = nsim,
                   conf = conf)
    bandX <- bands$bandX
    bandY <- bands$bandY
    bandZ <- bands$bandZ

    #Grafico de Z
    data_plot_Z <- data.frame(cbind(wh$Z, bandZ))
    colnames(data_plot_Z) <- c("Z.wh","Z.env.lo", "Z.env.up")
    data_plot_X <- data.frame(cbind(wh$X, bandX))
    colnames(data_plot_X) <- c("X.wh","X.env.lo", "X.env.up")
    data_plot_Y <- data.frame(cbind(wh$Y, bandY))
    colnames(data_plot_Y) <- c("Y.wh","Y.env.lo", "Y.env.up")

    plots <- list()
    # plots[["Z"]]  <- ggplot2::ggplot(data_plot_Z) +
    #   ggplot2::stat_qq(ggplot2::aes(sample = .data$Z.wh), size = 0.8) +
    #   ggplot2::geom_qq(ggplot2::aes(sample = .data$Z.env.lo), col = "red", geom = "path") +
    #   ggplot2::geom_qq(ggplot2::aes(sample = .data$Z.env.up), col =  "red", geom = "path") +
    #   ggplot2::geom_abline(slope = 1, intercept = 0, col = "gray40", size=0.3) +
    #   ggplot2::theme_bw() +
    #   ggplot2::theme(panel.grid.minor = ggplot2::element_blank()) +
    #   ggplot2::labs(title = "Z")
    #
    # plots[["X"]] <- ggplot2::ggplot(data_plot_X) +
    #   ggplot2::stat_qq(ggplot2::aes(sample = .data$X.wh), size = 0.8) +
    #   ggplot2::geom_qq(ggplot2::aes(sample = .data$X.env.lo), col = "red", geom = "path") +
    #   ggplot2::geom_qq(ggplot2::aes(sample = .data$X.env.up), col =  "red", geom = "path") +
    #   ggplot2::theme_bw() +
    #   ggplot2::labs(title = "X")
    #
    # plots[["Y"]] <- ggplot2::ggplot(data_plot_Y) +
    #   ggplot2::stat_qq(ggplot2::aes(sample = .data$Y.wh), size = 0.8) +
    #   ggplot2::geom_qq(ggplot2::aes(sample = .data$Y.env.lo), col = "red", geom = "path") +
    #   ggplot2::geom_qq(ggplot2::aes(sample = .data$Y.env.up), col =  "red", geom = "path") +
    #   ggplot2::theme_bw() +
    #   ggplot2::labs(title = "Y")

    plots[["Z"]]  <- ggplot(data_plot_Z) +
      stat_qq(aes(sample = .data$Z.wh), size = 0.8) +
      geom_qq(aes(sample = .data$Z.env.lo), col = "red", geom = "path") +
      geom_qq(aes(sample = .data$Z.env.up), col =  "red", geom = "path") +
      geom_abline(slope = 1, intercept = 0, col = "gray40", size=0.3) +
      theme_bw() +
      theme(panel.grid.minor = element_blank()) +
      labs(title = "Z")

    plots[["X"]] <- ggplot(data_plot_X) +
      stat_qq(aes(sample = .data$X.wh), size = 0.8) +
      geom_qq(aes(sample = .data$X.env.lo), col = "red", geom = "path") +
      geom_qq(aes(sample = .data$X.env.up), col =  "red", geom = "path") +
      theme_bw() +
      labs(title = "X")

    plots[["Y"]] <- ggplot(data_plot_Y) +
      stat_qq(aes(sample = .data$Y.wh), size = 0.8) +
      geom_qq(aes(sample = .data$Y.env.lo), col = "red", geom = "path") +
      geom_qq(aes(sample = .data$Y.env.up), col =  "red", geom = "path") +
      theme_bw() +
      labs(title = "Y")


    if (plot.all) {
      if (!requireNamespace("patchwork", quietly = TRUE)) {
        stop(
          "Package \"patchwork\" must be installed to use this function.",
          call. = FALSE
        )
      }
      (plots$Z | (plots$X / plots$Y)) +
        patchwork::plot_annotation(title = "QQ-plot of the transformed distances")
    } else {
      plots$Z  +
      patchwork::plot_annotation(title = "QQ-plot of the transformed distances")
    }
  }
