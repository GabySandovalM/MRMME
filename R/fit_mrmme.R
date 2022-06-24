#' Fitting a multivariate regression model with measurement errors
#'
#' @param Y
#' @param X
#' @param type
#' @param se.type
#'
#' @return
#' @export
#'
#' @examples
mrmme <- function(Y, X, type = "ChanMak", se.type = "empirical", ...) {

  Y = as.matrix(Y)
  X = as.matrix(X)


  # Validations -------------------------------------------------------------

  if(nrow(Y) != nrow(X)) stop("Number of rows in Y and X must coincide.")
  if(any(is.na(Y))) stop("There are some NAs values in Y matrix.")
  if(any(is.na(X))) stop("There are some NAs values in the X design matrix.")

  if(type == "ChanMak"){
    fit <- fit_ChanMak(X,Y)
  }else{
    if(type == "EM"){
      fit <- fit_EM(X,Y)
    }else{
      stop("The only two models available are 'ChanMak' and 'EM'")
    }
  }

  if(se.type == "empirical"){
    se.fit <- se_N(fit$theta,fit$X,fit$Y, type = "emp")
  }else{
    if(se.type == "expected"){
      se.fit <- se_N(fit$theta,fit$X,fit$Y, type = "exp")
    }else{
      if(se.type == "observed"){
        se.fit <- se_N(fit$theta,fit$X,fit$Y, type = "obs")
      }else{
        stop("The only three methods for SEs computation are...")
      }
    }
  }

  # Print -------------------------------------------------------------------

  cat('\n')
  call = match.call()
  cat("Call:\n")
  print(call)
  cat('\n')

  estim <- data.frame(Estimate = round(fit$theta,3),
                      "Std. Error" = paste0("(",round(se.fit$se,3),")"),
                      "z value" = round(fit$theta/se.fit$se,3),
                      "Pr(>|z|)" = pnorm(q = abs(fit$theta)/se.fit$se,
                                         lower.tail = FALSE),
                      check.names = FALSE)

  rownames(estim) <- nam2(dim(Y)[2],dim(X)[2])

  asteriscos = sapply(X = estim[,4],FUN = defast)
  estim      = data.frame(estim," " = asteriscos,
                          check.names = FALSE)

  cat('\n')
  # cat('---------\n')
  # cat('Estimates\n')
  cat('Coefficients:\n')
  # cat('\n')
  print(estim)
  cat('---\n')
  cat('Signif. codes:  0 "***" 0.001 "**" 0.01 "*" 0.05 "." 0.1 " " 1\n')
  cat('\n')

  cat('AIC:',fit$AIC, 'BIC:',fit$BIC,'\n')
  cat('Obs. loglik:', fit$logL, '\n')

  # Output ------------------------------------------------------------------


  #list_out = list(fit = fit) #elegir obj para el usuario
  #out = list(fit = list_out)
  list_out <- list("a" = fit$a,
                   "B" = fit$B,
                   "phi" = fit$phi,
                   "mu.x" = fit$mu_x,
                   "Sigma.x" = fit$Sigma_x,
                   "coef" = fit$coef,
                   "theta" = fit$theta,
                   "AIC" = fit$AIC,
                   "BIC" = fit$BIC,
                   "logL" = fit$logL,
                   "X" = fit$X,
                   "Y" = fit$Y)
  out <- list("a" = fit$a,
             "B" = fit$B,
             "phi" = fit$phi,
             "mu.x" = fit$mu_x,
             "Sigma.x" = fit$Sigma_x,
             "coef" = fit$coef,
             "theta" = fit$theta,
             "AIC" = fit$AIC,
             "BIC" = fit$BIC,
             "logL" = fit$logL,
             "X" = fit$X,
             "Y" = fit$Y)
  class(list_out) <- "mrmme"

  invisible(list_out)
}
