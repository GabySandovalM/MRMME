#' Fitting a multivariate regression model with measurement errors
#'
#' `mrmme` is used to fit an estructural multivariate regression model when the covariates are measured
#' with error.
#'
#' @param Y a matrix or a data frame of \emph{n} rows and \emph{q} columns, where \emph{n} is the number of
#' observations and \emph{q} is the number of response variables.
#' @param X a matrix or a data frame of \emph{n} rows and \emph{p} columns, where \emph{n} is the number of
#' observations and \eqn{p} is the number of covariates.
#' @param method the method to be used in parameter estimation. Two are available:
#' `'ChanMak'` (the default) corresponds to estimators based on Chan and Mak (1984),
#' and `'EM'` that run the EM algorithm for \strong{MRMMEs} proposed by.
#' @param se.type the type of estimator of the Fisher information matrix to compute the standard
#' errors of the parameter estimates. Three are available: `'empirical'` (the default), `'expected'`, and `'observed'`.
#'
#' @param ... For `mrmme()` other arguments could be passed (see below).
#'
#' @param crit a numeric value giving the convergence criterion when `'EM'` algorithm is used for parameter estimation.
#' The default is 1e-10.
#'
#' @details
#'  Agregar
#'
#' @return An object of class `mrmme`. The output contains the parameter estimates,
#' standard errors, metrics of the model, and the data.
#' @export
#'
#' @references
#'
#' @seealso
#' \code{\link{autoplot.mrmme}}, \code{\link{influence.mrmme}}, \code{\link{residuals.mrmme}}
#'
#' @examples
#' X <- cbind(lung$X1,lung$X2)
#' Y <- cbind(lung$Y1,lung$Y2)
#' fit <- mrmme(X=X, Y=Y)
#' class(fit)
#' names(fit)

mrmme <- function(Y, X, method = "ChanMak", se.type = "empirical", ...) {

  Y = as.matrix(Y)
  X = as.matrix(X)


  # Validations -------------------------------------------------------------

  if(nrow(Y) != nrow(X)) stop("Number of rows in Y and X must coincide.")
  if(any(is.na(Y))) stop("There are some NAs values in Y matrix.")
  if(any(is.na(X))) stop("There are some NAs values in the X design matrix.")

  if(method == "ChanMak"){
    fit <- fit_ChanMak(X,Y)
  }else{
    if(method == "EM"){
      fit <- fit_EM(X,Y,...)
    }else{
      stop("The only two methods of estimation are 'ChanMak' and 'EM'")
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
        stop("The only three methods for SEs computation are 'empirical',
             'expected' and 'observed'")
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
                      "Pr(>|z|)" = stats::pnorm(q = abs(fit$theta)/se.fit$se,
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
  cat('Number of iterations:', fit$iter, '\n')

  # Output ------------------------------------------------------------------


  #list_out = list(fit = fit) #elegir obj para el usuario
  #out = list(fit = list_out)
  list_out <- list("a" = fit$a,
                   "B" = fit$B,
                   "phi" = fit$phi,
                   "mu.x" = fit$mu_x,
                   "Sigma.x" = fit$Sigma_x,
                   #"coef" = fit$coef,
                   "theta" = fit$theta,
                   "se.fit" = se.fit,
                   "se.type" = se.type,
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
             #"coef" = fit$coef,
             "theta" = fit$theta,
             "se.fit" = se.fit,
             "se.type" = se.type,
             "AIC" = fit$AIC,
             "BIC" = fit$BIC,
             "logL" = fit$logL,
             "X" = fit$X,
             "Y" = fit$Y)
  class(list_out) <- "mrmme"

  invisible(list_out)
}
