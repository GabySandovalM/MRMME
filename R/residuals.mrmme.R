#' Residuals in MRMME
#'
#' @param model
#'
#' @return
#' @importFrom stats residuals
#' @export residuals.mrmme
#' @export
#' @examples
residuals.mrmme <- function(model){

# validations -------------------------------------------------------------

  if(!inherits(model,"mrmme")) stop("The argument model must be an object of class mrmme.")

  Y <- as.matrix(model$Y)
  X <- as.matrix(model$X)
  Z <- cbind(X,Y)
  p <- dim(X)[2]
  q <- dim(Y)[2]
  r <- p+q
  Lambda <- rbind(diag(p),model$B)
  alpha <- c(rep(0,p),model$a)
  Psi <- Lambda %*% model$Sigma.x %*% t(Lambda) + model$phi*diag(r)
  Psi_inv <- solve(Psi)
  eta <- c(alpha + Lambda%*%model$mu.x)
  hat_x <- model$mu.x + model$Sigma.x %*% t(Lambda) %*% Psi_inv %*% (t(Z)-eta)
  hat_y <- c(model$a) + model$B %*% hat_x
  hat_x <- t(hat_x)
  hat_y <- t(hat_y)
  res <- list("hat.x" = hat_x,
              "hat.y" = hat_y)
  return(res)
}
