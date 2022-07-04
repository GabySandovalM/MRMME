#' Generate a vector with parameter names.
#'
#' @param p a numeric value indicating the number of covariates.
#' @param q a numeric value indicating the number of response variables.
#' @param print if TRUE uses `[ ]` if it is FALSE uses `_` in the subscripts.
#'
#' @return a vector with the names of the parameters involved in the fitted MRMME.
nms <- function(p, q, print = FALSE) {
  if (print == TRUE) {
    name_theta <- NULL
    for (i in 1:q) {
      name_a <- paste0("a[", i, "]")
      name_theta <- c(name_theta, name_a)
    }

    for (k in 1:p) {
      for (l in 1:q) {
        name_b <- paste0("b[", l, k, "]")
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
  } else {
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
  }
  name_theta
}
