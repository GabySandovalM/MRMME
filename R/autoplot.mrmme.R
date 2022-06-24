#' Plot model assessment for a MRMME
#'
#' @param model The \code{\link{mrmme}} object
#'
#' @return A \code{\link{ggplot}} object.
#'
#' @importFrom patchwork
#' @importFrom ggplot2 autoplot
#' @export autoplot.mrmme
#' @export
#' @examples
autoplot.mrmme <- function(model,...) {

# Validations -------------------------------------------------------------


# qqplot ------------------------------------------------------------------

  p1 <- GGally::ggpairs(data.frame(cbind(mod$Y,mod$X)), axisLabels = "none") +
        ggplot2::theme_bw()
  p2 <- qqplot_Nmrmme(X=mod$X, Y=mod$Y, theta = mod$theta, plot.all = TRUE)
  p3 <- qqplot_Nmrmme(X=mod$X, Y=mod$Y, theta = mod$theta, plot.all = FALSE)
  devAskNewPage(TRUE)
  print(p1)
  print(p2)
  print(p3)
  devAskNewPage(options("device.ask.default")[[1]])
}

