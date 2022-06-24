#' Influence measures in MRMME
#'
#' @param model
#'
#' @return
#' @importFrom stats influence
#' @export influence.mrmme
#' @export
#' @examples
influence.mrmme <- function(model){

# Validations -------------------------------------------------------------



# Influencia local --------------------------------------------------------

  ip <- influence_mrmme(model$theta, X = model$X, Y = model$Y)

  # Print -------------------------------------------------------------------
  cat('\n')
  cat('INFLUENCE DIAGNOSTICS')
  cat('\n')
  cat('\n')
  cat('Method: Local influence')
  cat('\n')
  cat('Potentially influential observations under Scheme 1:',ip$influentials_SCH1,'\n')
  cat('Potentially influential observations under Scheme 2:',ip$influentials_SCH2)
  cat('\n')
  cat('\n')
  #cat('Method: Case-Deletion')
  cat('\n')
  devAskNewPage(TRUE)
  print(ip$plot)
  devAskNewPage(options("device.ask.default")[[1]])

  # Output ------------------------------------------------------------------
  #list_out = list(im = ip) #elegir obj para el usuario
  #out = list(im = list_out)
  list_out <- list("influentials.SCH1" = ip$influentials_SCH1,
                   "influentials.SCH2" = ip$influentials_SCH2,
                   "F.values.SCH1" = ip$F1,
                   "F.values.SCH2" = ip$F2,
                   "cut.point.SCH1" = ip$cut_point_SCH1,
                   "cut.point.SCH2" = ip$cut_point_SCH2)
  out <- list("influentials.SCH1" = ip$influentials_SCH1,
              "influentials.SCH2" = ip$influentials_SCH2,
              "F.values.SCH1" = ip$F1,
              "F.values.SCH2" = ip$F2,
              "cut.point.SCH1" = ip$cut_point_SCH1,
              "cut.point.SCH2" = ip$cut_point_SCH2)

  invisible(list_out)

}
