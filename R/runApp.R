#!/usr/bin env Rscript

#' @export
runExample <- function() {
  appDir <- system.file("application", package = "BIGDAWGv2")
  if (appDir == "") {
    stop("Could not find myapp. Try re-installing 'BIGDAWGv2'.", call. = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal")
}
