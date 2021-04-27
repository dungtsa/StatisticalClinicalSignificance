
#' @name StatClinSigApp

#' @title Two sample Log Rank Power Analysis

#' @description This function will run Log Rank Power Analysis Shiny application

#' @export



StatClinSigApp <- function() {

  appDir <- system.file("shiny-examples", "myapp", package = "StatisticalClinicalSignificance")

  if (appDir == "") {

    stop("Could not find example directory. Try re-installing `StatisticalClinicalSignificance`.", call. = FALSE)

  }



  shiny::runApp(paste(appDir,'/StatClinSigApp.R',sep=''), launch.browser =T,host = getOption( "127.0.0.1"))

}
