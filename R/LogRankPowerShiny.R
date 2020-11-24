
#' @name Log_Rank_Power

#' @title Two sample Log Rank Power Analysis

#' @description This function will run Log Rank Power Analysis Shiny application

#' @export



Log_Rank_Power <- function() {

  appDir <- system.file("shiny-examples", "myapp", package = "LogRankPower")

  if (appDir == "") {

    stop("Could not find example directory. Try re-installing `LogRankPower`.", call. = FALSE)

  }



  shiny::runApp(paste(appDir,'/LogRankPower_app.R',sep=''), launch.browser =T,host = getOption( "127.0.0.1"))

}
