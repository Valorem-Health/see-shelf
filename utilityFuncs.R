# Description: Utility functions for returning SHELF analyses for chips and bins data from multiple experts

filterOutInvalidData <- function(expertData) {
  
  exclNoChips <- 0
  exclNamesNoChips <- c()
  warningMsg <- ""

  ids_ <- names(expertData)
  for (id in ids_) {
    expertName <- expertData[[id]]$name
    if (is.null(expertData[[id]]$chips) || 
    length(expertData[[id]]$chips) == 0  || 
    sum(expertData[[id]]$chips) == 0) {
      exclNoChips <- exclNoChips + 1
      exclNamesNoChips <- c(expertName, exclNamesNoChips)
      expertData <- expertData[-which(names(expertData) == id)]
      next
    } 
  }

  if (exclNoChips > 0) {
    warningMsg <- paste0(warningMsg, exclNoChips, " expert(s) were excluded because chips data is missing: ", paste(exclNamesNoChips, collapse = ", "), "\n")
  }
  
  return(list(expertData = expertData, warningMsg = warningMsg))
}



# from best.fitting value to model name
translateDist <- function(bestFittingVal) {
  switch(bestFittingVal,
         "normal" = "Normal",
         "skewnormal" = "Skewnormal",
         "gamma" = "Gamma",
         "mirrorgamma" = "mirrorgamma",
         "lognormal" = "Log.normal",
         "mirrorlognormal" = "mirrorlognormal",
         "logt" = "Log.Student.t",
         "mirrorlogt" = "mirrorlogt",
         "beta" = "Beta",
         "t" = "Student.t",
  )
}