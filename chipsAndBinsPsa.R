# install.packages("SHELF")
# install.packages("jsonlite")
# library(SHELF)
# library(jsonlite)
# rm(list = ls())

source("./utilityFuncs.R")

chipsAndBinsPsa <- function(expertData = NULL, model = "best", psaSamples = 1000, xMinLimit = -Inf, xMaxLimit = Inf, json = TRUE) {
  
  require(SHELF)
  require(jsonlite)

  # filter out experts with no chips 
  filtered <- filterOutInvalidData(expertData)
  expertData <- filtered$expertData
  warningMsg <- filtered$warningMsg


  if(length(expertData) == 0){
    response <- list(success = FALSE, error = "All experts removed due to invalid data")
    return(jsonlite::toJSON(response, pretty = TRUE, auto_unbox = TRUE))
  }
  if(length(expertData) == 1){
    response <- list(success = FALSE, error = "Pooled estimation requires data from at least two experts")
    return(jsonlite::toJSON(response, pretty = TRUE, auto_unbox = TRUE))
  }

  psaSamples <- as.integer(psaSamples)

  # adjustment, making SHELF::fitdist() work for any chips allocation
  valuesAdj <- sapply(expertData, \(x){
    binWidth <- (x$xMax - x$xMin) / length(x$chips)
    xMinAdjusted <- x$xMin + binWidth/200
    xMaxAdjusted <- x$xMax - binWidth/200
    xTemp <- seq(xMinAdjusted, xMaxAdjusted, length.out = length(x$chips)+1)
    return(c(xTemp, x$xMax))
  })

  propsAdj <- sapply(expertData, \(expert){
    y <- expert$chips
    chipsTot <- sum(y) 
    temp <- c(0,y,0)
    chipSplit <- chipsTot/198 
    temp[which(temp > 0)[1]-1] <- chipSplit
    temp[which(temp > 0)[length(which(temp > 0))]+1] <- chipSplit
    return(cumsum(temp)/sum(temp))
  })

  weights <- sapply(expertData, \(x) ifelse(is.null(x$weight), 1, x$weight))

  # fit it!
  suppressWarnings({
    fitAdj <- SHELF::fitdist(
      vals = valuesAdj, 
      probs = propsAdj,
      lower = xMinLimit, 
      upper = xMaxLimit,
      weights = weights,
      expertnames = names(expertData)
    )

    sampledData <- rlinearpool(
      fit = fitAdj, 
      n = psaSamples,
      d = model, 
      w = weights
    )
    sampledData[sampledData < xMinLimit] <- xMinLimit
    sampledData[sampledData > xMaxLimit] <- xMaxLimit
  })
  
  # return response
  if (json) {
    response <- list(success = TRUE, data = sampledData)
    return(jsonlite::toJSON(response, pretty = TRUE, auto_unbox = TRUE))
  } else {
    return(sampledData)
  }
}


# # # # # Example usage:
# # Sample data
expertData <- list(
  alice = list(
    chips = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
    name = "expert1",
    xMin = 0,
    xMax = 1,
    weight = 1
  ),
  bob = list(
    chips = c(0, 0, 1, 3, 5, 3, 1, 0, 0, 0),
    name = "expert2",
    xMin = 0,
    xMax = 1,
    weight = 1
  ),
  claire = list(
    chips = c(0, 0, 0, 5, 5, 0, 0, 0, 0, 0),
    name = "expert3",
    xMin = 0,
    xMax = 1,
    weight = 1
  ),
  alex = list(
    chips = c(0, 0, 0, 0, 1, 2, 3, 4, 5, 5),
    name = "expert3",
    xMin = 0.5,
    xMax = 1,
    weight = 1
  ),
  george = list(
    chips = c(0, 0, 4, 5, 5, 4, 0, 0, 0, 0),
    name = "expert3",
    xMin = 0.5,
    xMax = 1,
    weight = 1
  ),
  william = list(
    chips = c(10, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    name = "expert3",
    xMin = 0,
    xMax = 1,
    weight = 1
  )
)


out2 <- chipsAndBinsPsa(expertData = expertData, model = "best", psaSamples = 10000, xMinLimit = 0, xMaxLimit = 1, json = FALSE)

# library(ggplot2)
ggplot() +
  geom_histogram(aes(out),bins = 100, alpha = 0.7, col = "lightgray", fill="cadetblue") +
  ggtitle("Sample Values for PSA") +
  xlab("Value") +
  ylab("Frequency") +
  theme_minimal()
