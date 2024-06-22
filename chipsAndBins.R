# install.packages("SHELF")
# install.packages("jsonlite")
# library(SHELF)
# library(jsonlite)
# rm(list = ls())

source("./utilityFuncs.R")

chipsAndBins <- function(
  expertData = NULL, 
  model = "best", 
  nx = 100, 
  xMinLimit = -Inf,
  xMaxLimit = Inf,
  json = TRUE
  ) {

    require(SHELF)
    require(jsonlite)


    filtered <- filterOutInvalidData(expertData)
    expertData <- filtered$expertData
    warningMsg <- filtered$warningMsg
  

    if(length(expertData) == 0){
      response <- list(success = FALSE, error = "All experts removed due to invalid data")
      return(jsonlite::toJSON(response, pretty = TRUE, auto_unbox = TRUE))
  }


    xMinGlobal <- min(sapply(expertData, \(x) x$xMin))
    xMaxGlobal <- max(sapply(expertData, \(x) x$xMax))
    xPredGlobal <- seq(xMinGlobal, xMaxGlobal, length.out = nx)
    binWidthGlobal <- (xMaxGlobal - xMinGlobal) / nx

  out <- list()
  for(i in 1:length(expertData)){
    expertId <- names(expertData)[i]
    d_ <- expertData[[i]]

    xBins <- length(d_$chips)
    edges <- seq(from = d_$xMin, to = d_$xMax, length.out = xBins + 1)
    upperBinEdges <- edges[-1]
    lowerBinEdges <- edges[-length(edges)]
    binWidth <- (d_$xMax - d_$xMin) / xBins  
    probs <- cumsum(d_$chips) / sum(d_$chips)
    
    histAdjFactor <- binWidth / binWidthGlobal
    yHist <- diff(c(0,probs)) / histAdjFactor
      
  
    out[[expertId]] <- list(
      expertId = expertId,
      expertName = d_$name,
      chips = d_$chips,
      xBins = xBins,
      binWidth = binWidth,
      xMin = d_$xMin,
      xMax = d_$xMax,
      binEdges = edges,
      upperBinEdges = upperBinEdges,
      lowerBinEdges = lowerBinEdges,
      rawResponseData = list(
        binNumber = 1:xBins,
        chips = d_$chips,
        binEdges = paste0(lowerBinEdges, "-", upperBinEdges),
        probs = probs
        ),
      hist = list(
        x = c(edges[1], edges),
        y = c(0, yHist,0)
      ),
      xout = c(edges, d_$xMax + binWidth / 2),
      weight = ifelse(is.null(d_$weight), 1, d_$weight)
    )
  }

  # adjustment, making SHELF::fitdist() work for any chips allocation
  valuesAdj <- sapply(expertData, \(x){
    binWidth <- (x$xMax - x$xMin) / length(x$chips)
    xMinAdjusted <- x$xMin + binWidth/200
    xMaxAdjusted <- x$xMax - binWidth/200
    xTemp <- seq(xMinAdjusted, xMaxAdjusted, length.out = length(x$chips)+1)
    
    # return(xTemp)
    c(xTemp, x$xMax)
    # c(xTemp)
  })

  propsAdj <- sapply(expertData, \(expert){
    
    # t_ <- expert$chips
    # return(cumsum(t_)/sum(t_))

    y <- expert$chips
    chipsTot <- sum(y) 
    temp <- c(0,y,0)
    chipSplit <- chipsTot/198 
    temp[which(temp > 0)[1]-1] <- chipSplit
    temp[which(temp > 0)[length(which(temp > 0))]+1] <- chipSplit
    return(cumsum(temp)/sum(temp))
  })

  weights <- sapply(expertData, \(expert) expert$weight)

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

quants <- c(0.05, 0.25, 0.5, 0.75, 0.95)

# get fitted probabilities + quantiles
  fittedAdj <- SHELF:::feedbackgroup(
    fit = fitAdj, 
    quantiles = quants,
    values = xPredGlobal, 
    dist = model, 
    sf = 10
  )
  })

  
  
  for(i in seq_along(expertData)){
    
    expertId <- names(expertData)[i]
    
    out[[expertId]]$fitBest <- fitAdj$best.fitting[i,]
    fitUsed <- fittedAdj$distribution[1,i]
    out[[expertId]]$fitUsed <- fitUsed
    out[[expertId]]$modelCoefs <- as.list(fitAdj[[translateDist(fitUsed)]][expertId,])

    out[[expertId]]$quantiles <- list(
      quantile = quants,
      value = fittedAdj$fitted.quantiles[,i]
    )
    out[[expertId]]$fit <- list(
      x = xPredGlobal,
      y = c(0,diff(fittedAdj$fitted.probabilities[,i]))
    )
  }

    
    # linear pool: hist + fitted
    yHistPoolCum <- SHELF::plinearpool(
      fit = fitAdj, 
      x = xPredGlobal, 
      d = "best", 
      w = weights
      )
    yHistPool <- c(0,diff(yHistPoolCum))  
    # # no need for plinearpool, because this is identical:
    yLinearPool <- apply(sapply(out, \(x) x$fit$y),1,weighted.mean, w = weights)

  out[['linearPool']] <- list(
    expertId = "linearPool",
    expertName = "Linear Pool",
    fitUsed = ifelse(model == "best", "linearPool", model),
    fit = list(
      x = xPredGlobal,
      y = yLinearPool
    ),
    hist = list(
      x = xPredGlobal,
      y = yHistPool
    )
  )    

#### returning response
  
  if (json) {
    response <- list(success = TRUE, data = out, warning = warningMsg)
    return(jsonlite::toJSON(response, pretty = TRUE, auto_unbox = TRUE))
  } else {
    # return the results as a list for testing
    return(out)
  }

}



# ## EXAMPLE USAGE
plotIt <- function(response, hist = T, fit = T, pooled = F, facet = F){
  linearPool <- response[names(response) == "linearPool"]
  out2 <- response[names(response) != "linearPool"]
  require(ggplot2)
  xPredGlobal<- seq(
    from = min(sapply(out2, \(x) x$xMin)),
    to = max(sapply(out2, \(x) x$xMax)),
    length.out = nx
  )
  dfPlot <- data.frame(
        x = unlist(lapply(out2, \(x) x$fit$x)),
        y = unlist(lapply(out2, \(x) x$fit$y)),
        group = rep(names(out2), each = length(xPredGlobal)),
        type = "fitted"
      )
    dfHist <- data.frame(
      x = unlist(lapply(out2, \(x) x$hist$x)),
      y = unlist(lapply(out2, \(x) x$hist$y)),
      group = rep(names(out2), each = length(out2[[1]]$hist$x)),
      type = "hist"
    )
    dfPool <- data.frame(
      x = unlist(lapply(linearPool, \(x) x$fit$x)),
      y = unlist(lapply(linearPool, \(x) x$fit$y)),
      group = rep(names(linearPool), each = length(xPredGlobal)),
      type = "fitted"
    )
  
    p1 = ggplot() 
    if(fit){
      p1 = p1 + geom_line(
          data = dfPlot,
          aes(x = x, y = y, color = group),
          size = ifelse(pooled, 0.7, 1.5)
      ) 
    }
    if(hist){
      p1 = p1 + geom_step(
          data = dfHist,
          aes(x = x, y = y, color = group),
          direction = "hv"
      ) 
    }
    if(pooled){
      p1 <- p1 + geom_line(
        data = dfPool,
        aes(x = x, y = y, color = group),
        colour = "black",
        size = 3
      ) 
    }
    if(facet){
      p1 = p1  + facet_wrap(~group)  #, ncol = 3) +
    }
    p1 + theme_minimal() +
    # scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
    theme(legend.position = "top")
}


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

  
out <- chipsAndBins(
   expertData = expertData, 
   model = "best", 
   nx = 100, 
   xMinLimit = 0, 
   xMaxLimit = 1, 
   json = FALSE
  )
 
plotIt(out, hist = T, fit = T, pooled = F, facet = T)
plotIt(out, hist = F, fit = T, pooled = T)

# individual expert data
str(out$alice)

# linear pool data
str(out$linearPool)
