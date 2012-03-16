#' @param horizons a vector of points in the future to select
#' @import plyr
embed.series <- function(series, dimensions, horizons=1, rem.NA = FALSE) {
  if(length(class(series))==1) attr(series, "dim") <- c(length(series), 1)
    output <- alply(series, 2, function(z) {
    embedding <- matrix(NA, nrow = length(z), ncol = 1+dimensions+length(horizons))
    embedding[,1] = seq(from = tsp(z)[1], to = tsp(z)[2], by = tsp(z)[3])
    for (i in 1:dimensions) {
      embedding[,i+1] <- c(rep(NA, i-1), z[1:(length(z)-i+1)])
      }
    for (i in 1:length(horizons)) {
      embedding[,1+dimensions+i] = c(z[(horizons[i]+1):length(z)], rep(NA, horizons[i]))
    }
    if(rem.NA) embedding <- embedding[complete.cases(embedding),]
    colnames(embedding) <- c("Time", paste(rep("D", dimensions), 1:dimensions,
        sep=""), paste(rep("H", length(horizons)), 1:length(horizons), sep=""))
    class(embedding) <- "embedding"
    attr(embedding, "dimensions") <- dimensions
    attr(embedding, "horizons") <-horizons
    return(embedding)
  })
  names(output) <- colnames(series)
  attr(output, "split_type") <- NULL
  attr(output, "split_labels") <- NULL
  if(length(output)==1) output <- output[[1]]
  return(output)
  }

#' Embed a time series in multiple dimensions and eliminate invalid vectrs
#'
#' This function converts the time-series into a matrix of vectors
#' embedded in dimension E whose elements are the data points and their
#' time-lagged points at (t, t-1, ..., t-E+1)
#'
#' It also eliminates vectors that include missing data or project to missing 
#' data, creates a list of prejection points
#' 
#' @param series the time series, as a vector of numeric values
#' @param E an integer number of dimensions into which embed the data
#' @param steps the number of time steps to project into the future.  
#' If steps = 0, no list of future values is returned.
#' @param cleanup an option of whether to remove vectors with missing data.  
#' If zero, no cleanup, if 1, removes vectors with missing values or with future 
#' projected missing values.  If 2, removes vectors that have missing values 
#' between themselves and projected futures.  (For instance, for time series 
#' with irregular sampling intervals)
#' @return a list consisting of matrix of of size \code{(length(series) - E + 1)}
#' by \code{E} ($cloud), and a vector of future values ($futures)
#'
#' @examples
#' data <- 1:100
#' embedded <- embed.series(data, 3)
#' embedded
embed.series <- function(series, E, steps = 0, cleanup = 0, changes=TRUE) {
  embedded <- matrix(NA, nrow = (length(series) - E + 1 - steps), ncol = (E+steps))
  for (i in 1:(dim(embedded)[1])) {
    embedded[i, ] = series[i:(i + E + steps - 1)]
  }
  embedded <- unique(embedded)
  if (cleanup == 1) {
    embedded <- as.matrix(embedded[which(!is.na(rowSums(embedded))), ])
  }
  if (cleanup == 2 ) {
    embedded <- as.matrix(embedded[which(!is.na(rowSums(cbind(embedded[, 1:E],
                                                              embedded[, E+steps])))), ])
  }
  
  training <- list(cloud = as.matrix(embedded[, 1:E]))
  if(steps > 0) training$futures <- embedded[, (E+steps)]
  if(changes) training$changes <- training$futures - training$cloud[,E]
  return(training)
}
