#' Forecast forward using neighbors
#'
#' This functions 
#'
#' @param indices a matrix of neighbors whose rows correspond to the target 
#' values, and whose indices are the rows of vectors in the training data
#' @param distances a matrix of the same size as \code{indices}, containing the 
#' distance from the target vector to the training vector referenced in
#' \code{indices}
#' @param futures a vector of for values of each training vector, of the value
#' in the future which is being estimated
#' @param weighting a parameter for exponential weighting of the importance of
#' neighbors for the projects.  \code{weighting = 0} means all neighbors have
#' equal value, while higher values favor closer neighbors over farther ones.

forecast <- function(indices, distances, futures, weighting = 1, scale=FALSE) {
  if (nrow(indices) != nrow(distances) || ncol(indices) != ncol(distances)) {
    stop("index and distance matrices must be of same dimension")
  }
  
  if(!scale) {
    weights <- t(apply(distances, 1, function(z) exp(-weighting * z)))
    }
  else {
    weights <- t(apply(distances, 1, function(z) exp(-weighting * z/scale)))
  }
  forecasts <- sapply(1:nrow(indices), 
                  function(z) weighted.mean(futures[indices[z,]], weights[z,]))
  return(list(forecasts=forecasts, weights=weights))
}

#'Forecast based on the whole cloud of points  

forecast.smaps <- function(targets, cloud, futures, theta=1) {

    distmat <- as.matrix(dist(rbind(targets, cloud),diag=TRUE, upper=TRUE ))
    allindices <- matrix(1:nrow(cloud),ncol=nrow(cloud),nrow=nrow(targets),byrow=TRUE)
    alldistances <- mat(distmat[1:nrow(targets), (nrow(targets) + 1):ncol(distmat)])
    scale <- mean(dist(cloud))

    changes <- futures-cloud[, ncol(cloud)]
    
    predict <- forecast(allindices, alldistances, changes, theta, scale)

    predict$forecasts <- predict$forecasts + targets[, ncol(targets)]
    
    return(predict)
  }
