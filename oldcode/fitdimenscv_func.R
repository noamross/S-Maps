fit.dimensions.cv <- function(series, dimensions, horizon = 1) {
  embeds <- alply(dimensions, .margins=1, function(z) embed.series(series, dimensions=z, horizons=horizon, rem.NA=TRUE))  #Create a list of embeddings in different dimensions
  distmats <- 
  lsqerrs <- vector("list",length=length(dimensions))
  for (i in 1:length(embeds)) {  # now for each embedding in the list we will do a cross-validation dimensional fit
      sqerrs <- rep(NA, nrow(embeds[[i]]))
      for(j in 1:nrow(embeds[[i]])) { # Do a leave-one-out fit for every point
        target <- embeds[[i]][j,]; attributes(target) <- append(list(dim=c(1,ncol(embeds[[i]]))), attributes(embeds[[i]])[-1])
        cloud <- embeds[[i]][-j,]; attributes(cloud) <- append(attributes(cloud), attributes(embeds[[i]])[-(1:2)])
        simplex <- neighborhood(target, cloud, method="simplex")  #find the nearest simplex of predictor points to each target
        predictions <- forecast(simplex$indices, simplex$distances, cloud[, 2+attr(cloud, "dimensions")])
        sqerrs[j] <- (predictions$forecasts - target[, 2+attr(cloud, "dimensions")])^2
    }
  lsqerrs[[i]] <- sqerrs
  names(lsqerrs) <- paste(paste("D",dimensions, sep=""))
  } 
  msqerrs <- laply(lsqerrs, mean)
  names(msqerrs) <- paste(paste("D",dimensions, sep=""))
  return(c(list(msqerrs=msqerrs), lsqerrs))
}

#This version uses correlation coefficients rather than mean square error to determine goodness of fit.
fit.dimensions.cv2 <- function(series, dimensions, horizon = 1) {
  embeds <- alply(dimensions, .margins=1, function(z) embed.series(series, dimensions=z, horizons=horizon, rem.NA=TRUE))  #Create a list of embeddings in different dimensions
  corrs <- vector("list",length=length(dimensions))
  for (i in 1:length(embeds)) {  # now for each embedding in the list we will do a cross-validation dimensional fit
    predicts <- rep(NA, nrow(embeds[[i]]))
    reals <- predicts
    for(j in 1:nrow(embeds[[i]])) { # Do a leave-one-out fit for every point
      target <- embeds[[i]][j,]; attributes(target) <- append(list(dim=c(1,ncol(embeds[[i]]))), attributes(embeds[[i]])[-1])
      cloud <- embeds[[i]][-j,]; attributes(cloud) <- append(attributes(cloud), attributes(embeds[[i]])[-(1:2)])
      simplex <- neighborhood(target, cloud, method="simplex")  #find the nearest simplex of predictor points to each target
      predictions <- forecast(simplex$indices, simplex$distances, cloud[, 2+attr(cloud, "dimensions")])
      predicts[j] <- predictions$forecasts
      reals[j] <- target[, 2+attr(cloud, "dimensions")]
    }
    corrs[[i]] <- cbind(predicts, reals); colnames(corrs[[i]]) <- c("predicts","reals")
    names(corrs) <- paste(paste("D",dimensions, sep=""))
  } 
  corrcoeff <- laply(corrs, corr)
  names(corrcoeff) <- paste("D",dimensions, sep="")
  return(c(list(corrcoeff=corrcoeff), corrs))
}

#This version uses correlation coefficients rather than mean square error to determine goodness of fit.  It will also only compute the distance matrix once per embedding
fit.dimensions.cv3 <- function(series, dimensions, horizon = 1) {
  embeds <- alply(dimensions, .margins=1, function(z) embed.series(series, dimensions=z, horizons=horizon, rem.NA=TRUE))  #Create a list of embeddings in different dimensions
  distmats <- llply(embeds, function(z) as.matrix(dist(z,diag=TRUE, upper=TRUE)))
  corrs <- vector("list",length=length(dimensions))
  for (i in 1:length(embeds)) {  # now for each embedding in the list we will do a cross-validation dimensional fit
    predicts <- rep(NA, nrow(embeds[[i]]))
    reals <- predicts
    for(j in 1:nrow(embeds[[i]])) { # Do a leave-one-out fit for every point
      target <- embeds[[i]][j,]; attributes(target) <- append(list(dim=c(1,ncol(embeds[[i]]))), attributes(embeds[[i]])[-1])
      cloud <- embeds[[i]][-j,]; attributes(cloud) <- append(attributes(cloud), attributes(embeds[[i]])[-(1:2)])
      simplex <- neighborhood(target, cloud, method="simplex", cvdistmat=distmats[[i]], cvn=j)  #find the nearest simplex of predictor points to each target
      predictions <- forecast(simplex$indices, simplex$distances, cloud[, 2+attr(cloud, "dimensions")])
      predicts[j] <- predictions$forecasts
      reals[j] <- target[, 2+attr(cloud, "dimensions")]
    }
    corrs[[i]] <- cbind(predicts, reals); colnames(corrs[[i]]) <- c("predicts","reals")
    names(corrs) <- paste(paste("D",dimensions, sep=""))
  } 
  corrcoeff <- laply(corrs, corr)
  names(corrcoeff) <- paste("D",dimensions, sep="")
  return(c(list(corrcoeff=corrcoeff), corrs))
}

