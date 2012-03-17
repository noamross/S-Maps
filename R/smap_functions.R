# Call relevant packages
require(yaImpute)        											# For nearest-neighbor search function
require(lpSolve); require(lpSolveAPI)					# For linear program solving
require(boot)																	# Basic statistical functions
require(utils)                                                                
require(scatterplot3d)												# For graphics
require(rgl)                                                                  
require(ggplot2)
require(rethinking)
require(proftools)                            #for code profiling
require(profr)
require(compiler)
require(foreach)
require(doMC)
require(optimx)
require(deSolve)
require(plyr) #for loops
require(zoo)
require(pastecs)
require(fractal)
require(tseriesChaos)

registerDoMC()

#` @param model the time series function to use
#' @param initial the initial value for the time series
#' @param netps the number of points desired
#' @parms a list of parameters.  The parameters may be vectors, which must be of equal length.  When this is the case, multiple time series are generated
#' @value a matrix, whose are values, time going down rows, and columns for different model runs.

time.series <- function(model, initial, nsteps, parms) {
  if(length(unique(c(length(initial), sapply(parms, length)))) != 1) {
    stop("Inital values and parameters must be vectors of identical length")
  }
  output = do.call(model, list(initial=initial, nsteps=nsteps, parms=parms))
  return(ts(output))
}

#' @param parms should be a list of vectors of parameter values `a` and `z`
#' @param initial can also be a vector
ricker.series <- function(initial, nsteps, parms) {
  output = matrix(initial,nrow=nsteps,ncol=length(initial), , byrow=TRUE)
  with(parms, {
    for (i in 1:(nsteps-1)) {
      output[i+1,] <- output[i,] * exp(a*(1 - output[i,] + 
        runif(ncol(output), min=-z, max=z)))
    }
    return(output)
  })
}

troph3nd.func <- function(time, states, parms){
  with(as.list(c(states, parms)), {
    fx <- (a1 * x) / (1 + (b1 * x))
    fy <- (a2 * y) / (1 + (b2 * y))
    dx <- (x * (1 - x)) - (fx * y)
    dy <- (fx * y) - (fy * z) - (d1 * y)
    dz <- (fy * z) - (d2 * z)
    return(list(derivs=c(dx,dy,dz)))
  })
}

troph3nd.series <- function(initials=c(x=0.75, y=0.125, z=9.8), nsteps, parms = c(a1 = 5, b1 = 5, a2 = 0.1, b2 = 2, d1 = 0.4, d2 = 0.01)) {
  times <- (1:nsteps)
  output <- lsoda(y=initials, times=times, func=troph3nd.func, parms=parms, rtol = 1e-8, atol = 1e-8)
  return(ts(output[,"z"]))
}

troph3s.func <- function(t, states, parms) {
  with (as.list(c(states, parms)), {
    dR = R * (1 - R/K + forceR(t)*sd) - (xc * yc * C * R)/(R + R0)
    dC = xc * C * (-1 + (yc * R)/(R + R0) + forceC(t)*sd) - (xp * yp * P * C) / (C + C0)
    dP = xp * P * (-1 + (yp * C)/(C + C0) + forceP(t)*sd)
    return(list(derivs=c(dR=dR, dC=dC, dP=dP)))
  })
}

#' Embed a time series in multiple dimensions and eliminate invalid vectrs
#' @param horizons a vector of points in the future to select
#' @import plyr
embed.series <- function(series, dimensions, horizons=1, rem.NA = TRUE, dist = TRUE) {
  z <- series
  times <- matrix(seq(from = tsp(z)[1], to = tsp(z)[2], by = tsp(z)[3]), nrow=length(z))
  dimensions <- dimensions
  vectors <- matrix(NA, nrow = length(z), ncol = dimensions)
  for (i in 1:dimensions) {
    vectors[,i] <- c(rep(NA, i-1), z[1:(length(z)-i+1)])
  }
  futures <- matrix(NA, nrow=length(z), ncol=length(horizons))
  for (i in 1:length(horizons)) {
    futures[,i] = c(z[(horizons[i]+1):length(z)], rep(NA, horizons[i]))
  }
  if(rem.NA) {
    ccases <- complete.cases(cbind(vectors, futures))
    times <- times[ccases,, drop=FALSE]
    vectors <- vectors[ccases, ,drop=FALSE]
    futures <- futures[ccases, ,drop=FALSE]
  }
  colnames(times) <- "Time"
  colnames(vectors) <- paste(rep("D", dimensions), 1:dimensions, sep="")
  colnames(futures) <- paste(rep("H", length(horizons)), horizons, sep="")
  index <- 1:nrow(vectors)
  size <- nrow(vectors)
  if(dist) {
    d <- dist(vectors)
    distmat = as.matrix(d)
    distscale = mean(d)
  } else {
    distmat = NULL
    distscale = NULL
  }
  return(list(vectors=vectors,times=times, index=index, futures=futures, dimensions=dimensions, horizons=horizons, distmat=distmat, distscale=distscale, size=size))
}

#This is a function to pull out individual vectors without fucking up dimensionality and keeping metadata
extract <- function(embedding, rows) {
  e <- embedding
  
  dimensions <- e$dimensions
  horizons <- e$horizons
  size <- length(rows)
  distscale <- ifelse(is.null(e$distscale), NULL, e$distscale)
  
  vectors <- e$vectors[rows, , drop=FALSE]
  futures <- e$futures[rows, , drop=FALSE]

  times <- e$times[rows, , drop=FALSE]
  index <- e$index[rows]

  if(is.null(e$distmat)) {
    distmat <- NULL
  } else {
    distmat <- e$distmat[rows, , drop=FALSE]
  }
  return(list(vectors=vectors,times=times, index=index, futures=futures, dimensions=dimensions, horizons=horizons, distmat=distmat, distscale=distscale, size=size))
}

#This function predicts forward with the simplex method
forecast.simplex <- function(targets, cloud, neighs=(cloud$dimension+1), scale=NULL) {
  if(is.null(scale)) {
    scale <- cloud$distscale
  }

  neighbors <- matrix(NA, nrow=targets$size, ncol=neighs)
  neigh_dist <- neighbors
  futures  <- array(NA, dim=c(targets$size, neighs, length(cloud$horizons)))
  weights <- neighbors
  forecasts <- matrix(NA, nrow=targets$size, ncol=length(cloud$horizons))
  
  for(i in 1:targets$size) {
    neigh_dist[i, ] <- sort(targets$distmat[i,-targets$index])[1:neighs]
    neighbors[i, ] <- match(neigh_dist[i, ], targets$distmat[i, ])
    futures[i, ,] <- cloud$futures[match(neighbors[i, ], cloud$index), ]
    weights[i, ] <- exp(-neigh_dist[i, ]/scale)
    forecasts[i, ] <- apply(futures[i, , , drop=FALSE], 3, function(z) weighted.mean(z, weights[i, ]))
  }

  return(list(neighbors=neighbors, neigh_dist=neigh_dist, futures=futures, weights=weights, forecasts=forecasts))
}


# #cross-validate an embedded series with simplex projection (an older version)
# simplex.cv_old <- function(embedding, allvalues=FALSE) {
#   reals <- matrix(NA, ncol = length(embedding$horizons), nrow = embedding$size )
#   predictions <- reals
#   for(i in 1:embedding$size) { # Do a leave-one-out fit for every point
#     targets <- extract(embedding, i)
#     cloud <- extract(embedding, embedding$index[-i])
#     forecasts <- forecast.simplex(targets, cloud, scale=embedding$distscale) 
#     reals[i, ] <- targets$futures
#     predictions[i, ] <- forecasts$forecasts
#   }
#   correlation <- sapply(1:length(embedding$horizons), function(z) corr(cbind(reals[, z], predictions[, z])))
#   if(allvalues==TRUE) {
#     return(list(reals=reals, predictions=predictions, correlation=correlation))
#   } else {
#     return(correlation)
#   }
# }

#An experimental version that should run faster, but breaks connectivity with the normal forecast.simplex and extract
# Now working!  This is the new version!
simplex.cv <- function(embedding, allvalues = FALSE) {
  predictions <- matrix(NA, ncol = length(embedding$horizons), nrow = embedding$size )
  neighs <- 2:(embedding$dimensions + 2)
  scale <- embedding$distscale
  for(i in 1:embedding$size) {
    neighbors <- order(na.last = TRUE, decreasing = FALSE, embedding$distmat[i, ])[neighs]
    neigh_dist <- embedding$distmat[i, neighs]
    predictions[i, ] <- apply(embedding$futures[neighbors, , drop=FALSE], 2, function(z) weighted.mean(z, exp(-neigh_dist/scale)))
  }
  correlation <- sapply(1:length(embedding$horizons), function(z) corr(cbind(embedding$futures[, z], predictions[, z])))
  if(allvalues==TRUE) {
    return(list(reals=embedding$futures, predictions=predictions, correlation=correlation))
  } else {
    return(correlation)
  }
}

#This function predicts forward with the s-map method
forecast.smap <- function(targets, cloud, theta = 1, scale=NULL) {
  if(is.null(scale)) {
    scale <- cloud$distscale
  }

  weights <- matrix(NA, nrow=targets$size, ncol=cloud$size)
  forecasts <- matrix(NA, nrow = targets$size,ncol = length(cloud$horizons))
  
  for(i in 1:targets$size) {
    weights[i, ] <- exp(-targets$distmat[i, cloud$index] * theta /scale)
    forecasts[i, ] <- apply(cloud$futures, 2, function(z) weighted.mean(z, weights[i, ]))
  }
  
  return(list(weights=weights, forecasts=forecasts))
}

#cross-validate an embedded series with s-map projection

smap.cv_old <- function(embedding, theta, horizon=embedding$horizons, allvalues=FALSE) {
  reals <- matrix(NA, ncol = length(embedding$horizons), nrow = embedding$size )
  predictions <- reals
  for(i in 1:embedding$size) { # Do a leave-one-out fit for every point
    targets <- extract(embedding, i)
    cloud <- extract(embedding, embedding$index[-i])
    forecasts <- forecast.smap(targets, cloud, theta=theta, scale=embedding$scale) 
    reals[i, ] <- targets$futures
    predictions[i, ] <- forecasts$forecasts
  }
  correlation <- sapply(1:length(embedding$horizons), function(z) corr(cbind(reals[, z], predictions[, z])))
  if(allvalues==TRUE) {
    return(list(reals=reals, predictions=predictions, correlation=correlation, theta=theta))
  } else {
  return(correlation[horizon])
  }
}

#experimental, faster, vectorized cross-validation
# Now the main function!
smap.cv <- function(embedding, theta, horizon=embedding$horizons, allvalues=FALSE) {
  weights <- exp(-embedding$distmat * theta/embedding$distscale)
  diag(weights) <- 0
  weights <- weights/rowSums(weights)
  predictions <- weights %*% embedding$futures
  correlation <- sapply(1:length(embedding$horizons), function(z) corr(cbind(embedding$futures[, z], predictions[, z])))
  if(allvalues==TRUE) {
    return(list(reals=embedding$futures, predictions=predictions, correlation=correlation, theta=theta))
  } else {
    return(correlation[horizon])
  }
}

##now fit the smap theta value
fit.smap <- function(embedding, theta=1, horizon=1) {
  minfun <- function(x) smap.cv(embedding=embedding, theta=x, allvalues=FALSE, horizon)
# opt <- optimx(par=c(x=theta), minfun, method="nlm", control=list(trace=1))
# unlist(nlm(minfun, p=25, iterlim=1000))  nlm isn't catching multiple wells!
  unlist(optimize(f=minfun, interval=c(0,5000), maximum=TRUE))
#  out <- optim(par=1, minfun, method="Brent", lower=0, upper=5000, control=list(fnscale=-1))
#  return(c(theta=out$par, correlation=out$value, counts=out$counts[1], convergence=out$convergence))
}

#an attempt to capture all the periodic maxims
# period <- function(x) {
#   a <- spectrum(x, plot=FALSE)
#   period <- 1/a$freq[which.max(a$spec)]
#   r <- rollapply(x, width=round(2*period), max)
#   maxs <- unique(r)
#   maxloc <- match(maxs, x)
#   plot(x); points(maxloc, maxs)
# }

#' Creates a poincare section of a multidimensional timeseries by linear interpolation
#' @param series a multidimensional time series, 
#' @param psvar the column in the matrix of the time series corresponding to the variable of fixed value in the section
#' @param psval the value of psvar at the section
poincaresec <- function(timeseries, psvar, psval, dir="both") {
  dwnpts <- rep(0, nrow(out))
  uppts <- dwnpts
  for (i in 1:(nrow(timeseries)-1)) {
    if((timeseries[i,psvar] > psval) & (timeseries[i+1,psvar] < psval )) dwnpts[i+1] = i+1
  }
  for (i in 1:(nrow(timeseries)-1)) {
    if((timeseries[i,psvar] < psval) & (timeseries[i+1,psvar] > psval )) uppts[i+1] = i+1
  }
  if (dir=="down") pspts <- dwnpts
  if (dir=="up") pspts <- uppts
  if (dir=="both") pspts <- dwnpts + uppts
  pspts=pspts[which(pspts !=0)]
  pcsec <- matrix(NA, nrow=length(pspts), ncol=ncol(timeseries))
  for(i in 1:length(pspts)) {
    pcsec[i, ] <- timeseries[pspts[i]-1, ] + ((timeseries[pspts[i], ] - timeseries[pspts[i]-1,])/(timeseries[pspts[i], psvar ] - timeseries[pspts[i]-1, psvar]))*(psval - timeseries[pspts[i]-1, psvar])
  }
  colnames(pcsec) <- colnames(timeseries)
  return(pcsec)
}

pointsused <- function(embedding, theta, threshold=0.95) {
  weights1 <- exp(-embedding$distmat* theta /embedding$distscale)
  diag(weights1) <- 0
  weights <- aaply(weights1, 1, function(z) z/sum(z))
  mat <- mean(aaply(weights, 1, function(z) which(cumsum(sort(z, decreasing=TRUE)) > threshold, useNames=FALSE)[1] - 1))
  return(mat)
}

firstminmut <- function(x) which(peaks(-mutual(x, plot=FALSE)))[1]
firstacc <- function(x)  which(acf(x, plot=FALSE, type="correlation", lag.max=100)$acf < 0)[1] - 1

FNNdim <- function(x) which(FNN(x, dimension=20, tlag=1)[3,]==0)[1]
