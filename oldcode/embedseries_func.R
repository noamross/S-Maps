# Call relevant packages
require(yaImpute)        											# For nearest-neighbor search function
require(lpSolve); require(lpSolveAPI)					# For linear program solving
require(boot)																	# Basic statistical functions
require(utils)                                                                
require(scatterplot3d)												# For graphics
require(rgl)                                                                  
require(ggplot2)
require(plyr)                                 #for loops
require(rethinking)
require(proftools)                            #for code profiling
require(profr)
require(compiler)
require(foreach)
require(doMC)

# A convensionce function for making sure matrices stay matrices
mat2 <- function(x, h=FALSE) {
  if(is.matrix(x)) {
    out <- x
  } else {
    if(h) {
      out <- t(as.matrix(x))
    } else {
      out <- as.matrix(x)
    }
  }
  return(out)
}
mat <- cmpfun(mat2)

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
  futures <- matrix(NA, nrow=length(z), ncol=horizons)
  for (i in 1:length(horizons)) {
    futures[,i] = c(z[(horizons[i]+1):length(z)], rep(NA, horizons[i]))
  }
  if(rem.NA) {
    ccases <- complete.cases(cbind(vectors, futures))
    times <- mat(times[ccases,])
    vectors <- mat(vectors[ccases,])
    futures <- mat(futures[ccases,])
  }
  colnames(times) <- "Time"
  colnames(vectors) <- paste(rep("D", dimensions), 1:dimensions, sep="")
  colnames(futures) <- paste(rep("H", length(horizons)), horizons, sep="")
  index <- matrix(1:nrow(vectors),ncol=1)
  if(dist) {
    d <- dist(vectors)
    distmat = as.matrix(d)
    distscale = mean(d)
  } else {
    distmat = NULL
    distscale = NULL
  }
  size <- nrow(index)
  return(list(vectors=vectors,times=times, index=index, futures=futures, dimensions=dimensions, horizons=horizons, distmat=distmat, distscale=distscale, size=size))
}

#This is a function to pull out individual vectors without fucking up dimensionality and keeping metadata
extract <- function(embedding, rows) {
  e <- embedding
  r1 <- length(rows) == 1
  
  dimensions <- e$dimensions
  horizons <- e$horizons
  distscale <- ifelse(is.null(e$distscale), NULL, e$distscale)
  
  vectors <- mat(e$vectors[rows,], r1)
  futures <- mat(e$futures[rows,], r1)

  times <- mat(e$times[rows,])
  index <- mat(e$index[rows,])

  if(is.null(e$distmat)) {
    distmat <- NULL
  } else {
    distmat <- mat(e$distmat[rows,], r1)
  }
  return(list(vectors=vectors,times=times, index=index, futures=futures, dimensions=dimensions, horizons=horizons, distmat=distmat, distscale=distscale))
}

#This function predicts forward with the simplex method
forecast.simplex <- function(targets, cloud, weighting = 1, scale=FALSE, neighs=(targets$dimension+1)) {
  if(!scale) {
    scale <- cloud$distscale
  }
    mout <- matrix(NA, nrow = nrow(targets$vectors), ncol = 4*neighs + 1 )
    for(z in 1:nrow(targets$vectors)) {
    neigh_dist <- sort(targets$distmat[z,-targets$index ])[1:neighs]
    neighbors <- match(neigh_dist, targets$distmat[z, ])
    futures <- cloud$futures[match(neighbors, cloud$index)]
    weights <- exp(-neigh_dist/scale)
    forecasts <- weighted.mean(futures, weights)
    mout[z, ] <- (c(neighbors, neigh_dist, futures, weights, forecasts))
  }
  neighbors <- mout[, 1:neighs]; storage.mode(neighbors) <- "integer"; colnames(neighbors) <- NULL
  neigh_dist <- mout[, (neighs+1):(2*neighs)]; colnames(neigh_dist) <- NULL
  futures  <- mout[, (2*neighs + 1):(3*neighs)]; colnames(futures) <- NULL
  weights <- mout[, (3*neighs + 1):(4*neighs)]; colnames(weights) <- NULL
  forecasts <- mat(mout[, (4*neighs + 1)]); colnames(forecasts) <- NULL
  return(list(neighbors=neighbors, neigh_dist=neigh_dist, futures=futures, weights=weights, forecasts=forecasts))
}

#cross-validate an embedded series with simplex projection
fit.cv.simplex <- function(embedding) {
  reals <- matrix(NA, ncol = 1, nrow = embedding$size )
  predictions <- reals
  for(i in 1:embedding$size) { # Do a leave-one-out fit for every point
    targets <- extract(embedding, i)
    cloud <- extract(embedding, embedding$index[-i])
    forecasts <- forecast.simplex(targets, cloud) 
    reals[i] <- targets$futures
    predictions[i] <- forecasts$forecasts
    }
    correlation <- corr(cbind(reals, predictions))
    return(correlation)
}

#the same as above, but return all values
fit.cv.simplex2 <- function(embedding) {
  reals <- matrix(NA, ncol = 1, nrow = embedding$size )
  predictions <- reals
  for(i in 1:embedding$size) { # Do a leave-one-out fit for every point
    targets <- extract(embedding, i)
    cloud <- extract(embedding, embedding$index[-i])
    forecasts <- forecast.simplex(targets, cloud) 
    reals[i] <- targets$futures
    predictions[i] <- forecasts$forecasts
  }
  correlation <- corr(cbind(reals, predictions))
  return(list(reals=reals, predictions=predictions, correlation=correlation))
}

##Now a script.  Let's create a matrix of parameters for time series and embedding
rickerparms <- expand.grid(a=c(0.5,1,2,3,4), z=c(0.001, 0.005, 0.01, 0.05, 0.1))
dimensions <- 1:10
initial = 0.8
nsteps=200
ricks <- alply(rickerparms, 1, function(x) time.series(ricker.series, initial, nsteps, parms=list(a=x[,"a"], z=x[,"z"])))
par(mfrow=c(5,5), mar=c(1,1,1,1))
for (i in 1:25) {
  plot(ricks[[i]], col="slateblue")
}
emlist <- expand.grid(E=1:10, ser=1:25)
Rprof("embedding.out")
ems <- alply(emlist, 1, function(x) embed.series(series=ricks[[x[,"ser"]]], dimensions=x[,"E"]))
Rprof(NULL)
summaryRprof("embedding.out")
#Wooo!  Faster! and I've done all the 'dist' calls!
Rprof("correlations.out")
corrs <- laply(ems, fit.cv.simplex, .progress="text")
Rprof(NULL)
summaryRprof("correlations.out")
#Wow!  That damn mat() function is taking all the space.  
par(mfrow=c(5,5))
a_ply(corrs, 2, function(x) plot(1:10, x, type="l", ylim=c(0,1), xlim=c(1,10), col="slateblue"))
#It appears that dimensionality measures are screwed up when you have cyclic behavior.
dims1 <- names(ems)[seq(1,250, by=10)]
emt <- ems[match(dims1, names(ems))]
par(mfrow=c(5,5), mar=c(1,1,1,1))
l_ply(emt, function(x) {
  v <- fit.cv.simplex2(x)
  plot(v$reals, v$predictions, pch=16, cex=0.5, col=col.alpha("slateblue", alpha=0.3), xlim=c(0,5), ylim=c(0,5))
}, .progress="text")
l_ply(emt, function(x) {
  v <- fit.cv.simplex2(x)
  plot(v$reals, v$predictions, pch=16, cex=0.5, col=col.alpha("slateblue", alpha=0.3))
  abline(a=0, b=1, col="grey")
}, .progress="text")

#Test speed.
system.time(corrs <- laply(ems, fit.cv.simplex, .progress="text"))
registerDoMC()
system.time(corrs <- laply(ems, fit.cv.simplex, .progress="text", .parallel=TRUE))
#Hey! doing it in parrallel halves the time!  But don't use the progress bar when doing parallel processing.