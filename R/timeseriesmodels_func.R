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