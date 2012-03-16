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



ricker.patch.series <- function(initial, nsteps, parms) {
  with(parms, {
  output <- matrix(initial, nrow=nsteps, ncol=patches)
  for (i in 1:(nsteps-1)) {
    z.t <- runif(1, min = -z, max = z)
    if (z.t < 0) {
      d <- matrix(runif(1, min = dl[1], max = dl[2]), nrow=patches, ncol=patches)
    } else {
      d <- matrix(runif(1, min = dh[1], max = dh[2]), nrow=patches, ncol=patches)
    }
    n.tilde <- matrix(output[i, ] * exp(a*(1 - output[i, ] + z.t)))
    output[i+1,] <- d %*% n.tilde
    }
  return(output)
  })
  }

parms = list(a=2.228, z=0.001, dl = c(0.01, 0.02), dh= c(0.02, 0.03))
patches = 8
nsteps = 250
initial = rep(0.2, 8)
runs = 5
r = matrix(NA, nrow=nsteps, ncol=runs)
for (i in 1:runs) {
  r[, i] = rowSums(ricker.patch.series(initial, nsteps, parms))
}
plotseq = seq(from = 1, to = nsteps, by=3)
matplot(plotseq, r[plotseq, ], type="l", lty=1, ylim=c(0,11), pch=3, cex=0.5)

r2 <- ricker.patch.series(a = 2.3, z = 0.001, patches = 8, initial = rep(0.125,8), nsteps = 1000)
matplot(t(r2), type="l")

troph3.series <- function(K, xc, yc, xp, yp, R0, C0)
  
troph3.func <- function(time, states, parms) {
  with(as.list(c(states, parms)), {
    dR = R * (1 - R/K) - (xc * yc * C * R)/(R + R0)
    dC = xc * C * (-1 + (yc * R)/(R + R0)) - (xp * yp * P * C) / (C + C0)
    dP = xp * P * (-1 + (yp * C)/(C + C0))
    return(list(states=c(dR=dR, dC=dC, dP=dP)))
  })
}

fr <- function(u, a, b) {
  f <- 
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
  times <- (1:nsteps)*10
  output <- lsoda(y=initials, times=times, func=troph3nd.func, parms=parms, rtol = 1e-8, atol = 1e-8)
  return(ts(output[,"z"]))
  }
                            
times=seq(from=0, to=5000, by=10)
parmsnd = c(a1 = 5, b1 = 3, a2 = 0.1, b2 = 2, d1 = 0.4, d2 = 0.01)
initsnd = c(x = 1, y = .1, z = 10)
out <- lsoda(y=initsnd, times=times, func=troph3nd.func, parms=parmsnd, rtol = 1e-6, atol = 1e-6)
plot3d(out[,2:4], type="l")

times=seq(from=0, to=5000, by=1)
parms = c(xc=0.4, yc=2.009, xp=0.08, yp=2.876, R0=0.16129, C0=0.5, K=0.997)
inits = c(R=0.5, C=0.4, P=0.8)
out <- lsoda(y=inits, times=times, func=troph3.func, parms=parms, rtol = 1e-8, atol = 1e-8)
#par(mfrow=c(1,1))
#plot(out[,2], type="l")
#plot(out[,3], type="l")
#plot(out[,4], type="l")
#scatterplot3d(out[,2:4], type="l")
scatterplot3d(out[-(4999:5000),4],out[-c(1,5000),4], out[-(1:2),4], type="l")
require(fractaldim)