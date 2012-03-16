# An attempt to replicate plot 1(b) in Hastings and Wysham(2010)
# Noam Ross, March 1, 2012 

# This function generates a time-series based on a multipatch ricker model,
ricker.patch.series <- function(initial, nsteps, parms) {
  
# inital - a vector of initial vectors in each patch
# nsteps - number of steps
# parms - a list of paramers, a and z are stochastic ricker parameters, dl and
#         dh define the low and high ranges of dispersal matrices
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

# set parameters to replicate original paper
parms = list(a=2.228, z=0.001, dl = c(0.01, 0.02), dh= c(0.02, 0.03))
patches = 8
nsteps = 250
#initial values are a guess, but within the range of values produced in output
initial = rep(0.2, 8)
#run the model 5 times and take the total population in each
runs = 5
r = matrix(NA, nrow=nsteps, ncol=runs)
for (i in 1:runs) {
  r[, i] = rowSums(ricker.patch.series(initial, nsteps, parms))
}
#plot, showing every third point
plotseq = seq(from = 1, to = nsteps, by=3)
matplot(plotseq, r[plotseq, ], type="l", lty=1, ylim=c(0,11), pch=3, cex=0.5,
        xlab = "Time (generations)", ylab="Population")