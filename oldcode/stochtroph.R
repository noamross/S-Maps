# Stochastic version of 3-species trophic model, with Poincare section sampling.  Note the timestep must be 1
require(deSolve)
require(rgl)
require(scatterplot3d)
#troph3.series <- function(initials, nsteps, parms) {
#  
#}




troph3.func <- function(t, states, parms) {
  with (as.list(c(states, parms)), {
  dR = R * (1 - R/K) - (xc * yc * C * R)/(R + R0)
  dC = xc * C * (-1 + (yc * R)/(R + R0)) - (xp * yp * P * C) / (C + C0)
  dP = xp * P * (-1 + (yp * C)/(C + C0))
  dK = 0
  if(t < 1) {
    dK = 0
  } else if(psup == 1) {
    if((states[psvar] > psval) & (lagvalue(t-1, psvar) < psval)) dK = runif(1, min=-sigma, max=sigma)
  } else if(psup == 0) {
    if((states[psvar] < psval) & (lagvalue(t-1, psvar) > psval)) dK = runif(1, min=-sigma, max=sigma)
  }
  return(list(derivs=c(dR=dR, dC=dC, dP=dP, dK=dK)))
  })
}



times=seq(from=0, to=6000, by=1)
parms = c(xc=0.4, yc=2.009, xp=0.08, yp=2.876, R0=0.16129, C0=0.5, psvar=3, psval=0.95, psup=2, sigma=0.0000)
inits = c(R=0.5, C=0.4, P=0.8, K=0.997)

#create a bifurcation diagram
Kvals <- seq(from=0.84, to=0.99, by=0.001)
trseries2 <- alply(Kvals, 1, function(z) dede(y=c(R=0.5, C=0.4, P=0.8, K=z), times=times, func=troph3.func, parms=parms), .progress="text")
trseries3 <- llply(trseries2, function(z) z[1001:6001, ])
#par(mfrow=c(1,1), mar=c(1,1,1,1))
par(mfrow=c(1,1), mar = c(5, 4, 4, 2) + 0.1)
#l_ply(trseries3, function(z) scatterplot3d(z[,2:4], type="l"))
psecs <- llply(trseries3, function(z) poincaresec(z, psvar=4, psval=0.8, dir="down" ))
#l_ply(psecs, function(z) plot(z[,2], z[, 3], ylim=c(0.15, 0.25), xlim=c(0.65, 0.9)))
plot(0,pch="", xlim=c(min(Kvals), max(Kvals)), ylim=c(0, 1))
for (i in 1:length(Kvals)) {
  y <- poincareMap(trseries3[[i]][,4])$amplitude
  x <- rep(Kvals[i], length(y)) 
  points(x,y, pch=16, cex=0.1)
}
plot(0,pch="", xlim=c(min(Kvals), max(Kvals)), ylim=c(0.1, 0.3), xlab="Carrying Capacity K", ylab="Consumer population at Poincare Section")
for (i in 1:length(Kvals)) {
  y <- psecs[[i]][,3]
  x <- rep(Kvals[i], length(y)) 
  points(x,y, pch=16, cex=0.3, col=col.alpha("blue", alpha=0.5))
}
dims <- laply(trseries3, function(z) FNN(z[,4])["combined",], .progress="text")  #estimated dimensions are 3 for all of them using FNN!
dims2 <- laply(trseries3, function(z) FNN(z[seq(from=1, to=5001, by=10),4])["combined",], .progress="text")
# now estimate dimensions using the simpex projection for a limited set
Ksels <- seq(from=1, to=151, by=15)
Kvals2 <- Kvals[Ksels]
trs <- trseries3[Ksels]
trs <- llply(trs, function(z) ts(z[seq(from=1, to=5001, by=10),4]))
par(mfrow=c(4,3), mar=c(1,1,1,1))
l_ply(trs, plot)
emlist <- expand.grid(E=1:10, ser=1:11)
trem <- alply(emlist, 1, function(z) embed.series(series=trs[[z[,"ser"]]], dimensions=z[,"E"]), .progress="text")
simplex.fits <- laply(trem, simplex.cv, .progress="text")
par(mfrow=c(4,3), mar=c(1,1,1,1))
aaply(simplex.fits,2, which.max)
a_ply(simplex.fits, 2, function(z) plot(1:10, z, type="l", ylim=c(0,1), xlim=c(1,10),col="slateblue"))
#Predictions seem to jump up to 2 and occaisionally 3, but they don't peak.  Beacuse of lack of stochasticity?  No penality for looking at the wrong data, too?
lnames <- as.character(seq(from=3, to=103, by=10))
trem3 <- trem[match(lnames, names(trem))]
thetafits <- laply(trem3, fit.smap, .progress="text")
theta.seq = seq(0.1,50, by=1) 
grd <- expand.grid(ser = 1:11, theta=theta.seq)
thetacors <- aaply(grd, 1, function(z) smap.cv2(trem3[[z[, "ser"]]], theta=z[, "theta"], horizon=1), .progress="text")
a_ply(1:11, 1, function(z) {
  plot(theta.seq, thetacors[z,], type="l")
  points(thetafits[z, "estimate"], -thetafits[z, "minimum"])
})
abline(v=Kvals2)