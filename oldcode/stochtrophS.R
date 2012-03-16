troph3s.func <- function(t, states, parms) {
  with (as.list(c(states, parms)), {
    dR = R * (1 - R/K + forceR(t)*sd) - (xc * yc * C * R)/(R + R0)
    dC = xc * C * (-1 + (yc * R)/(R + R0) + forceC(t)*sd) - (xp * yp * P * C) / (C + C0)
    dP = xp * P * (-1 + (yp * C)/(C + C0) + forceP(t)*sd)
    return(list(derivs=c(dR=dR, dC=dC, dP=dP)))
  })
}

times=seq(from=0, to=10000, by=1)
parms = c(xc=0.4, yc=2.009, xp=0.08, yp=2.876, R0=0.16129, C0=0.5,  K=0.997,sd=0.001)
inits = c(R=0.5, C=0.4, P=0.8)
Ksvals <- rep(seq(from=0.84, to=0.99, length.out=5), times=5)
stochs <- rep(c(0.0005, 0.001, 0.005, 0.01, 0.05), each=5)
tsseries2 <- alply(1:25, 1, function(z) {
  randR = cbind(runif(6002,-1,1))
  randC = cbind(runif(6002,-1,1))
  randP = cbind(runif(6002,-1,1))
  forceR <- cmpfun(approxfun(x=0:6001, y=randR, method="linear", rule=2))
  forceC <- cmpfun(approxfun(x=0:6001, y=randR, method="linear", rule=2))
  forceP <- cmpfun(approxfun(x=0:6001, y=randR, method="linear", rule=2))
  lsoda(y=c(R=0.5, C=0.4, P=0.8), times=times, func=troph3s.func, parms=replace(parms, c("K", "sd"), c(Ksvals[z], stochs[z])))
  }, .parallel=TRUE)
tsseries3 <- llply(tsseries2, function(z) z[1001:6001, ])
par(mfrow=c(5,5), mar=c(0,0,0,0))
l_ply(tsseries3, function(z) scatterplot3d(z[seq(from=1, to=5001, by=10),2:4], type="l", axis=FALSE, mar=c(0,0,0,0), box=FALSE))
#Huzzah!  It has worked!
save(tsseries2, file="3spstochastic.Rdata")
tss <- llply(tsseries3, function(z) ts(z[seq(from=1, to=5001, by=1),4]))
l_ply(tss, plot, type="l")
emlist <- expand.grid(E=1:10, ser=1:25)
Rprof("em.out")
tsem <- alply(emlist, 1, function(z) embed.series(series=tss[[z[,"ser"]]], dimensions=z[,"E"]), .progress="text")
Rprof()
summaryRprof("em.out")
ssimplex.fits <- laply(tsem, simplex.cv, .parallel=TRUE)
par(mfrow=c(5,5), mar=c(0,0,0,0))
a_ply(ssimplex.fits, 2, function(z) plot(1:10, z, type="l", ylim=c(0,1), xlim=c(1,10),col="slateblue"))
aaply(simplex.fits,2, which.max)
dims <- laply(tss, function(z) FNN(z)["combined",], .progress="text")
dims
tsem3 <- tsem[seq(from=3, to=250, by=10)]
sthetafits <- laply(tsem3, fit.smap, .progress="text")
sthetafits
thetafits <- laply(trem3, fit.smap, .progress="text")
stheta.seq = seq(0.1,50, by=1) 
sgrd <- expand.grid(ser = 1:25, theta=theta.seq)
sthetacors <- aaply(sgrd, 1, function(z) smap.cv(tsem3[[z[, "ser"]]], theta=z[, "theta"], horizon=1), .progress="text")
a_ply(1:25, 1, function(z) {
  plot(stheta.seq, sthetacors[z,], type="l")
  points(sthetafits[z, "estimate"], -sthetafits[z,"minimum"])
})


# New question:  What does the theta value correspond to?
plot(laply(tsem3, function(z) z$distscale))

a <- weightmat(tsem3[[1]], 25)
b <- weightmat(tsem[[4]], 55)

# New thing - plot the attractor in embedded space
# Use a fit of theta to color the points according to their weights, or draw a sphere around, a certain projection point.
# Then show these points on the original attractor
# Is this a circle of constant distance?  No.  Should it be?  How should these points be selected?  How about points outisde the attractor?
# Remove autocorrelation using the first minimum of mututal information:
firstminmut <- function(x) which(peaks(-mutual(x, plot=FALSE)))[1]
firstacc <- function(x)  which(acf(x, plot=FALSE, type="correlation", lag.max=100)$acf < 0)[1] - 1
laply(tss, firstminmut)
laply(tss, firstacc)