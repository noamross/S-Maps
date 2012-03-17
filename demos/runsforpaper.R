#Runs for ECL232 presentation
#Noam Ross
#This script will reproduce all data and plots for the presentation in ECL232, March 16, 2002
#It requires the 'smap_functions' file to run

#Import packages and functions
source('R/smap_functions.R')  #This will load more libraries that you actually need.

#Create a list of Ricker time series of 500 points, cutting off the initial values
rickerparms <- expand.grid(a=seq(0.5,1,2,3,4), z=c(0.001, 0.005, 0.01, 0.05, 0.1))
initial = 0.8
nsteps=1000
ricks <- alply(rickerparms, 1, function(z) time.series(ricker.series, initial, nsteps, parms=list(a=z[,"a"], z=z[,"z"])), .parallel=TRUE)
ricks <- llply(ricks, function(z) ts(z[980:1000]))
#Plot these outputs
par(mfrow=c(5,5), mar=c(1,0.2,0,0.2), oma=c(2,2,1,1))
l_ply(1:25, function(z) plot(1:500, ricks[[z]], type="l", col="slateblue", ylim=c(0,7),yaxt=if(z%%5 != 1) "n", xaxt=if(z<21) "n"))
par(mfrow=c(5,5), mar=c(1.2,2,1.2,1.2))
a_ply(1:25, 1, function(z) plot(1:500, ricks[[z]], type="l", col="slateblue", yaxt=if(z%%5 != 1) "n", xaxt=if(z<21) "n"))


#Create a list of 3-species model attractors

times=seq(from=0, to=10000, by=1)
parms = c(xc=0.4, yc=2.009, xp=0.08, yp=2.876, R0=0.16129, C0=0.5,  K=0.997,sd=0.001)
inits = c(R=0.5, C=0.4, P=0.8)
Ksvals <- rep(seq(from=0.84, to=0.99, length.out=5), times=5)
stochs <- rep(c(0.0005, 0.001, 0.005, 0.01, 0.05), each=5)
randR = cbind(runif(6002,-1,1))
randC = cbind(runif(6002,-1,1))
randP = cbind(runif(6002,-1,1))
forceR <- cmpfun(approxfun(x=0:6001, y=randR, method="linear", rule=2))
forceC <- cmpfun(approxfun(x=0:6001, y=randR, method="linear", rule=2))
forceP <- cmpfun(approxfun(x=0:6001, y=randR, method="linear", rule=2))
trophs <- alply(1:25, 1, function(z) {
  randR = cbind(runif(6002,-1,1))
  randC = cbind(runif(6002,-1,1))
  randP = cbind(runif(6002,-1,1))
  forceR <- cmpfun(approxfun(x=0:6001, y=randR, method="linear", rule=2))
  forceC <- cmpfun(approxfun(x=0:6001, y=randR, method="linear", rule=2))
  forceP <- cmpfun(approxfun(x=0:6001, y=randR, method="linear", rule=2))
  lsoda(y=c(R=0.5, C=0.4, P=0.8), times=times, func=troph3s.func, parms=replace(parms, c("K", "sd"), c(Ksvals[z], stochs[z])))
}, .parallel=TRUE)

#Plot the 3-species model attractors
par(mfrow=c(5,5), mar=c(0,0,0,0))
citation
l_ply(trophs, function(z) scatterplot3d(z[,2:4], type="l", color=col.alpha("black", alpha=0.3), axis=FALSE, mar=c(0,0,0,0), box=FALSE, xlim=c(0,1), ylim=c(0,1), zlim=c(0,1)))
#Extract independent points and plot
mut.lags <- laply(trophs, function(z) firstminmut(z[,4]))
trophps <- alply(1:25, 1, function(z) ts(trophs[[z]][(101:600)*mut.lags[z],4]))
par(mfrow=c(5,5), mar=c(0.5,0.5,0.5,0.5))
a_ply(1:25,1, function(z) plot(1:500, trophps[[z]], type="l", col="slateblue", ylim=c(0,1), yaxt=if(z%%5 != 1) "n", xaxt=if(z<21) "n"))
par(mfrow=c(5,5), mar=c(0.1,0.1,0.1,0.1), oma=c(0.5,0.5,0.5,0.5))
l_ply(trophps, function(z) plot(1:500, z, type="l", col="slateblue", yaxt="n",xaxt="n")) # Plots with different scales, shows transitions between regimes better
par(mfrow=c(5,5), mar=c(0,0,0,0))
l_ply(trophs, function(z) scatterplot3d(z[(101:600)*mut.lags[z],2:4], type="l", color=col.alpha("black", alpha=1), axis=FALSE, mar=c(0,0,0,0), box=FALSE)) #xlim=c(0,1), ylim=c(0,1), zlim=c(0,1)

#Save the original time series and remove from workspace for memory, along with some other big things
save(trophs, file="trophs.R")
rm(trophs, randC, randP, randR, forceC, forceP, forceR)


#Now let's see how each do with dimensionality
emlist <- expand.grid(E=1:10, ser=1:25)  #make a grid of embeddings in different dimensions
em.t <- alply(emlist,1, function(z) embed.series(series=trophps[[z[,"ser"]]], dimensions=z[,"E"]), .parallel=TRUE) #embed the trophic series
em.r <- alply(emlist, 1, function(z) embed.series(series=ricks[[z[,"ser"]]], dimensions=z[,"E"]), .parallel=TRUE) #embed the ricker series
t.simfits <- llply(em.t, function(z) simplex.cv(z, allvalues=TRUE), .parallel=TRUE)  #Do a simlex projection cv fit to estmate dimensions
r.simfits <- llply(em.r, function(z) simplex.cv(z, allvalues=TRUE), .parallel=TRUE)
t.simfit.table <- laply(em.t, function(z) simplex.cv(z, allvalues=FALSE), .parallel=TRUE)
r.simfit.table <- laply(em.r, function(z) simplex.cv(z, allvalues=FALSE), .parallel=TRUE)
t.simfit.peaks <- aaply(t.simfit.table, 2, function(z) which.max(z))   #Find the maximum fit across dimensions
r.simfit.peaks <- aaply(r.simfit.table, 2, function(z) which.max(z))
par(mfrow=c(5,5), mar=c(1,0.2,0,0.2), oma=c(2,2,1,1))
a_ply(1:20, 1, function(z) {   #Plot the dimensional fits for the trophic series
  plot(1:10, t.simfit.table[,z], type="l", ylim=c(0.8,1), xlim=c(1,10), yaxt=if(z%%5 != 1) "n", xaxt="n", col="slateblue", lwd=4)
  points(t.simfit.peaks[z], t.simfit.table[t.simfit.peaks[z], z], col="blue", pch=16, cex=2)
  })
a_ply(21:25, 1, function(z) {
  plot(1:10, t.simfit.table[,z], type="l", ylim=c(0.8,1), xlim=c(1,10), yaxt=if(z%%5 != 1) "n", col="slateblue", lwd=4)
  points(t.simfit.peaks[z], t.simfit.table[t.simfit.peaks[z], z], col="blue", pch=16, cex=2)
})

a_ply(1:20, 1, function(z) {   #Plot the dimensional fits for the ricker
  plot(1:10, r.simfit.table[,z], type="l", ylim=c(0,1), xlim=c(1,10), yaxt=if(z%%5 != 1) "n", xaxt="n", col="slateblue", lwd=4)
  points(r.simfit.peaks[z], r.simfit.table[r.simfit.peaks[z], z], col="blue", pch=16, cex=2)
})
a_ply(21:25, 1, function(z) {
  plot(1:10, r.simfit.table[,z], type="l", ylim=c(0,1), xlim=c(1,10), yaxt=if(z%%5 != 1) "n", col="slateblue", lwd=4)
  points(r.simfit.peaks[z], r.simfit.table[r.simfit.peaks[z], z], col="blue", pch=16, cex=2)
})

##Plot the projections against the real values for the correct dimensionality for each of these
rick1s <- seq(from=1, to=250, by=10)
troph3s <- seq(from=3, to=250, by=10)
a_ply(rick1s, 1, function(z) plot(r.simfits[[z]]$reals, r.simfits[[z]]$predictions, col=col.alpha("slateblue", 0.5), cex=0.8, pch=16, yaxt=if(((z-1)/10)%%5 != 0) "n", xaxt=if(z<200) "n"))

a_ply(troph3s, 1, function(z) plot(t.simfits[[z]]$reals, t.simfits[[z]]$predictions, col=col.alpha("slateblue", 0.5), cex=0.8, pch=16, yaxt=if(((z-3)/10)%%5 != 0) "n", xaxt=if(z<200) "n"))

#Now let's do S-Map fits
em.r1 <- em.r[rick1s]
em.t3 <- em.t[troph3s]
r.smapfit <- laply(em.r1, fit.smap, .parallel=TRUE)
t.smapfit <- laply(em.t3, fit.smap, .parallel=TRUE)
rtheta.seq = 10^(seq(0,3.7, by=0.1))
rgrd <- expand.grid(ser = 1:25, theta=rtheta.seq)
ttheta.seq = seq(0,150, by=3)
tgrd <- expand.grid(ser = 1:25, theta=ttheta.seq)
rthetacors <- aaply(rgrd, 1, function(z) smap.cv(em.r1[[z[, "ser"]]], theta=z[, "theta"], horizon=1), .parallel=TRUE)
tthetacors <- aaply(tgrd, 1, function(z) smap.cv(em.t3[[z[, "ser"]]], theta=z[, "theta"], horizon=1), .parallel=TRUE)
par(mfrow=c(5,5), mar=c(1,0.2,0,0.2), oma=c(2,2,1,1))
r.smapfit[7,2] <- max(rthetacors[7, ])  #This didn't get fit right by the solver
r.smapfit[7,1] <- rtheta.seq[which.max(rthetacors[7,])]
t.smapfit[25,2] <- max(tthetacors[25, ], na.rm=TRUE)  #This didn't get fit right by the solver
t.smapfit[25,1] <- ttheta.seq[which.max(tthetacors[25,])]
a_ply(1:25, 1, function(z){
  plot(log10(rtheta.seq), rthetacors[z,], type="l",yaxt=if(z%%5 != 1) "n", xaxt=if(z<21) "n", col="slateblue", lwd=4, ylim=c(0,1))
  points(log10(r.smapfit[z,1]), r.smapfit[z,2], , col="blue", pch=16, cex=2)
})

a_ply(1:25, 1, function(z){
  plot(ttheta.seq, tthetacors[z,], type="l",yaxt=if(z%%5 != 1) "n", xaxt=if(z<21) "n", col="slateblue", lwd=4, ylim=c(0.85,1))
  points(t.smapfit[z,1], t.smapfit[z,2], , col="blue", pch=16, cex=2)
})