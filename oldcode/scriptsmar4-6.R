# Scripts run Week of Mar7

##Now a script.  Let's create a matrix of parameters for time series and embedding
rickerparms <- expand.grid(a=c(0.5,1,2,3,4), z=c(0.001, 0.005, 0.01, 0.05, 0.1))
dimensions <- 1:10
initial = 0.8
nsteps=5000
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