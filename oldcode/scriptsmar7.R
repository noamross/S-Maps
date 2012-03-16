rickerparms <- expand.grid(a=c(0.5,1,2,3,4), z=c(0.001, 0.005, 0.01, 0.05, 0.1))
dimensions <- 1:10
initial = 0.8
nsteps=200
ricks <- alply(rickerparms, 1, function(x) time.series(ricker.series, initial, nsteps, parms=list(a=x[,"a"], z=x[,"z"])))
emlist <- expand.grid(E=1:10, ser=1:25)
ems <- alply(emlist, 1, function(x) embed.series(series=ricks[[x[,"ser"]]], dimensions=x[,"E"]))
targets <- extract(ems[[12]], 1:50)
cloud <- extract(ems[[12]], 51:198)
Rprof("mat.out")
corrs <- laply(ems, fit.cv.simplex, .progress="text")
Rprof(NULL)
summaryRprof("mat.out")
par(mfrow=c(5,5))
a_ply(corrs, 2, function(x) plot(1:10, x, type="l", ylim=c(0,1), xlim=c(1,10), col="slateblue"))
a <- forecast.smap(targets,cloud, theta=200)
a$forecasts
par(mfrow=c(1,1))
plot(targets$futures, a$forecasts)

#testing the smap fitting function
rickerparms <- expand.grid(a=c(0.5,1,2,3,4), z=c(0.001, 0.005, 0.01, 0.05, 0.1))
dimensions <- 1:10
initial = 0.8
nsteps=200
ricks <- alply(rickerparms, 1, function(x) time.series(ricker.series, initial, nsteps, parms=list(a=x[,"a"], z=x[,"z"])))
emt <- llply(ricks, function(x) embed.series(x, dimensions=1))
thetas = seq(1,50,by=1)
emlistt = expand.grid(theta = thetas, series = 1:25)
Rprof("a.out")
smapfits <- aaply(emlistt, 1, function(x) smap.cv(embedding=emt[[x[,"series"]]], theta=x[,"theta"]))
Rprof()
summaryRprof("a.out")
par(mfrow=c(5,5), mar=c(1,1,1,1))
a_ply(smapfits, 2, function(z) plot(x=thetas,y=z, type="l", xlim=c(0,50), ylim=c(0,1), col="slateblue"))
smapfits2 <- alply(emlistt, 1, function(x) fit.cv.smap(embedding=emt[[x[,"series"]]], theta=x[,"theta"], allvalues=TRUE), .parallel=TRUE)
maxs <- aaply(smapfits,2, which.max)
for (i in 1:25) {
  maxfit <- maxs[i] + (i-1)*50
  plot(smapfits2[[maxfit]]$reals, smapfits2[[maxfit]]$predictions, pch=16, col=col.alpha("slateblue", alpha=0.3))
}

# It looks like the same thing happens with S-maps that happens with simplexes, the max for theta is still often very large.


#Now time to incorporate more horizons.  A test
series <- time.series(ricker.series, initial, nsteps, parms=list(a=3, z=0.05))
mh <- embed.series(series, dimensions=1, horizons=c(1,3,4))
mh$futures
targets <- extract(mh, 1)
cloud <- extract(mh, 5:195)
v <- forecast.simplex(targets, cloud)
forecast.smap(targets, cloud, theta=1)
fit.cv.smap(mh, theta=1)
fit.cv.simplex(mh)

#It works! Let's see how simplex projection works across the ricker set.
rickerparms <- expand.grid(a=c(0.5,1,2,3,4), z=c(0.001, 0.005, 0.01, 0.05, 0.1))
dimensions <- 1:10
initial = 0.8
nsteps=200
ricks <- alply(rickerparms, 1, function(x) time.series(ricker.series, initial, nsteps, parms=list(a=x[,"a"], z=x[,"z"])))
par(mfrow=c(5,5), mar=c(1,1,1,1))
l_ply(ricks, function(z) plot(1:200, z, col="slateblue", type="l"))
emt <- llply(ricks, function(z) embed.series(z, dimensions=1, horizons=1:10))
fits <- laply(emt, function(z) fit.cv.simplex(z), .progress="text")
a_ply(fits, 1, function(z) plot(1:10,z,col="slateblue", type="l", ylim=c(0,1) ) )
fits2 <- llply(emt, function(z) fit.cv.simplex(z, allvalues=TRUE), .progress="text")
for (i in 1:10) l_ply(fits2, function(z) plot(z$reals[,i], z$predictions[,i], col=col.alpha("slateblue",alpha=0.3), pch=16, cex=0.5))

#Look at one closely

targets <- extract(embedding, 1:5)
cloud <- extract(embedding, 6:190)

#Now time to fit optimal thetas
emt <- llply(ricks, function(z) embed.series(z, dimensions=1, horizons=1))
embedding <- emt[[3]]
fit.smap(embedding)
Rprof("fit.out")
fits <- laply(emt, fit.smap)
Rprof()
summaryRprof("fit.out")
par(mfrow=c(5,5), mar=c(1,1,1,1))
for (i in 1:25) {
  plot(x=thetas,y=smapfits[, i], type="l", xlim=c(0,max(fits[,"estimate"])), ylim=c(0,1), col="slateblue")
  points(fits[i, "estimate"], -fits[i, "minimum"], col="red", pch=16)
}