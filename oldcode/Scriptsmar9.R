#  Initialize my smap functions:
source('~/Dropbox/Workspace/code/smap/R/smap_functions.R')

#Generate a bunch of 3-species time series:
bvals <- seq(from=2, to=5, length.out=25)
trseries <- alply(bvals, 1, function(z) troph3nd.series(nsteps=500, parms = c(a1 = 5, b1 = z, a2 = 0.1, b2 = 2, d1 = 0.4, d2 = 0.01)))
str(trseries)
# plot them
par(mfrow=c(1,1), mar=c(1,1,1,1))
for (i in 1:25) {
  plot(trseries[[i]], col="slateblue", ylim=c(0,15))
  text(400,2, paste("b1 =", bvals[i]))
}
grid <- expand.grid(series=1:25, dimensions=1:10)

#embed the series in 1-10 dimensions, plot the somplex projection fits by dimension
em3 <- alply(grid, 1, function(z) embed.series(series=trseries[[z[,"series"]]], dimensions=z[,"dimensions"], horizons=1:10), .parallel=TRUE)
dim3vc <- llply(em3, function(z) simplex.cv(z, allvalues=TRUE), .parallel=TRUE)
simpfits <- matrix(NA, nrow=25, ncol=10)
simpfits <- simpfits[2:25,]
for (i in 1:250) simpfits[i] <- dim3vc[[i]]$correlation[1]
par(mfrow=c(5,5), mar=c(1,1,1,1))
a_ply(simpfits, 1, function(z) plot(1:10, z, type="l", lty=1, col="slateblue"))
aaply(simpfits, 1, function(x) which(x==max(x)))
#It looks like it picks 3-dimensions most times, but not every time.  Need to return with a lot more series
#Look at projections into the future with 
dim3vc2 <- dim3vc[51:75]
l_ply(dim3vc2, function(z) plot(1:10, z$correlation, type="l", col="slateblue", ylim=c(0.9,1)))
#find optimum theta for each 3-dimensional series
em33 <- em3[51:75]
Rprof("fit2.out")
thetafits2 <- laply(em33, fit.smap, .parallel=TRUE)
Rprof()
summaryRprof("fit2.out")
plot(bvals, thetafits2[,2], type="l", col="slateblue")
thetafits3 <- laply(em33, function(z) fit.smap(z, horizon=10), .parallel=TRUE)
plot(bvals, thetafits3[,2], type="l", col="slateblue")
theta.seq = seq(0.1,5, by=0.1) 
grd <- expand.grid(ser = 1:25, theta=theta.seq)
thetacors <- thetas <- aaply(grd, 1, function(z) smap.cv2(em33[[z[, "ser"]]], theta=z[, "theta"], horizon=10), .progress="text")
par(mfrow=c(5,5), mar=c(1,1,1,1))
a_ply(thetacors, 1, function(z) plot(theta.seq, z, type="l", col="slateblue", ylim=c(0,1)))

thetas <- alply(grd, 1, function(z) smap.cv2(em33[[z[, "ser"]]], theta=z[, "theta"], horizon=1, allvalue=TRUE), .progress="text")

l_ply(thetas[1:25], function(z) plot(z$reals, z$predictions, pch=16, col=col.alpha("slateblue", 0.1)))
l_ply(thetas[626:650], function(z) plot(z$reals, z$predictions, pch=16, col=col.alpha("slateblue", 0.1)))
thetas[650]

scatterplot3d(em33[[25]]$vectors, type="l")


#create a bifurcation diagram
bvals <- seq(from=2.2, to=3.2, length.out=1000)
trseries2 <- alply(bvals, 1, function(z) troph3nd.series(nsteps=1000, parms = c(a1 = 5, b1 = z, a2 = 0.1, b2 = 2, d1 = 0.4, d2 = 0.01)), .parallel=TRUE)
trseries3 <- llply(trseries2, function(z) z[200:1000])
plot(0,pch="", xlim=c(min(bvals), max(bvals)), ylim=c(2, 12.7))
for (i in 1:length(bvals)) {
 y <- poincareMap(trseries3[[i]])$amplitude
 x <- rep(bvals[i], length(y)) 
 points(x,y, pch=16, cex=0.1)
}