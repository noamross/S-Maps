# Walkthough of update Smap functions

## Load relevant libraries

require(yaImpute)  														# For nearest-neighbor search function
require(lpSolve); require(lpSolveAPI)					# For linear program solving
require(boot)																	# Basic statistical functions
require(utils)                                                                
require(scatterplot3d)												# For graphics
require(rgl)                                                                  
require(ggplot2)

# Always test single-dimensional and multi-dimensional cases!

walk <- cumsum(rnorm(n=500, mean=0)) # generate a random walk


dit <- read.csv("data/SIOdiat.csv")
#dit <- dia[which(dia$Year >= 1920),]
dates <- as.Date(paste(dit$Year,dit$Month,dit$Day,sep="-"))
dit <- rowSums(dit[,4:dim(dit)[2]], na.rm=TRUE)
names(dit) <- NULL
dit[which(dit==0)] = NA
dit.diffs <- dit[1:(length(dit) - 1)] - dit[-1]
date.steps <- dates[-1] - dates[1:(length(dates)-1)]
dit.diffs[which(date.steps != 7)] <- NA

data <- dit.diffs

#plot(data, type="l")

data <- a
dimensions <- 1:10
corri = rep(0, length.out=10)
for (i in 1:length(corri)) {
E <- 2
steps <-1

cloud <- embed.series(data[1:(length(data)/2)], E, steps=steps, cleanup=0)$cloud
futures <- embed.series(data[1:(length(data)/2)], E, steps=steps, cleanup=0)$futures
targets <- embed.series(data[(length(data)/2 + 1):(length(data))], E, steps=steps, cleanup=0)$cloud
reals <- embed.series(data[(length(data)/2 + 1):(length(data))], E, steps=steps, cleanup=0)$futures
              
#xlim <- c(max(walk), min(walk))
#ylim <- xlim
#plot(cloud, xlim=xlim, ylim=ylim)
#points(targets, col="red")
neighbors <- neighborhood(targets, cloud, method="simplex")

predictions <- forecast(neighbors$indices, neighbors$distances, futures)

#predictions <- forecast.smaps(targets, cloud, futures, theta=theta[i])

plot(reals, predictions$forecasts,  pch=0)
#abline(a=0, b=1)
corri[i] <- corr(cbind(reals,predictions$forecasts))
}