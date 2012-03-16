# load up some data
require(yaImpute)																															# For nearest-neighbor search function
require(lpSolve); require(lpSolveAPI)																					# For linear program solving
require(boot)																																	# Basic statistical functions
require(scatterplot3d)
require(rgl)							

	
neighborhood <- function(targets, cloud) {									#predict forward, using simplex data for those that have it and nearest neighbor for those that don't
	dimension <- dim(targets)[2]
	neighbors <- containing.simplexes(targets, cloud)
	edgecases <- which(rowSums(neighbors) == 0)
	edgetargets <- if(is.matrix(targets[edgecases,])) targets[edgecases,] else as.matrix(targets[edgecases,])
	edgeneighbors <- ann(cloud, edgetargets, k=dimension+1, verbose=FALSE )
	neighbors[edgecases,] <- edgeneighbors$knnIndexDist
	nhood <- list(ind=neighbors[,1:(dimension+1)], dist=neighbors[,(dimension+2):(2*dimension+2)])
}
	
forecast <- function(neighbors, distances, futures) {
	weights <- t(sapply(1:(dim(distances)[1]), function(z) exp(distances[z,]/(sum(distances[z,])))))
	weights[which(is.nan(weights[,1])),] <- dim(weights)[2]
	sapply(1:dim(neighbors)[1], function(z) weighted.mean(futures[neighbors[z,]], weights[z,])) #the exponentially weighted mean future step of the neighbors
	}
	
# This function determines whether the target vector is contained in some simplex
# in the cloud of points, a matrix whose rows are vectors.  If it is inside a 
# simplex, it determines the simplex with the smallest average distance from 
# the target to the vertices
containing.simplexes <- function(targets, cloud) {															
		ncloud <- dim(cloud)[1]
		ntargets <- dim(targets)[1]
		distances <- as.matrix(dist(rbind(targets, cloud),diag=TRUE, upper=TRUE ))        								#a vector of distances between the target and the cloud points, TODO: Make this more efficient I might just call it once for all data
		vertices <- t(sapply(1:ntargets, function(z) {
											   containing.simplex(targets[z,], cloud, 
											 	 distances[z, (ntargets+1):(ntargets+ncloud)], ncloud)
													}))
		vertex.distances <- t(sapply(1:ntargets, function(z) distances[z,vertices[1,]+ntargets]))
		vertex.distances[which(rowSums(vertices)==0),] <- 0
		cbind(vertices, vertex.distances)
	}

containing.simplex <- function(target, cloud, distances, m) {
	objective <- distances																												# The objective function is of the form lambda1*distance1 + lambda2*distance2 + ...
	constraint <- rbind(t(cloud), rep(1,length.out=m))														#coefficients for the contraints note that the solver assumes that all the variables are >=0
	constraint.sign <- c(rep("=", length.out = dim(cloud)[2]), "<=")													#equalities/inequalities for the contraints
	constraint.rhs <- c(target, 1)																								#
	out.lp <- lp("min", objective, constraint, constraint.sign, constraint.rhs)
	vertices <- which(out.lp$solution!=0)
	if(length(vertices)!= (ncol(target)+1)) vertices <- rep(0,(dim(cloud)[2]+1))
	return(vertices)
}

predict.smap <- function(data, split, dimensions, steps, theta) {

training <- data[1:split]
prediction <- data[(split+1):(length(data))]
futures <- data[(dimensions+steps):(length(training) + steps)]

cloud <- transform.series(training, dimensions)
targets <- transform.series(prediction, dimensions)


}

forecast.smaps <- function(targets, cloud, horizon, theta) {
	ncloud <- nrow(cloud)
	ntargets <- nrow(targets)
  dimension <- attr(targets, "dimension")
  horizons <- attr(targets, "horizon")
	distances <- as.matrix(dist(rbind(targets, cloud),diag=TRUE, upper=TRUE ))
	distances.cloud <- distances[(ntargets+1):(ntargets+ncloud),(ntargets+1):(ntargets+ncloud)]
	distance.scalar <- mean(distances.cloud[which(lower.tri(distances.cloud))])
	forecast <- sapply(1:ntargets, function(q) {
    target.distances <- distances[q, (ntargets+1):(ntargets+ncloud)]
    weightings <- exp(-theta*(target.distances/distance.scalar))
    forecast <- weighted.mean(cloud[,1+dimension+horizon],weightings)
	  })
  return(forecast)
  }
                      
							
forecast.smap <- function(cloud, theta, distance.scalar, horizon, target.distances) {
	
	
}

predict.simplex <- function(data, split, dimensions, steps) {

	training <- data[1:split]
	prediction <- data[(split+1):(length(data))]
	futures <- data[(dimensions+steps):(length(training) + steps)]
	
	cloud <- transform.series(training, dimensions)
	targets <- transform.series(prediction, dimensions)
	
	nhood <- neighborhood(targets, cloud)
	neighbors <- nhood$ind
	distances <- nhood$dist
	
	forecasts <- forecast(neighbors, distances, futures)
	predicted <- prediction[dimensions:length(prediction)]
	correlation <- corr(cbind(predicted,forecasts))
	mae <- mean(abs(predicted - forecasts))
	simplex <- list(forecasts=forecasts, correlation=correlation, mae=mae, dimensions=dimensions, steps=steps)
}

data <- economics$unemploy

training <- data[1:round(length(economics$unemploy)/2)]								#set the first half as training data
prediction <- data[(length(d.t) + 1):length(economics$unemploy)]					#set the second half as prediction data
dimensions <- 3																																					#set the embedded dimension
steps <- 2																																			#set how far forward to forecast
futures <- data[(dimensions+steps):(length(training) + steps)]																							#define our training data outputs

cloud <- embed.series(training, dimensions)
targets <- embed.series(prediction, dimensions)

nhood <- neighborhood(targets, cloud)
neighbors <- nhood$ind
distances <- nhood$dist

forecasts <- forecast(neighbors, distances, futures)

theta = 0

smap.corr <- function(s) {
	dimensions <- 3																																					#set the embedded dimension
	steps <- 1																																			#set how far forward to forecast
	theta <- s
	
	futures <- data[(dimensions+steps):(length(training) + steps)]																							#define our training data outputs

	cloud <- embed.series(training, dimensions)
	targets <- embed.series(prediction, dimensions)
	sforecast <- forecast.smaps(targets, cloud, futures, theta)
	corr(cbind(prediction[(dimensions + steps):length(prediction)], sforecast[1:(length(sforecast)-steps)]))
}

optim(par=0.5, smap.corr, lower=0, control="trace")
s <- seq(from=1, to=20, by=1)
t <- sapply(s, smap.corr)

plot(prediction[(dimensions + steps):length(prediction)], sforecast[1:(length(sforecast)-steps)])

plot(prediction[-1],forecasts)
corr(cbind(prediction[-1],forecasts))
mean(abs(prediction[-1] - forecasts))

plot(prediction[(dimensions + steps):length(prediction)], sforecast[1:(length(sforecast)-steps)])


steps <- rep(1:50)
corr <- rep(0, length.out=length(steps))
mae <- corr
for(i in 1:length(steps)) {
	a <- predict.simplex(economics$unemploy, 239, 3, i)
	corr[i]=a$correlation
	mae[i]=a$mae
}
plot(steps, corr)
plot(steps, mae)

data <- runif(200)
training <- data[1:100]
prediction <- data[101:200]

### Analyzing diatom data following Sugihara and May (1990), using scripps 1920-39 data

dia <- read.csv("SIOdiat.csv")
dit <- dia[which(dia$Year >= 1920),]
dit <- rowSums(dit[,4:dim(dit)[2]], na.rm=TRUE)
dit <- dit[which(dit!=0)]
dit <- log(dit)
dit <-  dit[1:(length(dit) - 1)] - dit[-1]
names(dit) <- NULL
data <-dit

#data <- economics$unemploy[1:(length(economics$unemploy) -1)] - economics$unemploy[-1]
#require(ggplot2)															
#data(economics)  																															#economic example dataset from ggplot


training <- data[1:round(length(data)/2)]								#set the first half as training data
prediction <- data[(length(training) + 1):length(data)]

step <- 1
dim <- 2
i <- 1
theta <- 0

cs <- rep(0,(length(step)))

for (i in 1:length(dim)) { 
	dimensions <- dim[i]
	steps <- step[i]																															#set the embedded dimension
																																			#set how far forward to forecast
	futures <- data[(dimensions+steps):(length(training) + steps)]

	cloud <- unique(embed.series(training, dimensions))
	targets <- unique(embed.series(prediction, dimensions))

	nhood <- neighborhood(targets, cloud)
	neighbors <- nhood$ind
	distances <- nhood$dist
	
	#TODO: REMOVE DATA WITH IDENTICAL NUMBERS - OR FIND A BETTER WAY TO INCORPORATE THEM. THEY SHOULDN'T FAIL

	forecasts <- forecast(neighbors, distances, futures)

	cs[i] <- corr(cbind(prediction[(dimensions + steps):(length(prediction)-1)], forecasts[1:(length(forecasts)-steps-1)]))
}


plot(forecasts, type="l", ylog=TRUE)
lines(prediction, type="l", col="red")

plot(prediction[(dimensions + steps):length(prediction)], forecasts[1:(length(forecasts)-steps)], pch=0)
		 xlim=c(-1000000,1000000), ylim=c(-1000000,1000000))
		abline(a=0,b=1)
		
		plot(step, cs, type="b", pch=0)# xlim=c(0,7), ylim=c(-0.1,0.5))



theta <- seq(from=0, to=2, by = 0.01)
cs <- rep(0, length(theta))
for (i in 1:(length(theta))) {
	sf <- forecast.smaps(targets,cloud,futures,theta[i])
	cs[i] <- corr(cbind(prediction[(dimensions + steps):(length(prediction)-1)], sf[1:(length(sf)-steps-1)]))
	}

plot(prediction[(dimensions + steps):length(prediction)], sf[1:(length(sf)-steps)], pch=0,
		 xlim=c(-1000000,1000000), ylim=c(-1000000,1000000))
		abline(a=0,b=1)

optim(par=0.5, smap.corr, lower=0, control=list(trace=TRUE,REPORT=1))


#TODO - DVS plots, changing the number of nearest neighbors








