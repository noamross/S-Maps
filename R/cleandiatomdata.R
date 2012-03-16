#'Embed a time series in multiple dimensions and clean the data
#'
#' This function converts the time-series into a matrix of vectors
#' embedded in dimension E whose elements are the data points and their
#' time-lagged points at (t, t-1, ..., t-E+1)

#' This version eliminates vectors that include missing data or project to missing data
#' Missing data must be designated NaN
#' "span" designates whether to omit just vectors with missing data or vectors for which data is missing between the vector and the future value
embed.series2 <- function(series, E, steps = 0, span = FALSE) {
	embedded <- matrix(NA, nrow = (length(series) - E + 1 - steps), ncol = (E+steps))
	for (i in 1:(dim(embedded)[1])) {
		embedded[i, ] = series[i:(i + E + steps - 1)]
		}
		if (span == TRUE) {
			embedded <- as.matrix(embedded[which(!is.nan(rowSums(embedded))), ])
			}
		else {
			embedded <- as.matrix(embedded[which(!is.nan(rowSums(cbind(embedded[, 1:E], embedded[, E+steps])))), ])
			}
		
		training <- list(cloud = as.matrix(embedded[, 1:E]))
		if(steps > 0) training$futures <- embedded[, (E+steps)]
		return(training)
	}
	
	
	

dia <- read.csv("data/SIOdiat.csv")
dit <- dia[which(dia$Year >= 1920),]
dates <- as.Date(paste(dit$Year,dit$Month,dit$Day,sep="-"))
dit <- rowSums(dit[,4:dim(dit)[2]], na.rm=TRUE)
names(dit) <- NULL
#dit[which(dit==0)] = NaN
dit.diffs <- dit[1:(length(dit) - 1)] - dit[-1]
date.steps <- dates[-1] - dates[1:(length(dates)-1)]
dit.diffs[which(date.steps != 7)] <- NA
dit.diffs[which(dit.diffs==0)] <- NA
dit.diffs <- ts(dit.diffs)
dimensions <- 1:10
dems <- alply(1:10, 1, function(z) embed.series(dit.diffs, z))
Rprof("cv.out")
cvdems <- llply(dems, function(z) simplex.cv(z, allvalues=TRUE), .progress="text")
Rprof()
summaryRprof("cv.out")
Rprof("cv2.out")
cvdems3 <- laply(dems, function(z) simplex.cv2(z, allvalues=FALSE), .progress="text")
Rprof()
summaryRprof("cv2.out")

E <- 3
steps <- 1

training <- dit.diffs[1:(round(length(dit.diffs)/2))]
prediction <- dit.diffs[(length(training)+1):(length(dit.diffs))]



train <- embed.series2(training, E, steps = 1, span = TRUE)
cloud <- train$cloud
futures <- train$futures
targs <- embed.series2(prediction, E, steps = 1, span = TRUE)
targets <- targs$cloud
reals <- targs$futures

nhood <- neighborhood(targets,cloud)

neighbors <- nhood$ind
distances <- nhood$dist

forecasts <- forecast(neighbors, distances, futures)


plot(reals, forecasts, pch=0,
		 xlim=c(-1000000,1000000), ylim=c(-1000000,1000000))
		abline(a=0,b=1)

corr(cbind(reals, forecasts))

