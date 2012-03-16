
fit.smap.cv <- function(theta = 1, series, dimensions, horizons=1) {
  embedded <- embed.series(series, dimensions=dimensions, horizons=horizons, rem.NA=TRUE)  #turn the time series into an embedding
  cvtargets <- seq(from=1, to=nrow(embedded), by=1)
  sqerrs <- rep(NA, length.out=length(cvtargets))
  predictions <- sqerrs
  for (i in 1:length(cvtargets)) {
      target <- embedded[i,]; attributes(target) <- append(list(dim=c(1,ncol(embedded))), attributes(embedded)[-1])
      cloud <- embedded[-i,]; attributes(cloud) <- append(attributes(cloud), attributes(embedded)[-(1:2)])
      prediction <- forecast.smaps(target, cloud, horizons[1], theta)
      sqerrs[i] <- (target[,1+dimensions+horizons] - prediction)^2
      predictions[i] <- prediction
      }
  return(list(cvtargets=cvtargets, sqerrs=sqerrs, predictions=predictions, embedded=embedded))
  }

smap.cv.ssq <- function(theta = 1, series, dimensions, horizons=1) {
  a <- fit.smap.cv(theta = 1, series, dimensions, horizons=1)
  return(cor(a$embedded[,2], a$predictions))
}
  

fit.smap <- function(series, dimensions, horizons=1) {
  opt <- optimize(f = smap.cv.ssq, interval=c(0,1000), series=series, dimensions=dimensions, horizons=horizons, maximum=TRUE)
  return(opt)
}

fit.smap(series, E, horizon=1)
thetas = seq(from=0, to=50, by=0.25)
errs <- rep(NA, length.out=length(thetas))
for (i in 1:length(thetas)) {
  errs[i] <- fit.smap.cv(theta=thetas[i], series=a, E=1, horizon=1)
}

plot(thetas, errs, type="l")