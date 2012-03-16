x <- rep(1, length.out=100)
y <- x
for (i in 1:1000) {
  x[i+1] = 1/y[i]
  y[i+1] = 
  z[i+1] = y[i]*2
  
  }

plot(x,y, type="l")
cbind(x,y)
# See http://en.wikipedia.org/wiki/List_of_chaotic_maps
a = 1.2
b = 1
n = 500
x0 = pi+0.4
y0 = pi
x = rep(x0, length.out=n)
y = rep(y0, length.out=n)

for (i in 1:n){
     x[i+1] = (x[i] -a* sin(y[i]))%%pi
     y[i+1] = (x[i] + y[i])%%pi
}

plot(x,y,cex=0.1)
head(cbind(x,y))
plot(x,type="l")
plot(x[-1], x[-501], type="l")