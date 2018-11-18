#### James Holland Jones
#### Department of Anthropology
#### Stanford University
#### 7 January 2009
 
# re-create figures in May (1976)
 
# logistic map
lmap <- expression(a*x*(1-x))
 
## make x bound on [0,1] to avoid pathological behavior of negative pop size
x <- seq(0,1,length=1000)
 
### first plot the unstable recruitment curve
a <- 3.414
x1 <- eval(lmap)
plot(x,x1, type="l", xlab=expression(X[t]),ylab=expression(X[t+1]))
 
## equilibrium for logistic map
xstar1 <- 1-(1/a)
 
## numerical derivative
x1p <- diff(x1)
xp <- diff(x)
### the equilibrium x is approximately x[706]
m1 <- x1p[706]/xp[706]
 
### use point-slope eq for a line y - y_1 = m(x - x_1)
### we know the point (xstar,xstar) so solve for eq we can use to draw line
lines(x[550:850], m1*x[550:850]-m1*xstar1+xstar1,lty=2)
 
### now plot the stable recruitment curve
a <- 2.707
x2 <- eval(lmap)
lines(x,x2)
xstar2 <- 1-(1/a)
x2p <- diff(x2)
xp <- diff(x)
### the equilibrium x is approximately x[700]
m2 <- x2p[630]/xp[630]
lines(x[480:780], m2*x[480:780]-m2*xstar2+xstar2,lty=2)
 
lines(x,x)
## the figures are clearly mislabeled in May (1976)
lmap2 <- expression(a*(a*x*(1-x))*(1-(a*x*(1-x))))
 
## this is what figure 2 should be
x2 <- eval(lmap2)
plot(x,x2,type="l",xlab=expression(X[t]),ylab=expression(X[t+2]), ylim=c(0,0.8))
lines(x[480:780], m2^2*x[480:780]-m2^2*xstar2+xstar2,lty=2)
lines(x,x)
 
## this is what figure 3 should be
a <- 3.414
plot(x,eval(lmap2),type="l",xlab=expression(X[t]),ylab=expression(X[t+2]))
lines(x[550:850], m1^2*x[550:850]-m1^2*xstar1+xstar1,lty=2)
lines(x,x)
 
## chaotic region
a <- 3.571
plot(x,eval(lmap2),type="l",xlab=expression(X[t]),ylab=expression(X[t+2]))
lines(x,x)
a <- 3.414
lines(x,eval(lmap2))
 
## human (r=0.02)
a <- 1.64872
xh <- eval(lmap)
plot(x,xh, type="l", xlab=expression(X[t]),ylab=expression(X[t+1]))
lines(x,x)
 
xstarh <- 1-(1/a)
#[1] 0.3934689
xhp <- diff(xh)
 
### search for the point near the equlibrium
max(which(x<0.4))
x[394]
#[1] 0.3933934
mh <- xhp[393]/xp[393]
lines(x[300:500], mh*x[300:500]-mh*xstarh+xstarh,lty=2)
text(0.2,0.4, expression(paste(lambda, "=0.3532", sep="")))
