a = 1;
h = 1;
eps = .01;
r = 10;
m12 = 10;
m21 = 1;
e = .0001;
q = 1;
k1 = 50;
k2 = 100 ;

initstate <- cbind(100,100,10);
tottime = 100;
timepart = 100; #integer
parms <- cbind();
times <- cbind(0,tottime/timepart*(1:timepart)) 

ratdep <- function(t,y,parms)
{
dy1 <- a*y[1]*y[3]/(1 + a*h*y[1]) + eps*(r*y[1]*(1 - y[1]/k1) + m12*y[2] - m21*y[1])
dy2 <- eps*(r*y[2]*(1 - y[2]/k2) + m21*y[1] - m12*y[2])
dy3 <- (e*a*y[1]*y[3])/(1 + a*h*y[1]) - eps*q*y[3]
list(cbind(dy1,dy2,dy3))
}

out <- ode(y = initstate,times = times,func = ratdep, parms, method = "euler")

matplot(out[,1],cbind(out[,2],out[,3],out[,4]),col = cbind("red","blue","green"),type = "l",xlab = "Time", ylab = "Population",
main = "Ratio Dependent Predation?",lwd = 2)

