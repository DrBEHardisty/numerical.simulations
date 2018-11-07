## A version of the BIG model of Topp et al (2000), to model insulin dynamics & mTOR
## B = pancreatic Beta cell mass
## I = fasting insulin concentration in the blood
## G = blood glucose concentration

## Insulin and glucose dynamics happen on a timescale of minutes
## Beta cell dynamics happen on a timescale of days

## The simplest way to study how mTOR could drive the insulin system into dysregulation
## is to alter

## Basal rate of insulin secretion is 43.2
## Ran it at: 80
## Ran it at 120

d_0 = 0.06; r_1 = 0.00084; r_2 = 0.0000024; s = 120; a = 20000; k = 432; R_0 = 864; E_g0 = 1.44; S_i = 0.72

yini <-c(B = 0, I = 0, G = 0)

Ins <- function(t, y, parms) {
	with(as.list(y), {
  		B <- y[1];
  		I <- y[2];
  		G <- y[3];
  		
  		dB <- -(d_0 + r_1 * G - r_2 * G^2) * B
  		dI <- ((B * s * G^2) / (a + B^2)) - (k * I)
		dG <- R_0 - ((E_g0 + S_i * I) * G)
		list(c(dB, dI, dG)) })
}
 
times <- seq(from = 0, to = 50, by = 0.1)

library(deSolve)

## Use the Newton-Raphson Method to find the steady states
library(rootSolve)
StS <- stode(y = state, fun = Ins, parms = pars, pos=TRUE)

out1 <- ode(y = yini, times = times, func = Ins, parms = NULL)
plot(out1)

yini <- c(B = 0.5, I = 0.5, G = 0.5)
out2 <- ode(y = yini, times = times, func = Ins, parms = NULL)
plot(out2)

yini <- c(B = 1, I = 1, G = 1)
out3 <- ode(y = yini, times = times, func = Ins, parms = NULL)
plot(out3)

yini <- c(B = 5, I = 5, G = 5)
out4 <- ode(y = yini, times = times, func = Ins, parms = NULL)
plot(out4)

## Plot all runs to see how things look for different initial conditions, but same parameters
plot(out2, out3)