## A simulation of the Lopez-Caamal model of the coupled Akt/mTOR Auto-Activation Pathway
## Simulates behavior of an auto-activation system: the Akt/mTOR protein kinase pathway
## The reactions in the mTOR pathway are:
## P + A >k1> 2A
## A >k2> 0
## P >k3f> <k3b< 0
## A + I >k4> 0
## I >k5f> <k5b< 0

## The differential equations that constitute the model are:
## dC1/dt = -k3fC1-k1C1C2+k3b, C1 is the active form of Akt/mTOR
## dC2/dt = -k2C2+k1C1C2-k4C2C3, C2 is the inactive form of Akt/mTOR
## dC3/dt = -k5fC3-k4C2C3+k5b, C3 is the inhibitor of the auto-activation path

## parameters = k3f, k1, k3b, k2, k4, k5f, k5b
## state variables = C1, C2, C3
err <- 10^(-10) #accepted error
k3f=1; k1=1111; k3b = 1; k2=1; k4 = 1; k5f = 1; k5b = 1;
state <- c(C1 = 0.75, C2= 0, C3= 0.5)

mTOR <- function(Time, y, Parameters) {

  C1 <- y[1];
  C2 <- y[2]; 
  C3 <- y[3];
  dC1 <- -k3f * C1 -k1 * C1 * C2 +k3b ## concentration of Inactive mTOR 
  dC2 <- -k2 * C2 +k1 * C1 * C2 -k4 * C2 * C3 ## concentration of Active mTOR
  dC3 <- -k5f * C3 -k4 * C2 * C3 +k5b ## concentration Inhibited mTOR
return(list(c(dC1, dC2, dC3)))
}
	
times <- seq(0, 100, by = 0.01)

library(deSolve)
out <- ode(y = state, times = times, func = mTOR, parms = NULL, atol = err, rtol = err)

## You can inspect the output for all 3 states if you like
head(out, n=100)
## Tells you a couple of handy things
## diagnostics(out)

## Plot the outputs as 3-D scatter plots to see how parameter values affect levels
## inactive, active and bound Akt/mTOR floating around
## library(scatterplot3d)
## par(mar = c(0, 0, 0, 0))
## scatterplot3d(out[,-1], main="Akt/mTOR Protein Kinase Levels", xlab = "Inactive mTOR", ylab = "Active mTOR", zlab = "Inhibited mTOR", type="p", color="blue", label.tick.marks = TRUE)
matplot(times, out[, 2:4], type = 'p')
