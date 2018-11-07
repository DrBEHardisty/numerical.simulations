## "Insulin" is a simplified model of the PI3K-AkT-GLUT4 pathway
## Only features PI3K, AkT & GLUT4
## Inputs: concentration of insulin, reaction rates, decay rates
## Outputs: concentrations of PI3K, AkT and GLUT4

err <- 10^(-10) #accepted error
k1 = 0.1; ## rate of input of Insulin into the pathway
I = 5; ## rate of infusion of Insulin into the PI3K pathway
sig_P = 0.01; ## rate at which P is lost to binding with AkT

k2 = 0.1; ## rate of input of AkT into the system
K_m = 0.1; ## half-saturation constant of binding between IRS1-PI3K complex & AkT
sig_A = 0.01; ## rate at which phosphorylated AkT is lost
A_p = 0.25; ## rate at which phosphorylated AP is lost 

c = 0.1; ## rate
a_p = 0.10; ## rate of infusion of phosphorylated AkT into GLUT4
sig_G = 0.10; ## rate at which GLUT4 is lost

state <- c(P=10, A=10, G=20)

Insulin <- function(Time, y, Parameters) {

  P <- y[1];
  A <- y[2]; 
  G <- y[3];
  dP <- k1 * I - sig_P * P ## [concentration of PI3K]
  dA <- ((k2 * P * A) / (K_m + A)) - sig_A * A_p ## [concentration of AkT]
  dG <- c * a_p - sig_G * G ## [concentration of GLUT4]
return(list(c(dP, dA, dG)))
}
	
times <- seq(0, 100, by = 0.1)

library(deSolve)
out <- ode(y = state, times = times, func = Insulin, parms = NULL, atol = err, rtol = err)
## RS <- runsteady(y = state, fun = Insulin, parms = pars, times = c(0, 1e5))
library(rootSolve)
stode(y = state, fun = Insulin, parms = pars, pos = TRUE)

## You can inspect the output for all 3 states if you like
## head(out, n=1000)
## 'Diagnostics' (sometimes) tells you a couple of handy things
diagnostics(out)

## Plot the outputs to see how parameter values affect levels of PI3K, AkT & GLUT4
## matplot(times, out[, 2:4], type = 'l', ylab="P13K-AkT-GLUT4")
matplot(times, out[,2], ylab = "Insulin", xlab="IRS1-PI3K Levels", main = "A Simple Model of the IRS1 Pathway",type = 'l')
## library(scatterplot3d)
## scatterplot3d(out[,-1], main="Protein Levels", xlab = "PI3K", ylab = "AkT", zlab = "GLUT4", type="p", color="blue", label.tick.marks=TRUE)