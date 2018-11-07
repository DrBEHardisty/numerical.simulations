## Johnston et al. (2007) CRC model with saturating feedback
## differentiation rates a_2, b_2 will be multiplied by a small factor, which is
## qualitatively the same effect as activation by mTOR or growth factors
## Additionally, g, the death rate of fully matured cells (N_2) could be
## reduced by a range of fractions to simulate the same effect, that of
## a reduction in policing.

err <- 10^(-10) #accepted error

Colon <- function(Time, y, pars) {
	with(as.list(c(state, pars)), {
  		N_0 <- y[1];
  		N_1 <- y[2];
  		N_2 <- y[3];
  		
  		dN_0 <- (a_3 - a_1 - a_2) * N_0 - (k_0 * N_0^2) / (1 + m_0 * N_0) 
  		dN_1 <- (b_3 - b_1 - b_2) * N_1 - (k_1 * N_1^2) / (1 + m_1 * N_1) + (a_2 * N_0) + (k_0 * N_0^2) / (1 + m_0 * N_0)
		dN_2 <- - (g * N_2) + (b_2 * N_1) + (k_1 * N_1^2) / (1 + m_1 * N_1)
		return(list(c(dN_0, dN_1, dN_2)))
	})
}
 
pars <- c(a_3 = 0.596,
		a_1 = 0.01,
		a_2 = 0.3,
		k_0 = 0.1,
		m_0 = 0.1,
		b_3 = 0.743,
		b_1 = 0.01,
		b_2 = 0.3,
		k_1 = 0.01,
		m_1 = 0.01,
		g = 0.323)

state <- c(N_0=50, N_1=250, N_2=250) ## initial values
times <- seq(0, 100, by = 0.1)

library(deSolve)
out <- ode(state, times, Colon, pars)

## Use the Newton-Raphson Method to find the steady states
library(rootSolve)
StS <- stode(y = state, fun = Colon, parms = pars, pos=TRUE)
 
## Plot outputs & see what they tell you:
plot(out)

## Matplot
matplot(x = out[,1], y = out[,-1], type = "l", lwd = 2,
lty = "solid", col = c("red", "blue", "black"),
xlab = "time", ylab = "State Variables", main = "Colon Crypt Cell Dynamics")
legend("topright", col = c("red", "blue", "black"),
legend = c("Stem cells", "Semi-differientiated Stem Cells", "Mature Cells"), lwd = 2)

