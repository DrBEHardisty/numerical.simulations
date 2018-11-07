## model of mTOR and insulin activity, S = input of amino acids to the insulin/mTOR coupled pathway

err <- 10^(-10) ## accepted error

mTORPA <- function(Time, y, pars) {
	with(as.list(c(state, pars)), {
  		mTOR <- y[1];
  		insul <- y[2];

  		dmTOR <- (K_1 * S) - (K_2 * insul * mTOR) ##
  		dinsul <- (K_3 * S) - (K_4 * insul) ##	
	
return(list(c(dmTOR, dinsul)))
	})
}

pars <- c(K_1 = 0.01,
		K_2 = 0.001,
		K_3 = 0.1,
		S = 0.2,
		K_4 = 0.001)

state <- c(mTOR=1, insul=1) ## initial values
times <- seq(0, 2500, by = 0.1)

library(deSolve)
out <- ode(state, times, mTORPA, pars)

## Use the Newton-Raphson Method to find the steady states
library(rootSolve)
StS <- stode(y = state, fun = mTORPA, parms = pars, pos=TRUE)

plot(out)

matplot(out[ , 1], out[ , 2:3], type = "l", xlab = "Time", ylab = "[Protein Kinase Pathway Element]", main = "PAS: Insulin & mTOR Adapt to Levels of Amino Acids", lwd = 2)

