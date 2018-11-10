## Does sensitivity analysis/numerical perturbation
## of whatever parameters one chooses
library(deSolve)
library(rootSolve)

## There are 2 species, C and S, in this universe.
## C and S compete to eat as much of the resource, Re, as possible.
## They eat Re with efficiency eta. C spends P time eating Re; S spends 1-P time eating Re.
## The dead bodies of C and S decay at rate phi and
## are recycled back into Re at rate Q_S and/or Q_C if they are not eaten.
## F_c is the response of C to the level of Re and D_S, the number of dead S.
## F_s is hte response of S to the level of Re.
## In the beginning, we assume that while species C sometimes eats bodies of C and S,
## S only eats Re.

## The mathematical model
SvsC <- function (Time, y, parms) {
  with(as.list(c(state, parms)), {
	Re <-y[1];
	C <- y[2];
	S <- y[3];
	D_C <- y[4];
	D_S <- y[5];
	dyRe <- (delta - (eta_C*P*C*Re) - (eta_S*S*Re) - (A_R*Re+ sig_S*Q_S*D_S + sig_C*Q_C*D_C))
	dyC <- C*(phi_S - A_C)
	dyS <- S*(phi_S - A_S)
	dyD_C <- (A_C*C) - (sig_C*D_C)
	dyD_S <- (A_S*S) - (sig_S*D_S) - (eta_D*(1-P)*C*D_S)
    return(list(c(dyRe, dyC, dyS, dyD_C, dyD_S)))
  })
}

n <- 10

## Set the model parameters and the simulation initial conditions:
parms <- c(Q_S=0.01, Q_C=0.01, A_R=0.01, A_S=0.01, A_C=0.01, P=0.25, delta=10, eta_C=10, eta_S=10, eta_D=20, sig_S=0.1, sig_C=0.1, phi_C=0.5, phi_S=0.5)

y <- c(Re=1,C=1,S=1,D_C=1,D_S=1)

## How long should you run the simulation for?:
Time <- seq(0, 10, length = n)

## Calculate the solutions, summarize them and plot them:
out <- ode(y, Time, SvsC, parms)
summary(out)
plot(out)

## Use the Newton-Raphson Method to find the steady states,
# given the initial conditions:
StS <- stode(y = y, fun = SvsC, parms = parms, pos = TRUE)

