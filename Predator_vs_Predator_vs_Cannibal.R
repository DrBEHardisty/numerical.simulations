library(deSolve)
library(rootSolve)

## There are 2 species, C and S, in this universe.
## C and S compete to eat as much of the resource, Re, as possible.
## They eat Re with efficiency eta. C spends P time eating Re; S spends 1-P time eating Re.
## The dead bodies of C and S decay at rate phi and
## are recycled back into Re at rate Q_S and/or Q_C if they are not eaten.
## F_c is the response of C to the level of Re and D_S, the number of dead S.
## F_s is tte response of S to the level of Re.
## In the beginning, we assume that while species C sometimes eats bodies of C and S,
## S only eats Re.

## The mathematical model
SvsC <- function (Time, state, parms) {
  with(as.list(c(state, parms)), {
	dyRe <- (delta - (eta_C*P*C*Re) - (eta_S*S*Re) - (A_R*Re+ sig_S*Q_S*D_S + sig_C*Q_C*D_C))
	dyC <- C*(phi_S - A_C)
	dyS <- S*(phi_S - A_S)
	dyD_C <- (A_C*C) - (sig_C*D_C)
	dyD_S <- (A_S*S) - (sig_S*D_S) - (eta_D*(1-P)*C*D_S)
    return(list(c(dyRe, dyC, dyS, dyD_C, dyD_S)))
  })
}

n <- 10

## Set the model parameters:
parms <- c(delta = 0.01, Q_S=0.01, Q_C=0.01, A_R=0.01, A_S=0.01, A_C=0.01, P=0.55, eta_C=10, eta_S=0.01, eta_D=0.01, sig_S=0.01, sig_C=0.01, phi_C=0.1, phi_S=0.1)

## Set the initial conditions:
state <- c(Re=5,C=1,S=1,D_C=1,D_S=1)

## How long should you run the simulation for?:
Time <- seq(0, 100, length = n)

## Calculate the solutions, summarize them and plot them:
out <- ode(y, Time, fun = SvsC, parms = parms)
summary(out)
diagnostics(out)
plot(out)

## Use the Newton-Raphson Method to find the steady states,
# given the initial conditions.
## pos = TRUE because only positive steady states are relevant:
StS <- stode(state, fun = SvsC, parms = parms, pos = TRUE)

## Calculate the Jacobian of the model given your initial conditions:
cannibal_jacobian <- jacobian.full(y = c(Re=0,C=0,S=0,D_C=0,D_S=0), func = SvsC , parms = parms)
cannibal_jacobian.2 <- jacobian.full(y = c(Re=0.5,C=0.5,S=0.5,D_C=0.5,D_S=0.5), func = SvsC , parms = parms)
cannibal_jacobian.3 <- jacobian.full(y = c(Re=1,C=1,S=1,D_C=1,D_S=1), func = SvsC , parms = parms)

## Find the eigen values of the Jacobian to evaluate the stability given the Initial Values:
eigen(cannibal_jacobian)$values
eigen(cannibal_jacobian.2)$values
eigen(cannibal_jacobian.3)$values


### Numerical perturbations and their outcomes:
n <- 100
param.name <- "delta" # choose parameter to perturb
param.seq <- seq(0,1,length = 10) # choose range of parameters

# Choose the settings for your bifurcation hunt: 
Pars <- c(delta = 0.01, Q_S=0.01, Q_C=0.01, A_R=0.01, A_S=0.01, A_C=0.01, P=0.55, eta_C=10, eta_S=0.01, eta_D=0.01, sig_S=0.01, sig_C=0.01, phi_C=0.1, phi_S=0.1)
Time <- seq(0, 10, length = n)
State <- c(Re=1,C=1,S=1,D_C=1,D_S=1)
 
param.index <- which(param.name == names(Pars))
out <- list()
for (i in 1:length(param.seq))
  out[[i]] <- matrix(0, n, length(State))
 
for (i in 1:length(param.seq)) {
  # set params
  Pars.loop <- Pars
  Pars.loop[param.index] <- param.seq[i]
  # converge
  init <- ode(State, Time, SvsC, Pars.loop)
  # get converged points
  out[[i]] <- ode(init[n,-1], Time, SvsC, Pars.loop)[,-1]
}
 
range.lim <- lapply(out, function(x) apply(x, 2, range))
range.lim <- apply(do.call("rbind", range.lim), 2, range)
plot.variable <- "Re" # choose which variable to show
plot(0, 0, pch = "", xlab = param.name, ylab = plot.variable,
     xlim = range(param.seq), ylim = range.lim[,plot.variable])
for (i in 1:length(param.seq)) {
  points(rep(param.seq[i], n), out[[i]][,plot.variable])
}
