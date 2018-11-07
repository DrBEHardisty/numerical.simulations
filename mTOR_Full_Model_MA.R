## A mass action kinetics model of the IRS1/PI3K/AKT/mTOR pathway
## from the moment when an Insulin molecule docks at an IRS1 receptor, to the end of the pathway, 
## when cell division occurs because S6 has been activated.

err <- 10^(-10) #accepted error

Insulin <- function(Time, y, pars) {
	with(as.list(c(state, pars)), {
  		dIRS1 <- y[1];
  		dIp <- y[2]; 
  		dPI3k <- y[3];
  		dPIP3 <- y[4];
  		dPIP2 <- y[5];
  		dPTEN <- [y6];
  		dIP <- [y7];
  		dAKT <- [y8];
  		dAKTp <- [y9];
  		dE3 <- [y10];
  		dTSC1 <- [y11];
  		dTSC2 <- [y12];
  		dTSCp <- [y13];
  		dGTPRheb <- [y14];
  		dRAPTOR <- [y15];
  		dmTOR <- [y16];
  		dmTOR_RAPTOR <- [y17];
  		dmTORa <- [y18];
  		dmTORa_AKTp <- [y19];
  		dS6K1 <- [y20];
  		dS6 <- [y21];
  		dMAS <- [y22];
  		dG <- [y23];
  		dAKTp_G <- [y24];
  		dGmem <- [y25];
  		dE1 <- [y26];
  		dE2 <- [y27];
        
  		dIRS1 <- Ki - (K1 * I * IRS1) ## [unbound IRS1]
  		dIp <- (K1 * I * IRS1) + (K2 * I * PKC) + (K19 * AKTp * E1) - (K3 * Ip * PI3K) - (K20 * Ip) ## [phosphorylated IRS1]
  		dPI3K <- K_PI3K + (K4 * IP) - (K5 * PI3K * PIP2) ## [PI3 kinase]
  		dPIP3 <- (K5 * PI3K * PIP2) - (K6 * PIP3 * PTEN) - (K7 * PIP3 * AKT) ## [PIP3, a phospholipid]
  		dPIP2 <- (K6 * PIP3 * PTEN)- (K5 * PI3K * PIP2) ## [PIP2, a phospholipid]
  		dPTEN <- K_PTEN - (K6 * PIP3 * PTEN) ## [enzyme that catalyzes the creation of PIP2]
  		dIP <- (K3 * Ip * PI3K) - (K_IP * IP) ## [phosphorylated IRS1-PI3K complex]
  		dAKT <- K_A + (K8 * AKTp) - (K7 * PIP3 * AKT) ## [available AKT]
  		dAKTp <- (K7 * PIP3 * AKT) + (K14 * mTORa_AKTp) + (K18 * AKTp_G) - (K8 * AKTp) - (K13 * AKTp * mTORa) - (K17 * AKTp * G) - (K19 * AKTp * E1) - (K21 * AKTp * E2) ## [phosphorylated AKT] 
  		dE3 <- K_E3 - (K9 * (TSC1 * TSC2 * E2)) ## ["Enzyme 3": turns TSC1 & TSC2 into TSCp]
  		dTSC1 <- K_T1 - (K9 * (TSC1 * TSC2 * E3)) ## [TSC1]
  		dTSC2 <- K_T2 - (K9 * (TSC1 * TSC2 * E3)) ## [TSC2]
  		dTSCp <- (K9 * (TSC1 * TSC2 * E3)) - (K_T * TSCp) ## [phosphorylated TSC1-TSC2 complex]
  		dGTPRheb <- K_GTP + (K9 * TSC1 * TSC2 * E3) - (K11 * TSCp * GTPRheb * mTOR_RAPTOR) ## VERIFY IF THIS EQN IS CORRECT! 
  		dRAPTOR <- K_R + (K11 * (TSCp * GTPRheb * mTOR_R)) - (K10 * RAPTOR * mTOR) ## [RAPTOR protein]
  		dmTOR <- K_M - (K10 * RAPTOR * mTOR) ## [mTOR protein kinase]
  		dmTOR_RAPTOR <- (K10 * RAPTOR * mTOR) - (K11 * (mTOR_RAPTOR * TSCp * GTPRheb)) ## [mTOR-RAPTOR complex]
  		dmTORa <- (K11 * (GTPRheb * mTOR_RAPTOR * TSCp)) - (K13 * AKTp * mTORa) ## [activated mTOR, or, mTOR-GTPRheb
  		## which forces RAPTOR off of mTOR-RAPTOR and leaves mTOR's active binding sites ready for action]
  		dmTORa_AKTp <- (K13 * AKTp * mTORa) + (K16 * MAS) - (K14 * mTORa_AKTp) ## [mTORa-AKTp complex]
  		dS6K1 <- K_S6 - (K21 * AKTp * E2) - (K15 * mTORa_AKTp * S6K1) - (K22 * S6K1) ## [S6K1 kinase]
  		dS6 <- (K16 * MAS) - (K_S * S6) ## [concentration of S6, a protein that activates cell division & other processes]
  		dMAS <- (K15 * mTORa_AKTp * S6K1) - (K16 * MAS) ## [mTORa-AKTp-S6K1 complex]
  		dG <- K_G - (K18 * AKTp_G) - (K17 * AKTp * G) ## [GLUT4]
  		dAKTp_G <- (K17 *AKTp * G) - (K18 * AKTp_G) ## [AKTp-GLUT4 complex, which moves GLUT4 to membrane to be taken in & used by a cell needing energy in the form of glucose]
  		dGmem <- (K18 * AKTp_G) - (K_G * Gmem) ## [GLUT4 available at cell membrane]
  		dE1 <- K_E1 + (K20 * Ip) - (K19 * AKTp * E1) ## ["Enzyme 1": turns AKTp into Ip]
  		dE2 <- K_E2 + (K22 * S6K1) - (K21 * AKTp * E2)   ## ["Enzyme 2": turns AKTp into S6K1]		

return(list(c(dIRS1, dIp, dPI3K, dPIP3, dPIP2, dPTEN, dIP, dAKT, dAKTp, dE3, dTSC1, dTSC2, dTSCp, dGTPRheb, dRAPTOR, dmTOR, dmTOR_RAPTOR, dmTORa, AKTp, dS6K1, dS6, dMAS, dG, dAKTp_G, dGmem, dE1, dE2)))
	})
}

pars <- c(Ki = 0.1,
		K1 = 0.1,
		K2 = 0.1,
		K3 = 0.1,
		K_PI3K = 0.1,
		K4 = 0.1,
		K5 = 0.1,
		K6 = 0.1,
		K_PTEN = 0.1,
		K_IP = 0.1,
		K7 = 0.1,
		K8 = 0.1,
		K9 = 0.1,
		K10 = 0.1,
		K_A = 0.1,
		K_E3 = 0.1,
		K_T1 = 0.1,
		K_T2 = 0.1,
		K11 = 0.1,
		K_T = 0.1,
		K_GTP = 0.1,
		K_R = 0.1,
		K_M = 0.1,
		K12 = 0.1,
		K13 = 0.1,
		K14 = 0.1,
		K_S6 = 0.1,
		K_S = 0.1,
		K15 = 0.1,
		K16 = 0.1,
		K_G = 0.1,
		K17 = 0.1,
		K18 = 0.1,
		K19 = 0.1,
		K_E1 = 0.1,
		K_E2 = 0.1,
		K20 = 0.1,
		K21 = 0.1)

state <- c(dIRS1=0, dIp=0, dPI3K=0, dPIP3=0, dPIP2=0, dPTEN=0, dIP=0, dAKT=0, dAKTp=0, dE3=0, dTSC1=0, dTSC2=0, dTSCp=0, dGTPRheb=0, dRAPTOR=0, dmTOR=0, dmTOR_RAPTOR=0, dmTORa=0, AKTp=0, dS6K1=0, dS6=0, dMAS=0, dG=0, dAKTp_G=0, dGmem=0, dE1=0, dE2=0) ## a vector of initial values

times <- seq(0, 100, by = 0.2)

library(deSolve)
out <- ode(state, times, Insulin, pars)

## Find the steady states using the Newton-Raphson Method
library(rootSolve)
StS <- stode(y = state, fun = Insulin, parms = pars, pos=TRUE)

## Plot the outputs of the simulation & see what they tell you:
## 2 plots are cool:
plot(out)

## The fancy schmancy colored plot:
## plot_col <- c("black","blue","green","red","purple")
## plot(P ~ time,out, col=plot_col[1], ylab="Protein Kinase Concentration", xlab="Time", main="A Model of the IRS1-mTOR/GLUT4 Pathway")
## lines(Ap ~ time,out,col=plot_col[2],lwd=2)
## lines(Mp ~ time,out,col=plot_col[3],lwd=2)
## lines(Sp ~ time,out,col=plot_col[4],lwd=2)
## lines(G ~ time,out,col=plot_col[5],lwd=2)
## legend("topright",c("[P]","[Ap]","[Mp]","[Sp]","[G]"),col=plot_col,lty=1,lwd=2)

## Matplot also works too
matplot(out[ , 1], out[ , 2:5], type = "l", xlab = "Time", ylab = "[Protein Kinase]", main = "IRS1/PI3K/AKT1/mTOR Pathway", lwd = 2)


## Use a loop to run over parameter values to see what happens:
outlist <- list()
plist <- cbind(Ki = runif(100, min = 0.01, max = 0.99),
				I = runif(100, min = 0.1, max = .99),
			  	a_i = runif(100, min = 0.01, max = 0.99),
			  	sig_P = runif(100, min = 0.01, max = 0.99),
			  	Kap = runif(100, min = 0.01, max = 0.99),
			  	Ta = runif(100, min = 0.01, max = 0.99),
			  	sig_Ap = runif(100, min = 0.01, max = 0.99),
			  	Kam = runif(100, min = 0.01, max = 0.99),
			  	Tm = runif(100, min = 0.01, max = 0.99),
			  	sig_Mp = runif(100, min = 0.01, max = 0.99),
			  	Kms = runif(100, min = 0.01, max = 0.99),
			  	Ts = runif(100, min = 0.01, max = 0.99),
			  	sig_Sp = runif(100, min = 0.01, max = 0.99),
			  	c = runif(100, min = 0.01, max = 0.99),
			  	a_p = runif(100, min = 0.01, max = 0.99),
			  	sig_G = runif(100, min = 0.01, max = 0.99))

for (i in 1:nrow(plist))
	outlist[[i]] <- ode(y = state, parms = plist[i,], func = Insulin, times = times)
plot(out,outlist)

## Find the eigenvalues to determine the stability of my system:
# L <- matrix(c(), nrow=, byrow=TRUE)