## A mass action kinetics model of the IRS1/AkT/mTOR pathway
## with both positive and negative feedbacks incorporated

err <- 10^(-10) #accepted error

Insulin <- function(Time, y, pars) {
	with(as.list(c(state, pars)), {
  		IRS1P <- y[1];
  		PI3K <- y[2];
  		AkTP <- y[3];
  		TSC1 <- y[4];
  		TSC2 <- y[5];
  		TSC <- y[6];
  		E3 <- y[7];
  		GTP <- y[8];
  		TGTP <- y[9];
  		mTOR <- y[10];
  		dRAPT <- y[11];
  		mRAPT <- y[12];
  		TmG <- y[13];
  		TmGA <- y[14];
  		S6K1 <- y[15];
  		S6 <- y[16];
  		GLUT4 <- y[17];
  		GMEM <- y[18];
  		E1 <- y[19];
  		E2 <- y[20];
  		
  		dIRS1P <- S_IP + (K_AE1 * AkTP * E1) - (K_IPP * IRS1P * PI3K) - (K_IA * IRS1P) - (delta_IRS1P * IRS1P) ## [phosphorylated IRS1]
  		dPI3K <- S_PI3K - (K_IPP * IRS1P * PI3K) - (delta_PI3K * PI3K) ## [PI3K]
		dAkTP <- S_AkTP + (K_IA * IRS1P) + (K_SA * S6K1) - (K_Am * AkTP * TmG) - (K_AG * AkTP * GLUT4) - (K_AE1 * AkTP * E1) - (K_AE2 * AkTP * E2) - (delta_AkTP * AkTP)
		## [phosphorylated AkT]
		dTSC1 <- S_TSC1 - (K_tsc * TSC1 * TSC2) - (delta_TSC1 * TSC1)
		dTSC2 <- S_TSC2 - (K_tsc * TSC1 * TSC2) - (delta_TSC2 * TSC2)
		dTSC <- S_TSC - (K_TG * TSC * GTP) - (delta_TSC * TSC) ## [S_TSCP = K_TE3 * TSC * E3][phosporylated TSC1-TSC2 heterodimer]
		dE3 <- S_E3 - (K_te * TSC * E3) - (delta_E3 * E3)
		dGTP <-  (K_E3 * TSC * E3)- (K_TG * TSC * GTP) - (delta_GTP * GTP) ## [GTPRheb]
		dTGTP <- (K_TG * TSC * GTP) - (K_TGM * TGTP * mRAPT) - (delta_TGTP * TGTP) ## [TSCP-GTPRheb complex]
		dmTOR <- S_mTOR - (K_mR * mTOR * RAPT) - (delta_mTOR * mTOR) ## [mTOR]
		dRAPT <- (K_TGM * TGTP * mRAPT) - (K_mR * mTOR * RAPT) - (delta_RAPT * RAPT) ## [RAPTOR]
		dmRAPT <- (K_mR * mTOR * RAPT) - (K_TGM * TGTP * mRAPT) - (delta_mRAPT * mRAPT) ## [mTOR-RAPTOR complex]
		dTmG <- (K_TGM * TGTP * mRAPT) - (K_Am * AkTP * TmG) - (delta_TmG * TmG) ## [active mTOR]
		dTmGA <- (K_Am * AkTP * TmG) - (K_SAm * TmGA * S6K1) - (delta_TmGA * TmGA) ## [mTOR-AKT complex]
		dS6K1 <- S_S6K1 + (K_AE2 * AkTP * E2) - (K_SAm * TmGA * S6K1) - (K_SA * S6K1) - (delta_S6K1 * S6K1) ## [S6K1 protein kinase]
		dS6 <- (K_SAm * TmGA * S6K1) - (delta_S6 * S6) ## [S6 translation initiation factor]
		dGLUT4 <- S_G - (K_AG * AkTP * GLUT4) - (delta_GLUT4 * GLUT4)
		## [GLUT4 protein] GLUT4 is permanently lost when AkTP-GLUT4 complex arrives at membrane
		dGMEM <- (K_AG * AkTP * GLUT4) - (delta_GMEM * GMEM) ## [GLUT4 at membrane] GLUT4 at membrane can then enter the cell and be used for growth or movement
		dE1 <- S_E1 + (K_IA * IRS1P) - (K_AE1 * AkTP * E1) - (delta_E1 * E1)## [enzyme E1]
		dE2 <- S_E2 + (K_SA * S6K1) - (K_AE2 * AkTP * E2) - (delta_E2 * E2) ## [enzyme E2]
return(list(c(dIRS1P, dPI3K, dAkTP, dTSC1, dTSC2, dTSC, dE3, dGTP, dTGTP, dmTOR, dRAPT, dmRAPT, dTmG, dTmGA, dS6K1, dS6, dGLUT4, dGMEM, dE1, dE2)))
	})
}
 
pars <- c(S_IP = 1,
		I = 20.0,
		delta_IRS1P = 0.01,
		K_IPP = 0.1,
		K_IA = 0.1,
		S_PI3K = 2,
		K_IPP = 0.1,
		delta_PI3K = 0.01,
		S_AkTP = 2,
		K_IA = 0.1,
		K_SA = 0.1,
		K_Am = 0.1,
		K_AG = 0.1,
		K_AE1 = 0.1,
		K_AE2 = 0.1,
		delta_AkTP = 0.1,
		S_TSC = 1,
		K_tsc = 0.1,
		S_TSC1 = 1,
		delta_TSC1 = 0.1,
		S_TSC2 = 1,
		delta_TSC2 = 0.1,
		K_TG = 0.1,
		delta_TSC = 0.1,
		K_te = 0.1,
		S_E3 = 5,
		delta_E3 = 0.1,
		K_E3 = 0.1,
		S_GTP = 1,
		delta_GTP = 0.1,
		delta_TGTP = 0.1,
		S_mTOR = 10,
		delta_mTOR = 0.1,
		delta_RAPT = 0.01,
		K_TGM = 0.1,
		K_mR = 0.1,
		delta_TmG = 0.1,
		delta_TmGA = 0.1,
		delta_mRAPT = 0.1,
		S_S6K1 = 1,
		K_SAm = 0.1,
		delta_S6K1 = 0.1,
		delta_S6 = 0.1,
		S_G = 1,
		delta_GLUT4 = 0.1,
		delta_GMEM = 0.1,
		S_E1 = 1,
		delta_E1 = 0.1,
		S_E2 = 1,
		delta_E2 = 0.1)

state <- c(IRS1P=0, PI3K=0, AkTP=0, TSC1=0, TSC2=0, TSC=0, E3=0, GTP=0, TGTP=0, mTOR=0, RAPT=0, mRAPT=0, TmG=0, TmGA=0, S6K1=0, S6=0, GLUT4=0, GMEM=0, E1=0, E2=0) ## initial values
times <- seq(0, 1000, by = 0.1)

library(deSolve)
out <- ode(state, times, Insulin, pars)

## Use the Newton-Raphson Method to find the steady states
library(rootSolve)
StS <- stode(y = state, fun = Insulin, parms = pars, pos=TRUE)
 
## Plot the outputs of the simulation & see what they tell you:
plot(out)

## Matplot also works too
matplot(out[ , 1], out[ , 2:18], type = "l", xlab = "Time", ylab = "[Protein Kinase Pathway Element]", main = "IRS1/PI3K/AKT1/mTOR Pathway", lwd = 2)

##Run parameters as a loop to see what they look like:
outlist <- plist()
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

for (i in 1:nrow(list))
	outlist[[i]] <- ode(y = state, parms = plist[i,], func = Insulin, times = times)
plot(out,outlist)