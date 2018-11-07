## A mass action model of the IRS1/AkT/mTOR pathway

err <- 10^(-10) #accepted error

Insulin <- function(Time, y, pars) {
	with(as.list(c(state, pars)), {
  		IRS1P <- y[1];
  		PI3K <- y[2];
  		AkTP <- y[3];
  		TP <- y[4];
  		GTP <- y[5];
  		TGTP <- y[6];
  		mTOR <- y[7];
  		dRAPT <- y[8];
  		mRAPT <- y[9];
  		TmG <- y[10];
  		TmGA <- y[11];
  		GLUT <- y[12];
  		GMEM <- y[13];
  		
  		dIRS1P <- S_IP - (K_IP * IRS1P * PI3K) - (K_IA * IRS1P) - (delta_IRS1P * IRS1P) ## [phosphorylated IRS1]
  		dPI3K <- S_PI3K - (K_IP * IRS1P * PI3K) - (delta_PI3K * PI3K) ## [PI3K protein kinase]
		dAkTP <- S_AkTP + (K_IA * IRS1P) - (K_Am * AkTP * TmG) - (K_AG * AkTP * GLUT) - (delta_AkTP * AkTP)
		## [phosphorylated AkT]
		dTP <- S_TP - (K_TG * TP * GTP) - (delta_TP * TP) ##[phosporylated TSC1-TSC2 heterodimer]
		dGTP <- S_G - (K_TG * TP * GTP) - (delta_GTP * GTP) ##[GTPRheb]
		dTGTP <- (K_TG * TP * GTP) - (K_TGM * TGTP * mRAPT) - (delta_TGTP * TGTP) ## [TSCP-GTPRheb complex]
		dmTOR <- S_mTOR - (K_mR * mTOR * RAPT) - (delta_mTOR * mTOR) ## [mTOR]
		dRAPT <- S_RAPT - (K_mR * mTOR * RAPT) - (delta_RAPT * RAPT) ## [RAPTOR]
		dmRAPT <- (K_mR * mTOR * RAPT) - (K_TGM * TGTP * mRAPT) - (delta_mRAPT * mRAPT) ## [mTOR-RAPTOR complex]
		dTmG <- (K_TGM * TGTP * mRAPT) - (K_Am * AkTP * TmG) - (delta_TmG * TmG) ## [active mTOR]
		dTmGA <- (K_Am * AkTP * TmG) - (delta_TmGA * TmGA) ## [mTOR-AKT^P complex]
		dGLUT <- S_G - (K_AG * AkTP * GLUT) - (delta_GLUT * GLUT) ## [GLUT4 protein]
		dGMEM <- (K_AG * AkTP * GLUT) - (delta_GMEM * GMEM) ## [GLUT4 at membrane] GLUT4 at membrane can then enters cell & is used for growth or movement
return(list(c(dIRS1P, dPI3K, dAkTP, dTP, dGTP, dTGTP, dmTOR, dRAPT, dmRAPT, dTmG, dTmGA, dGLUT, dGMEM)))
	})
}
 
pars <- c(S_IP = 1,
		delta_IRS1P = 0.01,
		K_IP = 0.1,
		K_IA = 0.1,
		S_PI3K = 1,
		delta_PI3K = 0.01,
		S_AkTP = 1,
		K_SA = 0.1,
		K_Am = 0.1,
		K_AG = 0.1,
		delta_AkTP = 0.01,
		S_TP = 1,
		K_TG = 0.1,
		delta_TP = 0.1,
		S_GTP = 1,
		delta_GTP = 0.01,
		delta_TGTP = 0.01,
		S_mTOR = 1,
		delta_mTOR = 0.01,
		S_RAPT = 1,
		delta_RAPT = 0.01,
		K_TGM = 0.1,
		K_mR = 0.1,
		delta_TmG = 0.01,
		delta_TmGA = 0.01,
		delta_mRAPT = 0.01,
		S_G = 1,
		delta_GLUT = 0.01,
		delta_GMEM = 0.01)

state <- c(IRS1P=0, PI3K=0, AkTP=0, TP=0, GTP=0, TGTP=0, mTOR=0, RAPT=10, mRAPT=0, TmG=0, TmGA=0, GLUT=0, GMEM=0) ## initial values
times <- seq(0, 1000, by = 0.1)

library(deSolve)
out <- ode(state, times, Insulin, pars)

## Use the Newton-Raphson Method to find the steady states
library(rootSolve)
StS <- stode(y = state, fun = Insulin, parms = pars, pos=TRUE)

## Use expand.grid to loop over a range of values of Insulin to see how the concentration of Insulin affects state variables
## Use outputs from out.df, a dataframe made of all the outputs of the numerical simulation
out.df <- data.frame(out)
out.sub <- subset(out.df, select = c(IRS1P, PI3K, AkTP, mTOR))
Insul <- seq(0,1,by=0.01)
abmat <- expand.grid(Insul, out.sub)
names(abmat) <- c("IRS1P", "PI3K", "AkTP", "mTOR")
datasave <- abmat
for (i in 1:nrow(abmat)) {
  print(c(i,abmat$IRS1P[i],abmat$PI3K[i],abmat$AkTP[i],abmat$mTOR[i]))
}

2,3,5,8
## Plot the outputs of the simulation & see what they tell you:
plot(out)

## Matplot also works too
matplot(out[ , 1], out[ , 2:14], type = "l", xlab = "Time", ylab = "[Protein Kinase Pathway Element]", main = "IRS1/PI3K/AKT1/mTOR Pathway", lwd = 2)
