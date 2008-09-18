`sogliapotenzaRHO` <-
function(sigma2x0,sigma2y0,rho0,mux0,muy0,sigma2x1,sigma2y1,rho1,mux1,muy1 = 0.1,alpha,n,underH0){
Z0        <- 0.5 * log((1 + rho0) / (1 - rho0))
Z1        <- 0.5 * log((1 + rho1) / (1 - rho1))
sigmaZ0   <- sqrt(1 / (n - 2))
sigmaZ1   <- sqrt(1 / (n - 2))
argomento <- ((Z0 - Z1) + qnorm(1 - alpha) * sigmaZ0) / sigmaZ1
Thr0      <- tanh(Z0 + qnorm(1 - alpha) * sigmaZ0)
Thr1      <- tanh(Z1 + qnorm(1 - alpha) * sigmaZ1)
Prob      <- 1 - pnorm(argomento)
if(underH0 == TRUE)
list(Stat = rho0,Prob = 0.05,Thr = Thr0)
else
list(Stat = rho1,Prob = Prob,Thr = Thr1)
}

