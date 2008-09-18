`sogliapotenzaACCURACY` <-
function(sigma2x0,sigma2y0,rho0,mux0,muy0,sigma2x1,sigma2y1,rho1,mux1,muy1,alpha,n,underH0){
sigmax0    <- sqrt(sigma2x0)
sigmay0    <- sqrt(sigma2y0)
sigmax1    <- sqrt(sigma2x1)
sigmay1    <- sqrt(sigma2y1)
v0         <- abs(muy0 - mux0) / sqrt(sigmay0 * sigmax0)
w0         <- sigmay0 / sigmax0
v1         <- abs(muy1 - mux1) / sqrt(sigmay1 * sigmax1)
w1         <- sigmay1 / sigmax1
Cb0        <- 2 / (w0 + 1 / w0 + v0^2)
Cb1        <- 2 / (w1 + 1 / w1 + v1^2)          
L0         <- log(Cb0/(1 - Cb0))
L1         <- log(Cb1/(1 - Cb1))
prima0     <- Cb0^2 * v0^2 * (w0 + 1 / w0 - 2 * rho0) 
seconda0   <- Cb0^2 * (w0^2 + 1 / w0^2 + 2 * rho0^2) / 2
terza0     <- (1 + rho0^2) * (Cb0 * v0^2 - 1)
argomento0 <- (prima0 + seconda0 + terza0) / ((n - 2) * (1 - Cb0)^2)
sigmaL0    <- sqrt(argomento0)
prima1     <- Cb1^2 * v1^2 * (w1 + 1 / w1 - 2 * rho1) 
seconda1   <- Cb1^2 * (w1^2 + 1 / w1^2 + 2 * rho1^2) / 2
terza1     <- (1 + rho1^2) * (Cb1 * v1^2 - 1)
argomento1 <- (prima1 + seconda1 + terza1) / ((n - 2) * (1 - Cb1)^2)
sigmaL1    <- sqrt(argomento1)
Thr0       <- exp(L0 + qnorm(1 - alpha) * sigmaL0) / (1 + exp(L0 + qnorm(1 - alpha) * sigmaL0))
Thr1       <- exp(L1 + qnorm(1 - alpha) * sigmaL1) / (1 + exp(L1 + qnorm(1 - alpha) * sigmaL1))
Prob       <- 1 - pnorm(((L0 - L1) + qnorm(1 - alpha) * sigmaL0) / sigmaL1)
if(underH0 == TRUE)
list(Stat = Cb0,Prob = 0.05,Thr = Thr0)
else
list(Stat = Cb1,Prob = Prob,Thr = Thr1)
}

