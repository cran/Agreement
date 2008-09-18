`sogliapotenzaCCC` <-
function(sigma2x0,sigma2y0,rho0,mux0,muy0,sigma2x1,sigma2y1,rho1,mux1,muy1,alpha,n,underH0){
sigmax0   <- sqrt(sigma2x0)
sigmay0   <- sqrt(sigma2y0)
sigmax1   <- sqrt(sigma2x1)
sigmay1   <- sqrt(sigma2y1)
v0        <- abs(muy0 - mux0) / sqrt(sigmay0 * sigmax0)
w0        <- sigmay0 / sigmax0
v1        <- abs(muy1 - mux1) / sqrt(sigmay1 * sigmax1)
w1        <- sigmay1 / sigmax1
Cb0       <- 2 / (w0 + 1 / w0 + v0^2)
rhoc0     <- rho0 * Cb0
Cb1       <- 2 / (w1 + 1 / w1 + v1^2)
rhoc1     <- rho1 * Cb1
Z0        <- 0.5 * log((1 + rhoc0) / (1 - rhoc0))
Z1        <- 0.5 * log((1 + rhoc1) / (1 - rhoc1))
prima0    <- (1 - rho0^2) * rhoc0^2 / ((1 - rhoc0^2) * rho0^2)
seconda0  <- 2 * rhoc0^3 * (1 - rhoc0)* v0^2 / (rho0 * (1 - rhoc0^2)^2)
terza0    <- rhoc0^4 * v0^4 / (2 * rho0^2 * (1 - rhoc0^2)^2)
eta0      <- prima0 + seconda0 - terza0
prima1    <- (1 - rho1^2) * rhoc1^2 / ((1 - rhoc1^2) * rho1^2)
seconda1  <- 2 * rhoc1^3 * (1 - rhoc1)* v1^2 / (rho1 * (1 - rhoc1^2)^2)
terza1    <- rhoc1^4 * v1^4 / (2 * rho1^2 * (1 - rhoc1^2)^2)
eta1      <- prima1 + seconda1 - terza1
sigmaZ0   <- sqrt(eta0 / (n - 2))
sigmaZ1   <- sqrt(eta1 / (n - 2))
argomento <- ((Z0 - Z1) + qnorm(1 - alpha) * sigmaZ0) / sigmaZ1
Thr0      <- tanh(Z0 + qnorm(1 - alpha) * sigmaZ0)
Thr1      <- tanh(Z1 + qnorm(1 - alpha) * sigmaZ1)
Prob      <- 1 - pnorm(argomento)
nz        <- ((qnorm(Prob) * sqrt(eta1) + qnorm(1 - alpha) * sqrt(eta0)) / (Z0 - Z1))^2 + 2
if(underH0 == TRUE)
list(Stat = rhoc0,Prob = 0.05,Thr = Thr0)
else
list(Stat = rhoc1,Prob = Prob,Thr = Thr1)
}

