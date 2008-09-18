`sogliapotenzaCPk3` <-
function(sigma2x0,sigma2y0,rho0,mux0,muy0,sigma2x1,sigma2y1,rho1,mux1,muy1,alpha,n,underH0){
mud0      <- muy0 - mux0
mud1      <- muy1 - mux1
sigmax0   <- sqrt(sigma2x0)
sigmay0   <- sqrt(sigma2y0)
sigmax1   <- sqrt(sigma2x1)
sigmay1   <- sqrt(sigma2y1)
v0        <- abs(muy0 - mux0) / sqrt(sigmay0 * sigmax0)
w0        <- sigmay0 / sigmax0
v1        <- abs(muy1 - mux1) / sqrt(sigmay1 * sigmax1)
w1        <- sigmay1 / sigmax1
sigma2d0  <- sigmax0 * sigmay0 * (w0 + 1 / w0 - 2 * rho0)
sigma2d1  <- sigmax1 * sigmay1 * (w1 + 1 / w1 - 2 * rho1)
k         <- 2.5 * sqrt(sigma2d0)
delta0    <- sign(muy0 - mux0) * v0 / sqrt(w0 + 1 / w0 - 2 * rho0)
delta1    <- sign(muy1 - mux1) * v1 / sqrt(w1 + 1 / w1 - 2 * rho1)
pi0       <- pnorm(k / sqrt(sigma2d0) - delta0) - pnorm(-k / sqrt(sigma2d0) - delta0)
pi1       <- pnorm(k / sqrt(sigma2d1) - delta1) - pnorm(-k / sqrt(sigma2d1) - delta1)
T0        <- log(pi0 / (1 - pi0))
T1        <- log(pi1 / (1 - pi1))
prima0    <- (dnorm(k / sqrt(sigma2d0) - delta0) - dnorm(-k / sqrt(sigma2d0) - delta0))^2
seconda0  <- ((k / sqrt(sigma2d0) - delta0) * dnorm(k / sqrt(sigma2d0) - delta0) + (k / sqrt(sigma2d0) + delta0) * dnorm(k / sqrt(sigma2d0) + delta0))^2
psi0      <- (prima0 + 0.5 * seconda0) / (pi0^2 * (1 - pi0)^2)
prima1    <- (dnorm(k / sqrt(sigma2d1) - delta1) - dnorm(-k / sqrt(sigma2d1) - delta1))^2
seconda1  <- ((k / sqrt(sigma2d1) - delta1) * dnorm(k / sqrt(sigma2d1) - delta1) + (k / sqrt(sigma2d1) + delta1) * dnorm(k / sqrt(sigma2d1) + delta1))^2
psi1      <- (prima1 + 0.5 * seconda1) / (pi1^2 * (1 - pi1)^2)
sigmaT0   <- sqrt(psi0 / (n - 3))
sigmaT1   <- sqrt(psi1 / (n - 3))
argomento <- ((T0 - T1) + qnorm(1 - alpha) * sigmaT0) / sigmaT1
Thr0      <- exp(T0 + qnorm(1 - alpha) * sigmaT0) / (1 + exp(T0 + qnorm(1 - alpha) * sigmaT0))
Thr1      <- exp(T1 + qnorm(1 - alpha) * sigmaT1) / (1 + exp(T1 + qnorm(1 - alpha) * sigmaT1))
Prob      <- 1 - pnorm(argomento)
nt        <- ((qnorm(Prob) * sqrt(psi1) + qnorm(1 - alpha) * sqrt(psi0)) / (T0 - T1))^2 + 3
if(underH0 == TRUE)
list(Stat = pi0,Prob = 0.05,Thr = Thr0)
else
list(Stat = pi1,Prob = Prob,Thr = Thr1)
}

