`sogliapotenzaTDI` <-
function(sigma2x0,sigma2y0,rho0,mux0,muy0,sigma2x1,sigma2y1,rho1,mux1,muy1,alpha = 0.05,n,underH0,la_CP1){
sigmax0    <- sqrt(sigma2x0)
sigmay0    <- sqrt(sigma2y0)
sigmax1    <- sqrt(sigma2x1)
sigmay1    <- sqrt(sigma2y1)
v0         <- abs(muy0 - mux0) / sqrt(sigmay0 * sigmax0)
w0         <- sigmay0 / sigmax0
v1         <- abs(muy1 - mux1) / sqrt(sigmay1 * sigmax1)
w1         <- sigmay1 / sigmax1
epssquare0 <- sigmax0 * sigmay0 * (v0^2 + w0 + 1 / w0 - 2 * rho0)
epssquare1 <- sigmax1 * sigmay1 * (v1^2 + w1 + 1 / w1 - 2 * rho1)   
W0         <- log(epssquare0)
W1         <- log(epssquare1)
gamma0     <- 2 * (1 - v0^4 / (v0^2 + w0 + 1 / w0 - 2 * rho0)^2)
gamma1     <- 2 * (1 - v1^4 / (v1^2 + w1 + 1 / w1 - 2 * rho1)^2)
sigmaW0    <- sqrt(gamma0 / (n - 2))
sigmaW1    <- sqrt(gamma1 / (n - 2))
argomento  <- ((W0 - W1) - qnorm(1 - alpha) * sigmaW0) / sigmaW1
Thr0       <- exp(W0 - qnorm(1 - alpha) * sigmaW0)
Thr1       <- exp(W1 - qnorm(1 - alpha) * sigmaW1)
Prob       <- pnorm(argomento)
nw         <- ((qnorm(Prob) * sqrt(gamma1) + qnorm(1 - alpha) * sqrt(gamma0)) / (W0 - W1))^2 + 2
if(underH0 == TRUE)
list(Stat = qnorm(1 - (1 - la_CP1) / 2) * sqrt(epssquare0),Prob = 0.05,Thr = Thr0)
else
list(Stat = qnorm(1 - (1 - la_CP1) / 2) * sqrt(epssquare1),Prob = Prob,Thr = Thr1)
}

