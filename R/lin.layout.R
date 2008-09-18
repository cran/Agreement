`lin.layout` <-
function(
lin_z_r,
lin_L_ca,
lin_W,
lin_z_rc,
lin_T1,
lin_T3,
lin_SZ_r,
lin_SL_ca, 
lin_s_W,
lin_SZ_rc,
lin_s_T1,
lin_s_T3,
lin_r_l,
lin_ca_l,
lin_e2_u,  
lin_rc_l,
lin_cp_l1,
lin_cp_l3,
la_CP1,
underH0,
matH0,
matH1,
n,
alpha     
){
sigma2x0 <- matH0[1,1]
sigma2y0 <- matH0[2,2]
rho0     <- matH0[1,2] / sqrt(matH0[1,1] * matH0[2,2])
mux0     <- matH0[1,3]
muy0     <- matH0[2,3]
sigma2x1 <- matH1[1,1]
sigma2y1 <- matH1[2,2]
rho1     <- matH1[1,2] / sqrt(matH1[1,1] * matH1[2,2])
mux1     <- matH1[1,3]
muy1     <- matH1[2,3]
matrice  <- matrix(0,nrow = 6,ncol = 7)
alpha    <- alpha
n        <- n
underH0  <- underH0
rownames(matrice) <- c("Precision","Accuracy","TDI","CCC","Cpk1","CPk3")
colnames(matrice) <- c("Th val","Thr","Th prob","Mean of est","Std of est","Mean of std","Prop rej")
matrice[1,1]      <- sogliapotenzaRHO(sigma2x0,sigma2y0,rho0,mux0,muy0,sigma2x1,sigma2y1,round(rho1,4),mux1,muy1,alpha,n,underH0)$Stat
matrice[2,1]      <- sogliapotenzaACCURACY(sigma2x0,sigma2y0,rho0,mux0,muy0,sigma2x1,sigma2y1,round(rho1,4),mux1,muy1,alpha,n,underH0)$Stat
matrice[3,1]      <- sogliapotenzaTDI(sigma2x0,sigma2y0,rho0,mux0,muy0,sigma2x1,sigma2y1,round(rho1,4),mux1,muy1,alpha,n,underH0,la_CP1)$Stat
matrice[4,1]      <- sogliapotenzaCCC(sigma2x0,sigma2y0,rho0,mux0,muy0,sigma2x1,sigma2y1,round(rho1,4),mux1,muy1,alpha,n,underH0)$Stat
matrice[5,1]      <- sogliapotenzaCPk1(sigma2x0,sigma2y0,rho0,mux0,muy0,sigma2x1,sigma2y1,round(rho1,4),mux1,muy1,alpha,n,underH0)$Stat
matrice[6,1]      <- sogliapotenzaCPk3(sigma2x0,sigma2y0,rho0,mux0,muy0,sigma2x1,sigma2y1,round(rho1,4),mux1,muy1,alpha,n,underH0)$Stat
matrice[1,2]      <- sogliapotenzaRHO(sigma2x0,sigma2y0,rho0,mux0,muy0,sigma2x1,sigma2y1,rho1,mux1,muy1,alpha,n,underH0)$Thr
matrice[2,2]      <- sogliapotenzaACCURACY(sigma2x0,sigma2y0,rho0,mux0,muy0,sigma2x1,sigma2y1,rho1,mux1,muy1,alpha,n,underH0)$Thr
matrice[3,2]      <- sogliapotenzaTDI(sigma2x0,sigma2y0,rho0,mux0,muy0,sigma2x1,sigma2y1,rho1,mux1,muy1,alpha,n,underH0,la_CP1)$Thr
matrice[4,2]      <- sogliapotenzaCCC(sigma2x0,sigma2y0,rho0,mux0,muy0,sigma2x1,sigma2y1,rho1,mux1,muy1,alpha,n,underH0)$Thr
matrice[5,2]      <- sogliapotenzaCPk1(sigma2x0,sigma2y0,rho0,mux0,muy0,sigma2x1,sigma2y1,rho1,mux1,muy1,alpha,n,underH0)$Thr
matrice[6,2]      <- sogliapotenzaCPk3(sigma2x0,sigma2y0,rho0,mux0,muy0,sigma2x1,sigma2y1,rho1,mux1,muy1,alpha,n,underH0)$Thr
if(underH0 == TRUE)
	matrice[1:6,3]    <- 0.05
if(underH0 == FALSE){
	matrice[1,3]      <- sogliapotenzaRHO(sigma2x0,sigma2y0,rho0,mux0,muy0,sigma2x1,sigma2y1,round(rho1,4),mux1,muy1,alpha,n,underH0)$Prob
	matrice[2,3]      <- sogliapotenzaACCURACY(sigma2x0,sigma2y0,rho0,mux0,muy0,sigma2x1,sigma2y1,round(rho1,4),mux1,muy1,alpha,n,underH0)$Prob
	matrice[3,3]      <- sogliapotenzaTDI(sigma2x0,sigma2y0,rho0,mux0,muy0,sigma2x1,sigma2y1,rho1,mux1,muy1,alpha,n,underH0,la_CP1)$Prob
	matrice[4,3]      <- sogliapotenzaCCC(sigma2x0,sigma2y0,rho0,mux0,muy0,sigma2x1,sigma2y1,rho1,mux1,muy1,alpha,n,underH0)$Prob
	matrice[5,3]      <- sogliapotenzaCPk1(sigma2x0,sigma2y0,rho0,mux0,muy0,sigma2x1,sigma2y1,rho1,mux1,muy1,alpha,n,underH0)$Prob
	matrice[6,3]      <- sogliapotenzaCPk3(sigma2x0,sigma2y0,rho0,mux0,muy0,sigma2x1,sigma2y1,rho1,mux1,muy1,alpha,n,underH0)$Prob
}
matrice[1,4]      <- tanh(mean(lin_z_r))
matrice[2,4]      <- exp(mean(lin_L_ca)) / (1 + exp(mean(lin_L_ca)))
matrice[3,4]      <- qnorm(1 - (1 - la_CP1) / 2) * sqrt(exp(mean(lin_W)))
matrice[4,4]      <- tanh(mean(lin_z_rc))
matrice[5,4]      <- exp(mean(lin_T1)) / (1 + exp(mean(lin_T1)))
matrice[6,4]      <- exp(mean(lin_T3)) / (1 + exp(mean(lin_T3)))
matrice[1,5]      <- sd(lin_z_r)
matrice[2,5]      <- sd(lin_L_ca)
matrice[3,5]      <- sd(lin_W)
matrice[4,5]      <- sd(lin_z_rc)
matrice[5,5]      <- sd(lin_T1)
matrice[6,5]      <- sd(lin_T3)
matrice[1,6]      <- mean(lin_SZ_r)
matrice[2,6]      <- mean(lin_SL_ca)
matrice[3,6]      <- mean(lin_s_W)
matrice[4,6]      <- mean(lin_SZ_rc)
matrice[5,6]      <- mean(lin_s_T1)
matrice[6,6]      <- mean(lin_s_T3)
matrice[1,7]      <- 100 * mean(lin_r_l > sogliapotenzaRHO(sigma2x0,sigma2y0,rho0,mux0,muy0,sigma2x1,sigma2y1,rho1,mux1,muy1,alpha,n,underH0)$Thr)
matrice[2,7]      <- 100 * mean(lin_ca_l > sogliapotenzaACCURACY(sigma2x0,sigma2y0,rho0,mux0,muy0,sigma2x1,sigma2y1,rho1,mux1,muy1,alpha,n,underH0)$Thr)
matrice[3,7]      <- 100 * mean(lin_e2_u < sogliapotenzaTDI(sigma2x0,sigma2y0,rho0,mux0,muy0,sigma2x1,sigma2y1,rho1,mux1,muy1,alpha,n,underH0,la_CP1)$Thr)
matrice[4,7]      <- 100 * mean(lin_rc_l > sogliapotenzaCCC(sigma2x0,sigma2y0,rho0,mux0,muy0,sigma2x1,sigma2y1,rho1,mux1,muy1,alpha,n,underH0)$Thr)
matrice[5,7]      <- 100 * mean(lin_cp_l1 > sogliapotenzaCPk1(sigma2x0,sigma2y0,rho0,mux0,muy0,sigma2x1,sigma2y1,rho1,mux1,muy1,alpha,n,underH0)$Thr)
matrice[6,7]      <- 100 * mean(lin_cp_l3 > sogliapotenzaCPk3(sigma2x0,sigma2y0,rho0,mux0,muy0,sigma2x1,sigma2y1,rho1,mux1,muy1,alpha,n,underH0)$Thr)
round(matrice,5)
}
