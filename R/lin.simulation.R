`lin.simulation` <-
function(NUM_CAMP = 5000,NUM = 30,matH0,matH1,underH0 = TRUE,ALPHA_CI = 0.05,la_CP1 = 0.9){
sigma2x    <- sigma2x0 <- matH0[1,1]
sigma2y    <- sigma2y0 <- matH0[2,2]
covxy      <- covxy0   <- matH0[1,2]
mux        <- mux0     <- matH0[1,3]
muy        <- muy0     <- matH0[2,3]

if(underH0 == FALSE){
sigma2x    <- sigma2x1 <- matH1[1,1]
sigma2y    <- sigma2y1 <- matH1[2,2]
covxy      <- covxy1   <- matH1[1,2]
mux        <- mux1     <- matH1[1,3]
muy        <- muy1     <- matH1[2,3]
}

sigmax0    <- sqrt(sigma2x0)
sigmay0    <- sqrt(sigma2y0)
rho0       <- covxy0 / (sigmax0 * sigmay0)
w0         <- sigmay0 / sigmax0
valore     <- sqrt(sigmay0 * sigmax0 * (w0 + 1 / w0 - 2 * rho0))

results    <- .C(NAME = "agreementgenerator",
as.integer(NUM_CAMP),
as.integer(NUM),
lin_z_r   = double(NUM_CAMP),
lin_L_ca  = double(NUM_CAMP),
lin_W     = double(NUM_CAMP),
lin_z_rc  = double(NUM_CAMP),
lin_T1    = double(NUM_CAMP),
lin_T3    = double(NUM_CAMP),
lin_SZ_r  = double(NUM_CAMP),
lin_SL_ca = double(NUM_CAMP),
lin_s_W   = double(NUM_CAMP),
lin_SZ_rc = double(NUM_CAMP),
lin_s_T1  = double(NUM_CAMP),
lin_s_T3  = double(NUM_CAMP),
lin_r_l   = double(NUM_CAMP),
lin_ca_l  = double(NUM_CAMP),
lin_e2_u  = double(NUM_CAMP),
lin_rc_l  = double(NUM_CAMP),
lin_cp_l1 = double(NUM_CAMP),
lin_cp_l3 = double(NUM_CAMP),
alpha     = as.double(ALPHA_CI),
matpar    = as.double(c(sigma2x,covxy,mux,covxy,sigma2y,muy)),
valore    = as.double(valore),
PACKAGE = "agreement"
)

results2  <-  lin.layout(
lin_z_r   = results$lin_z_r,   
lin_L_ca  = results$lin_L_ca,
lin_W     = results$lin_W,
lin_z_rc  = results$lin_z_rc, 
lin_T1    = results$lin_T1,
lin_T3    = results$lin_T3, 
lin_SZ_r  = results$lin_SZ_r,
lin_SL_ca = results$lin_SL_ca,
lin_s_W   = results$lin_s_W,
lin_SZ_rc = results$lin_SZ_rc,
lin_s_T1  = results$lin_s_T1,
lin_s_T3  = results$lin_s_T3, 
lin_r_l   = results$lin_r_l, 
lin_ca_l  = results$lin_ca_l,
lin_e2_u  = results$lin_e2_u,
lin_rc_l  = results$lin_rc_l,
lin_cp_l1 = results$lin_cp_l1,
lin_cp_l3 = results$lin_cp_l3,
la_CP1    = la_CP1,
underH0   = underH0,
matH0     = matH0,
matH1     = matH1,
n         = NUM,
alpha     = ALPHA_CI  
)

list(
table     = results2,
underH0   = underH0,
matH0     = matH0,
matH1     = matH1,
NUM_CAMP  = NUM_CAMP,  
NUM       = NUM,
alpha     = ALPHA_CI,
rho       = covxy / sqrt(sigma2x * sigma2y)
)

}
