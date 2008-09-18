#include <Rmath.h>
#include <R.h>

double calcolamedia(int n,double *vettore){
  int j = 1;
  double result = 0;
  for(j = 1;j < n + 1;j++)
  result = result + vettore[j - 1] / n;
  return(result);
}

double calcolavarianza(int n,double *vettore){
  int j = 1;
  double result = 0;
  for(j = 1;j < n + 1;j++)
  result = result + vettore[j - 1] * vettore[j - 1];
  result = result / n - calcolamedia(n,vettore) * calcolamedia(n,vettore);
  return(result);
}

double calcolacovarianza(int n,double *vettore1,double *vettore2){
  int j = 1;
  double result = 0;
  for(j = 1;j < n + 1;j++)
  result = result + vettore1[j - 1] * vettore2[j - 1];
  result = result / n - calcolamedia(n,vettore1) * calcolamedia(n,vettore2);
  return(result);
}

double calcolacorr(int n,double *vettore1,double *vettore2){
  double result = 0;
  result = calcolacovarianza(n,vettore1,vettore2) / sqrt(calcolavarianza(n,vettore1) * calcolavarianza(n,vettore2));
  return(result);
}

double calcolaaccuracy(int n,double *vettore1,double *vettore2){
  double u = 0;
  double v = 0;
  double result = 0;
  u = (calcolamedia(n,vettore2) - calcolamedia(n,vettore1)) / sqrt(sqrt(calcolavarianza(n,vettore1)) * sqrt(calcolavarianza(n,vettore2)));
  v = sqrt(calcolavarianza(n,vettore2) / calcolavarianza(n,vettore1));
  result = 2 / (u * u + v + 1 / v);
  return(result);
}

double calcolaccc(int n,double *vettore1,double *vettore2){
  double result = 0;
  result = calcolaaccuracy(n,vettore1,vettore2) * calcolacorr(n,vettore1,vettore2);
  return(result);
}

double calcolaMSD(int n,double *vettore1,double *vettore2){
  int j = 1;
  double result = 0;
  for(j = 1;j < n + 1;j++)
  result = result + (vettore1[j - 1] - vettore2[j - 1]) * (vettore1[j - 1] - vettore2[j - 1]);
  result = result / (n - 1);
  return(result);
}

double calcolavarianzadifferenza(int n,double *vettore1,double *vettore2){
  double result = 0;
  result = (calcolavarianza(n,vettore1) + calcolavarianza(n,vettore2) - 2 * calcolacovarianza(n,vettore1,vettore2)) * n / (n - 3);
  return(result);
}

double calcoladifferenza(int n,double *vettore1,double *vettore2){
  double result = 0;
  result = (calcolamedia(n,vettore2) - calcolamedia(n,vettore1)) * (calcolamedia(n,vettore2) - calcolamedia(n,vettore1)) / calcolavarianzadifferenza(n,vettore1,vettore2);
  return(result);
}

double calcolaprobchisqnc(double quantile,double df,double nc){
  double result = 0;
  result = pnchisq(quantile,df,nc,1,0);
  return(result);
}

double calcolaquantchisqnc(double probab,double df,double nc){
  double result = 0;
  result = qnchisq(probab,df,nc,1,0);
  return(result);
}

double calcolaCP(int n,double *vettore1,double *vettore2,double k,double valore){
  double result = 0;
  double linkp = (1 + 0.5 * k) * valore;
  result = ((linkp / calcolavarianzadifferenza(n,vettore1,vettore2) * linkp) > calcolaquantchisqnc(0.9999,1,calcoladifferenza(n,vettore1,vettore2))) ? calcolaquantchisqnc(0.9999,1,calcoladifferenza(n,vettore1,vettore2)) : (linkp / calcolavarianzadifferenza(n,vettore1,vettore2) * linkp);
  result = calcolaprobchisqnc(result,1,calcoladifferenza(n,vettore1,vettore2));
  return(result);
}

double calcolaSZr(int n){
  double result = 0;
  result = 1 / sqrt(n - 3);
  return(result);
}

double calcolaSZrc(int n,double *vettore1,double *vettore2){
  double result = 0;
  double PEARSON_R = 0;
  double lin_rc = 0;
  double lin_U = 0;
  PEARSON_R = calcolacorr(n,vettore1,vettore2);
  lin_rc = calcolaccc(n,vettore1,vettore2);
  lin_U = (calcolamedia(n,vettore2) - calcolamedia(n,vettore1)) / sqrt(sqrt(calcolavarianza(n,vettore1)) * sqrt(calcolavarianza(n,vettore2)));
  result = sqrt(((1 - PEARSON_R * PEARSON_R) * lin_rc * lin_rc / ((1 - lin_rc * lin_rc) * PEARSON_R * PEARSON_R) + 2 * lin_U * lin_U * (1 - lin_rc) * lin_rc * lin_rc * lin_rc / ((1 - lin_rc * lin_rc) * (1 - lin_rc * lin_rc) * PEARSON_R) - lin_U * lin_U * lin_U * lin_U * lin_rc * lin_rc * lin_rc * lin_rc / (2 * (1 - lin_rc * lin_rc) * (1 - lin_rc * lin_rc) * PEARSON_R * PEARSON_R)) / (n - 2));
  return(result);
}

double calcoladifferenzam(int n,double *vettore1,double *vettore2){
  double result = 0;
  result = (calcolamedia(n,vettore2) - calcolamedia(n,vettore1)) * (calcolamedia(n,vettore2) - calcolamedia(n,vettore1));
  return(result);
}

double calcolaSW(int n,double *vettore1,double *vettore2){
  double result = 0;
  double lin_d2 = 0;
  double lin_e2_UE = 0;
  lin_d2 = calcoladifferenzam(n,vettore1,vettore2);
  lin_e2_UE = calcolaMSD(n,vettore1,vettore2);
  result = sqrt((2 * (1 - lin_d2 * lin_d2 / (lin_e2_UE * lin_e2_UE))) / (n - 2));
  return(result);
}

double calcolaSLca(int n,double *vettore1,double *vettore2){
  double result = 0;
  double lin_ca = 0;
  double lin_U = 0;
  double lin_V = 0;
  double lin_r = 0;
  lin_ca = calcolaaccuracy(n,vettore1,vettore2);
  lin_U = (calcolamedia(n,vettore2) - calcolamedia(n,vettore1)) / sqrt(sqrt(calcolavarianza(n,vettore1)) * sqrt(calcolavarianza(n,vettore2)));
  lin_V = sqrt(calcolavarianza(n,vettore2) / calcolavarianza(n,vettore1));
  lin_r = calcolacorr(n,vettore1,vettore2);
  result = sqrt(((lin_ca * lin_ca * lin_U * lin_U * (lin_V + 1 / lin_V - 2 * lin_r) + lin_ca * lin_ca * (lin_V * lin_V + 1 / (lin_V * lin_V) + 2 * lin_r * lin_r) / 2 + (1 + lin_r * lin_r) * (lin_ca * lin_U * lin_U - 1)) / ((n - 2) * (1 - lin_ca) * (1 - lin_ca))));
  return(result);
}

double calcolaST(int n,double *vettore1,double *vettore2,double k,double lin_cp,double valore){
  double result = 0;
  double kk = (1 + 0.5 * k) * valore;
  double lin_kpm = 0;
  double lin_kmm = 0;
  lin_kpm = (kk + (calcolamedia(n,vettore2) - calcolamedia(n,vettore1))) / sqrt(calcolavarianzadifferenza(n,vettore1,vettore2));
  lin_kmm = (kk - (calcolamedia(n,vettore2) - calcolamedia(n,vettore1))) / sqrt(calcolavarianzadifferenza(n,vettore1,vettore2));
  result = sqrt(((dnorm(-lin_kpm,0,1,0) - dnorm(lin_kmm,0,1,0)) * (dnorm(-lin_kpm,0,1,0) - dnorm(lin_kmm,0,1,0)) + (lin_kmm * dnorm(lin_kmm,0,1,0) + lin_kpm * dnorm(-lin_kpm,0,1,0)) * (lin_kmm * dnorm(lin_kmm,0,1,0) + lin_kpm * dnorm(-lin_kpm,0,1,0)) / 2) / ((n - 3) * lin_cp * lin_cp * (1 - lin_cp) * (1 - lin_cp)));
  return(result);
}

double calcolaquantile(double alpha){
  double result = 0;
  result = qnorm(1 - alpha,0,1,1,0);
  return(result);
}

void agreementgenerator(int *NUM_CAMP,int *NUM,double *result1,double *result2,double *result3,double *result4,double *result5,double *result6,double *result7,double *result8,double *result9,double *result10,double *result11,double *result12,double *result13,double *result14,double *result15,double *result16,double *result17,double *result18,double *alpha,double matpar[2][3],double *valore){
int VOLTE = 1;
int nr = *NUM;
int nc = 2;
int i = 0;
int j = 0;
double zmatrice[nr][nc];
double nmatrice[nr][nc];
double vettore1[nr];
double vettore2[nr];
double matchol[2][2];

/* matrice di Choleskyi */

matchol[0][0] = sqrt(matpar[0][0]);
matchol[1][0] = 0;
matchol[0][1] = sqrt(matpar[1][1]) * matpar[0][1] / sqrt(matpar[0][0] * matpar[1][1]);  
matchol[1][1] = sqrt(matpar[1][1]) * sqrt(1 - (matpar[0][1] / sqrt(matpar[0][0] * matpar[1][1])) * (matpar[0][1] / sqrt(matpar[0][0] * matpar[1][1])));

for(VOLTE = 0;VOLTE < *NUM_CAMP;VOLTE++){

/* generazione dei numeri casuali  da una normale standard */

GetRNGstate();
for(j = 0;j < 2;j++)
for(i = 0;i < nr;i++)
zmatrice[i][j] = rnorm(0,1);
PutRNGstate();

/* generazione dei numeri casuali  da una normale qualunque */

for(i = 0;i < nr;i++){
for(j = 0;j < 2;j++){
nmatrice[i][j] = zmatrice[i][0] * matchol[0][j] + zmatrice[i][1] * matchol[1][j] + matpar[j][2];
}
}

for(i = 0;i < nr;i++){
vettore1[i] = nmatrice[i][0];
vettore2[i] = nmatrice[i][1];
}

result1[VOLTE]  = atanh(calcolacorr(*NUM,&vettore1[0],&vettore2[0]));
result2[VOLTE]  = calcolaaccuracy(*NUM,&vettore1[0],&vettore2[0]);
result2[VOLTE]  = log(result2[VOLTE] / (1 - result2[VOLTE]));
result3[VOLTE]  = log(calcolaMSD(*NUM,&vettore1[0],&vettore2[0]));
result4[VOLTE]  = atanh(calcolaccc(*NUM,&vettore1[0],&vettore2[0]));
result5[VOLTE]  = calcolaCP(*NUM,&vettore1[0],&vettore2[0],1,*valore);
result11[VOLTE] = calcolaST(*NUM,&vettore1[0],&vettore2[0],1,result5[VOLTE],*valore);
result5[VOLTE]  = log(result5[VOLTE] / (1 - result5[VOLTE]));
result6[VOLTE]  = calcolaCP(*NUM,&vettore1[0],&vettore2[0],3,*valore);
result12[VOLTE] = calcolaST(*NUM,&vettore1[0],&vettore2[0],3,result6[VOLTE],*valore);
result6[VOLTE]  = log(result6[VOLTE] / (1 - result6[VOLTE]));
result7[VOLTE]  = calcolaSZr(*NUM);
result8[VOLTE]  = calcolaSLca(*NUM,&vettore1[0],&vettore2[0]);
result9[VOLTE]  = calcolaSW(*NUM,&vettore1[0],&vettore2[0]);
result10[VOLTE] = calcolaSZrc(*NUM,&vettore1[0],&vettore2[0]);
result13[VOLTE] = tanh(result1[VOLTE] - calcolaquantile(*alpha) * result7[VOLTE]);
result14[VOLTE] = result2[VOLTE] - calcolaquantile(*alpha) * result8[VOLTE];
result14[VOLTE] = exp(result14[VOLTE]) / (1 + exp(result14[VOLTE]));
result15[VOLTE] = exp(result3[VOLTE] + calcolaquantile(*alpha) * result9[VOLTE]);
result16[VOLTE] = tanh(result4[VOLTE] - calcolaquantile(*alpha) * result10[VOLTE]);
result17[VOLTE] = result5[VOLTE] - calcolaquantile(*alpha) * result11[VOLTE];
result17[VOLTE] = exp(result17[VOLTE]) / (1 + exp(result17[VOLTE]));
result18[VOLTE] = result6[VOLTE] - calcolaquantile(*alpha) * result12[VOLTE];
result18[VOLTE] = exp(result18[VOLTE]) / (1 + exp(result18[VOLTE]));
}
}

void agreementgenerator2(int *NUM_CAMP,int *NUM,double *result1,double *result2,double *result3,double *result4,double *result5,double *result6,double *result7,double *result8,double *result9,double *result10,double *result11,double *result12,double *alpha,double matpar[2][3],double *valore){
int VOLTE = 1;
int nr = *NUM;
int nc = 2;
int i = 0;
int j = 0;
double zmatrice[nr][nc];
double nmatrice[nr][nc];
double vettore1[nr];
double vettore2[nr];
double matchol[2][2];

/* matrice di Choleskyi */

matchol[0][0] = sqrt(matpar[0][0]);
matchol[1][0] = 0;
matchol[0][1] = sqrt(matpar[1][1]) * matpar[0][1] / sqrt(matpar[0][0] * matpar[1][1]);  
matchol[1][1] = sqrt(matpar[1][1]) * sqrt(1 - (matpar[0][1] / sqrt(matpar[0][0] * matpar[1][1])) * (matpar[0][1] / sqrt(matpar[0][0] * matpar[1][1])));

for(VOLTE = 0;VOLTE < *NUM_CAMP;VOLTE++){

/* generazione dei numeri casuali  da una normale standard */

GetRNGstate();
for(j = 0;j < 2;j++)
for(i = 0;i < nr;i++)
zmatrice[i][j] = rnorm(0,1);
PutRNGstate();

/* generazione dei numeri casuali  da una normale qualunque */

for(i = 0;i < nr;i++){
for(j = 0;j < 2;j++){
nmatrice[i][j] = zmatrice[i][0] * matchol[0][j] + zmatrice[i][1] * matchol[1][j] + matpar[j][2];
}
}

for(i = 0;i < nr;i++){
vettore1[i] = nmatrice[i][0];
vettore2[i] = nmatrice[i][1];
}

result1[VOLTE]  = atanh(calcolacorr(*NUM,&vettore1[0],&vettore2[0]));
result2[VOLTE]  = calcolaaccuracy(*NUM,&vettore1[0],&vettore2[0]);
result2[VOLTE]  = log(result2[VOLTE] / (1 - result2[VOLTE]));
result3[VOLTE]  = log(calcolaMSD(*NUM,&vettore1[0],&vettore2[0]));
result4[VOLTE]  = atanh(calcolaccc(*NUM,&vettore1[0],&vettore2[0]));
result5[VOLTE]  = calcolaCP(*NUM,&vettore1[0],&vettore2[0],1,*valore);
result11[VOLTE] = calcolaST(*NUM,&vettore1[0],&vettore2[0],1,result5[VOLTE],*valore);
result5[VOLTE]  = log(result5[VOLTE] / (1 - result5[VOLTE]));
result6[VOLTE]  = calcolaCP(*NUM,&vettore1[0],&vettore2[0],3,*valore);
result12[VOLTE] = calcolaST(*NUM,&vettore1[0],&vettore2[0],3,result6[VOLTE],*valore);
result6[VOLTE]  = log(result6[VOLTE] / (1 - result6[VOLTE]));
result7[VOLTE]  = calcolaSZr(*NUM);
result8[VOLTE]  = calcolaSLca(*NUM,&vettore1[0],&vettore2[0]);
result9[VOLTE]  = calcolaSW(*NUM,&vettore1[0],&vettore2[0]);
result10[VOLTE] = calcolaSZrc(*NUM,&vettore1[0],&vettore2[0]);
}
}
