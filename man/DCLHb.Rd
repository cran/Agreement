\name{DCLHb}
\alias{DCLHb}
\docType{data}
\title{
DCLHb dataset
}
\description{
DCLHb is a treatment solution contain oxygen carrying hemoglobin. The DCLHb level is routinely measured by the Sigma method, and the simpler HemoCue method was modified to reproduce the DCLHb values of the Sigma method. The DCLHb values of each sample were measured by both methods twice.
}
\usage{data(DCLHb)}
\format{
  A data frame with 299 observations on 4 variables, 2 variables for Sigma method and 2 variables for HemoCue method.
}
\references{
Lawrence Lin, A. S Hedayat, Bikas Sinha, Min Yang. Journal of the American Statistical Association. March 1, 2002, 97(457): 257-270.
}
\examples{
data(DCLHb);
HemocueAve <- apply(DCLHb[,c("HEMOCUE1","HEMOCUE2")],1,mean);
SigmaAve <- apply(DCLHb[,c("SIGMA1","SIGMA2")],1,mean);

agr=agreement(y=HemocueAve,x=SigmaAve,V_label="Hemacue",H_label="Sigma",min=0,max=2000,by=250,CCC_a=0.9775,CP_a=0.9,TDI_a=150,error="const",target="random",dec=3,alpha=0.05);
summary(agr);
}
