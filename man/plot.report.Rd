\name{plot.report}
\alias{plot.report}
\title{
Plotting an Agreement Object
}
\description{
Perform the x-y plot for the agreement data. The x, y label, plotting range and axis increments are in the \code{report} object.
}
\usage{
\method{plot}{report}(x, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
The name of an object of class \code{report}. It can be get from function \code{agreement}.
}
  \item{\dots}{
Arguments passed to the \code{plot.method}.
}
}
\author{
Yue Yu <yyu@imyy.net>
}
\seealso{
\code{\link{agreement}}
}
\examples{
data(DCLHb);
HemocueAve=apply(DCLHb[,c("HEMOCUE1","HEMOCUE2")],1,mean);
SigmaAve=apply(DCLHb[,c("SIGMA1","SIGMA2")],1,mean);

agr=agreement(y=HemocueAve,x=SigmaAve,V_label="Hemacue",H_label="Sigma",min=0,max=2000,by=250,CCC_a=0.9775,CP_a=0.9,TDI_a=150,error="const",target="random",dec=3,alpha=0.05);
plot(agr);
}
