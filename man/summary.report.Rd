\name{summary.report}
\alias{summary.report}
\title{
Summarize an Agreement Object
}
\description{
This function is a method for \code{summary} for a unified agreement object class \code{report}, which can be get from function \code{agreement}. 
}
\usage{
\method{summary}{report}(object, \dots)
}
\arguments{
  \item{object}{
The name of an object of class \code{report}.
}
  \item{\dots}{
Arguments passed to the \code{print} method.
}
}
\author{
Yue Yu <yyu@imyy.net>
}
\seealso{
\code{\link{summary.agreement}}, \code{\link{agreement}}
}
\examples{
data(DCLHb);
HemocueAve <- apply(DCLHb[,c("HEMOCUE1","HEMOCUE2")],1,mean);
SigmaAve <- apply(DCLHb[,c("SIGMA1","SIGMA2")],1,mean);

agr=agreement(y=HemocueAve,x=SigmaAve,V_label="Hemacue",H_label="Sigma",min=0,max=2000,by=250,CCC_a=0.9775,CP_a=0.9,TDI_a=150,error="const",target="random",dec=3,alpha=0.05);
summary(agr);
}

