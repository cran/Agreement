\name{summary.unified_agreement}
\alias{summary.unified_agreement}
\title{
Summarize a Unified Agreement Object
}
\description{
This function is a method for \code{summary} for a unified agreement object class \code{unified_agreement}, which can be get from function \code{unified.agreement}. 
}
\usage{
\method{summary}{unified_agreement}(object, \dots)
}
\arguments{
  \item{object}{
The name of an object of class \code{unified_agreement}.
}
  \item{\dots}{
Arguments passed to the \code{print} method.
}
}
\author{
Yue Yu <yyu@imyy.net>
}
\seealso{
\code{\link{unified.agreement}}
}
\examples{
data(DCLHb);
ua <- unified.agreement(dataset=DCLHb, var=c("HEMOCUE1","HEMOCUE2","SIGMA1","SIGMA2"), k=2, m=2, CCC_a_intra=0.9943, CCC_a_inter=0.9775, CCC_a_total=0.9775, CP_a=0.9, tran=1, TDI_a_intra=75, TDI_a_inter=150, TDI_a_total=150, error="const", dec=1, alpha=0.05);
summary(ua);
}
