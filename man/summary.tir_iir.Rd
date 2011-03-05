\name{summary.tir_iir}
\alias{summary.tir_iir}
\title{
Summarize an TIR and IIR Object
}
\description{
This function is a method for \code{summary} for a unified agreement object class \code{tir_iir}, which can be get from function \code{TIR_IIR}. 
}
\usage{
\method{summary}{tir_iir}(object, \dots)
}
\arguments{
  \item{object}{
The name of an object of class \code{tir_iir}.
}
  \item{\dots}{
Arguments passed to the \code{print} method.
}
}
\author{
Yue Yu <yyu@imyy.net>
}
\seealso{
\code{\link{TIR_IIR}}
}
\examples{
data(DCLHb);

TI <- TIR_IIR(dataset=DCLHb, var=c("HEMOCUE1","HEMOCUE2","SIGMA1","SIGMA2"), k=2, m=2, TIR_test="1,2", TIR_ref="All", IIR_test="1", IIR_ref="2", error="const", alpha=0.05, TIR_a=2.25);
summary(TI);
}
