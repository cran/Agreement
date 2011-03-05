\name{IIR}
\alias{IIR}
\title{
Intra-Intra Ratio (IIR)
}
\description{
Intra-Intra Ratio (IIR) compares the intra-rater precision among selected raters.
}
\usage{
IIR(data, k, m, testindex, refindex, error, alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Name of your dataset.
}
  \item{k}{
Number of methods/raters/instruments/assay, etc.
}
  \item{m}{
Number of replications for each methods/raters/instruments/assay, etc.
}
  \item{testindex}{
The index for the test raters, it must be a numerical vector. For example, if \code{k=3} and \code{testindex=c(1,2)}, the first and second rater are selected as the test raters for calculating IIR.
}
  \item{refindex}{
The index for the reference raters, it must be a numerical vector. For example, if \code{k=3} and \code{refindex=3}, the third rater is selected as the reference raters for calculating TIR.
}
  \item{error}{
\code{error="const"} for the constant error structure for continuous data. \code{error="prop"} for the proportional error structure for continuous data. Here, log transformation to data will be applied for the continuous data. For categorical data, use \code{error="const"}.
}
  \item{alpha}{
100(1-\code{alpha})% two-tailed confidence interval for IIR. The default is \code{alpha=0.05}.
}
}
\note{
The selected two sets of test and reference raters should be mutually exclusive.
}
\value{
Return a list of values
\item{IIR}{The estimator of IIR.}
\item{IIR_upper}{The upper 1-\code{alpha}/2 confidence limit for IIR.}
\item{IIR_lower}{The lower 1-\code{alpha}/2 confidence limit for IIR.}
}
\author{
Yue Yu <yyu@imyy.net>
}
\seealso{
\code{\link{TIR_IIR}}, \code{\link{TIR}}
}
\examples{
data(DCLHb);
IIR(DCLHb, k=2,m=2,testindex=1,refindex=2,error="const",alpha=0.05);
}

