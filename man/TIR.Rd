\name{TIR}
\alias{TIR}
\title{
Total-Intra Ratio (TIR)
}
\description{
Total-Intra Ratio (TIR) is the ratio of MSD(total) and MSD(intra), it compares the total-rater agreement to intra-rater precision.
}
\usage{
TIR(data, k, m, testindex, refindex, error, alpha = 0.05)
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
The index for the test raters, it can be a numerical vector or a character value \code{"all"}. For example, if \code{testindex=c(1,2)}, the first and second rater are selected as the test raters for calculating TIR.
}
  \item{refindex}{
The index for the reference raters, it can be a numerical vector or a character value \code{"all"}. For example, if \code{refindex="all"}, all raters are selected as the reference raters for calculating TIR.
}
  \item{error}{
\code{error="const"} for the constant error structure for continuous data. \code{error="prop"} for the proportional error structure for continuous data. Here, log transformation to data will be applied for the continuous data. For categorical data, use \code{error="const"}.
}
  \item{alpha}{
100(1-\code{alpha})% one-tailed upper confidence limit for TIR. The default is \code{alpha=0.05}.
}
}
\value{
Return a list of values
\item{TIR}{The estimator of TIR.}
\item{TIR_upper}{The upper 1-\code{alpha} confidence limit for TIR.}
\item{TIR_lower}{The lower 1-\code{alpha} confidence limit for TIR.}
}
\author{
Yue Yu <yyu@imyy.net>
}
\seealso{
\code{\link{TIR_IIR}}, \code{\link{IIR}}
}
\examples{
data(DCLHb);
TIR(DCLHb, k=2,m=2,testindex=c(1,2),refindex="all",error="const",alpha=0.05);
}

