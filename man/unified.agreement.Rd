\name{unified.agreement}
\alias{unified.agreement}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A Unified Approach for Assessment of Agreement
}
\description{
The unified approach calculates the agreement statistics for both continuous and categorical data to cover multiple readings from each of the n subjects.
}
\usage{
unified.agreement(dataset, var = NA, k, m, CCC_a = 0.9775, CCC_a_intra = 0.995, CCC_a_inter = 0.9775, CCC_a_total = 0.9775, CP_a = 0.9, TDI_a = 150, TDI_a_intra = 75, TDI_a_inter = 150, TDI_a_total = 150, tran = 0, error = "prop", dec = 4, alpha = 0.025, iter = 35, toe = 1e-10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
Name of your dataset.
}
  \item{var}{
Name of the dependent variables used in the dataset, if it is not specified, all the variables in the dataset will be used.
}
  \item{k}{
Number of methods/raters/instruments/assay, etc.
}
  \item{m}{
Number of replications for each methods/raters/instruments/assay, etc.
}
  \item{CCC_a}{
CCC allowance when there is no replications, the default is \code{CCC_a=0.9775}.
}
  \item{CCC_a_intra}{
Intra CCC allowance.
}
  \item{CCC_a_inter}{
Inter CCC allowance.
}
  \item{CCC_a_total}{
Total CCC allowance.
}
  \item{CP_a}{
CP allowance for continuous data.
}
  \item{TDI_a}{
TDI allowance when there is no replications for continuous data. The default is \code{TDI_a=150}.
}
  \item{TDI_a_intra}{
Intra TDI allowance for continuous data.
}
  \item{TDI_a_inter}{
Inter TDI allowance for continuous data.
}
  \item{TDI_a_total}{
Total TDI allowance for continuous data.
}
  \item{tran}{
If transformation is used for the dataset. \code{tran=0}: no transformation will be used for statistical inference. \code{tran=1} transformation such as Z, logit and log will be used for statistical inference. \code{tran=1} can be used for categorical data but the TDI and CP outputs would become irrelevant. Therefore, \code{tran=0} is recommended for all categorical data and \code{tran=1} is recommended for all continuous data.
}
  \item{error}{
\code{error="const"}: constant error structure for continuous data. When \code{error="const"}, TDI is an absolute difference with the same measurement unit as for the original data. \code{error="prop"}: proportional error structure for continuous data. When \code{error="prop"}, TDI is the percent change. Log transformation will be applied to original data. For categorical data, use \code{error="const"}.
}
  \item{dec}{
Significant digits after the decimal point printed for TDI when the \code{error="const"} is specified. The default is \code{dec=2}.
}
  \item{alpha}{
100(1-\code{alpha})% one-tailed confidence limit. The default is \code{alpha=0.05}.
}
  \item{iter}{
Number of iterations in the GEE step. The default is \code{iter=35}.
}
  \item{toe}{
The error tolerance of the GEE iteration. The default is \code{toe = 1e-10}.
}
}
\value{
An object of class \code{unified_agreement} is returned, which contains all the estimated statistics and their confidence limits.
}
\author{
Yue Yu <yyu@imyy.net>
}

\seealso{
\code{\link{agreement}}, \code{\link{summary.unified_agreement}}, \code{\link{html.unified_agreement}}
}
\examples{
data(DCLHb);
ua <- unified.agreement(dataset=DCLHb, var=c("HEMOCUE1","HEMOCUE2","SIGMA1","SIGMA2"), k=2, m=2, CCC_a_intra=0.9943, CCC_a_inter=0.9775, CCC_a_total=0.9775, CP_a=0.9, tran=1, TDI_a_intra=75, TDI_a_inter=150, TDI_a_total=150, error="const", dec=1, alpha=0.05);
summary(ua);
}
