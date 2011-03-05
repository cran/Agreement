\name{TIR_IIR}
\alias{TIR_IIR}
\title{
A Comparative Model
}
\description{
This function is used for comparing the total-rater agreement to intra-rater precision, and for comparing intra-rater precision among selected raters.
}
\usage{
TIR_IIR(dataset, var, k, m, TIR_test, TIR_ref, IIR_test, IIR_ref, error = "const", alpha = 0.05, TIR_a)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
Name of your dataset.
}
  \item{var}{
Name of the dependent variables used in the dataset. The variables should be sorted by methods. For example, \code{var=c("Method1_Rep1","Method1_Rep2","Method2_Rep1","Method2_Rep2")}.
}
  \item{k}{
Number of methods/raters/instruments/assay, etc.
}
  \item{m}{
Number of replications for each methods/raters/instruments/assay, etc.
}
  \item{TIR_test}{
The selected test raters for calculating TIR, and must input in the format of a string \code{"1,2,...,k"}, where 1 represents rater 1, 2 represents rater 2, and k represents for rater k, etc. For example, \code{TIR_test="1,2"} means the selected test raters for calculating TIR are the first and second rater. When calculating multiple TIRs, the test raters for calculating each TIR must be written in a vector. For example, when k=3, we could specify \code{TIR_test=c("1,3","1,2,3","3","2","1,2")}. 
}
  \item{TIR_ref}{
The selected reference raters for calculating TIR that correspond to the \code{TIR_test}. It has the same format as \code{TIR_test}. If \code{TIR_ref="all"} is specified, then the intra-raters of all raters will be used as the denominator.
}
  \item{IIR_test}{
The selected test raters for calculating IIR. It has the same format as \code{TIR_test}, but can not be \code{"all"}.
}
  \item{IIR_ref}{
The selected reference raters for calculating IIR. It has the same format as \code{TIR_test}, but can not be \code{"all"}.
}
  \item{error}{
\code{error="const"} for the constant error structure for continuous data. \code{error="prop"} for the proportional error structure for continuous data. Here, log transformation to data will be applied for the continuous data. For categorical data, use \code{error="const"}.
}
  \item{alpha}{
1-\code{alpha} one-tailed upper confidence limit for TIR or two-tailed confidence interval for IIR. The default is \code{alpha=0.05}.
}
  \item{TIR_a}{
Allowance for TIR.
}
}
\details{
When \code{TIR_ref} is not specified as \code{"all"}, each TIR is computed as the total MSD of test vs selected reference raters relative to the intra MSD of the selected reference raters. When \code{TIR_ref} is specified as \code{"all"}, the function would assess the average of the total MSD of all raters relative to the average of intra MSD of all raters. For example, if \code{TIR_test=c("1,3","1,2,3","3","2","1,2")} and \code{TIR_ref=c("2","all","1,2","1","1,2")}. For the first pair in \code{TIR_test} and \code{TIR_ref}, the function would assess the average of the total MSD of "raters 1 vs 2 and raters 3 vs 2" relative to the intra MSD of "rater 2". For the second pair, the function would assess the average of the total MSD of all raters relative to the average of intra MSD of all raters. For the third pair, the function would assess the average of the total MSD of "raters 3 vs 1 and raters 3 vs 2" relative to the average of intra MSD of "raters 1 and 2". For the fourth pair, the function would assess the total MSD of "raters 2 and 1" relative to the intra MSD of "rater 1". For the fifth pair, the function would assess the total MSD of "raters 1 and 2" relative to the average of intra MSD of "raters 1 and 2".
}
\value{
An object of class \code{tir_iir} is returned, which contains all the estimated statistics and their confidence limits.
}
\author{
Yue Yu <yyu@imyy.net>
}
\seealso{
\code{\link{summary.tir_iir}}, \code{\link{TIR}}, \code{\link{IIR}}
}
\examples{
data(DCLHb);

TI <- TIR_IIR(dataset=DCLHb, var=c("HEMOCUE1","HEMOCUE2","SIGMA1","SIGMA2"), k=2, m=2, TIR_test="1,2", TIR_ref="All", IIR_test="1", IIR_ref="2", error="const", alpha=0.05, TIR_a=2.25);
summary(TI);
}
