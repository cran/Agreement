\name{agreement}
\alias{agreement}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Assessment of Agreement
}
\description{
Compute the concordance correlation coefficient (CCC), precision, accuracy, total deviation index (TDI), coverage probability (CP) and relative biased square (RBS) for the paired observations (test and target) based on the model proposed by Lin, Hedayat, Sinha and Yang (2002).
}
\usage{
agreement(x, y, error, target, CCC_a = 0.95, TDI_a, alpha = 0.05, CP_a = 0.9, H_label = "x", V_label = "y", min = NA, max = NA, by = NA, dec = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A vector of the readings for the target assay or rater, it will be shown in the horizontal axis of the agreement plot
}
  \item{y}{
A vector of the readings for the test assay or rater, it will be shown in the vertical axis of the agreement plot
}
  \item{error}{
Constant or proportional error structure. Use \code{error="const"} if the error structure is constant. Then, TDI is expressed as an absolute difference with the same measurement unit as the original data. Use \code{error="prop"} if the error structure is proportional, Then, TDI is expressed as a percent change. The natural log transformation to the data will be applied.
}
  \item{target}{
\code{target="random"} or \code{"fixed"}.
}
  \item{CCC_a}{
A CCC allowance, the default if \code{CCC_a=0.95}.
}
  \item{TDI_a}{
A TDI allowance which must be specified for computing CP, must be a percent value when \code{error="prop"} is specified or an absolute difference when \code{error="const"} is specified.
}
  \item{alpha}{
100(1-\code{alpha})% one-tailed confidence limit. The default is \code{alpha=0.05}.
}
  \item{CP_a}{
A CP allowance which must be specified for computing TDI. The default is \code{CP_a=0.9}.
}
  \item{H_label}{
Label for the horizontal axis (target) of the agreement plot. The default is \code{H_label="x"}.
}
  \item{V_label}{
Label for the vertical axis (test) of the agreement plot. The default is \code{V_label="y"}.
}
  \item{min}{
Minimum of the plotting range, if it is not specified, a range will be chosen automatically to show all the data point.
}
  \item{max}{
Maximum of the plotting range, if it is not specified, a range will be chosen automatically to show all the data point.
}
  \item{by}{
The increment of the plotting range. If \code{error="const"} is specified, it is the value of increment, for example, \code{by=5}. If \code{error="prop"} is specified, it is the log scale increments in between min and max. For example, if \code{min=1} and \code{max=64}, then \code{by=c(2,4,8,16,32)}. If \code{by} is not specified, it will be chosen automatically based on the error structure.
}
  \item{dec}{
Significant digits for the output. The default is \code{dec=2}.
}
}
\details{
The functions \code{summary} and \code{plot} can be used to print a summary of the results and plot the observations.
}
\value{
An object of class \code{report} is returned, which contains all the estimated statistics and their confidence limits.
}
\references{
Lawrence Lin, A. S Hedayat, Bikas Sinha, Min Yang. Journal of the American Statistical Association. March 1, 2002, 97(457): 257-270.
}
\author{
Yue Yu <yyu@imyy.net>
}

\seealso{
\code{\link{summary.agreement}}, \code{\link{summary.report}}, \code{\link{plot.report}}, \code{\link{html.report}}
}
\examples{
data(DCLHb);
HemocueAve <- apply(DCLHb[,c("HEMOCUE1","HEMOCUE2")],1,mean);
SigmaAve <- apply(DCLHb[,c("SIGMA1","SIGMA2")],1,mean);

agr=agreement(y=HemocueAve,x=SigmaAve,V_label="Hemacue",H_label="Sigma",min=0,max=2000,by=250,CCC_a=0.9775,CP_a=0.9,TDI_a=150,error="const",target="random",dec=3,alpha=0.05);
summary(agr);
plot(agr);
}
