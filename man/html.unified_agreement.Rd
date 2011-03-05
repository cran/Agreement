\name{html.unified_agreement}
\alias{html.unified_agreement}
\title{
Generate a HTML Report for Unified Agreement Statistics
}
\description{
Generate a HTML report for an agreement object class \code{unified_agreement}. Contain the plot of data and all agreement statistics. Require package \code{R2HTML}.
}
\usage{
\method{html}{unified_agreement}(result, dir = getwd(), file = "report", CSS = "R2HTML", \dots)
}
\arguments{
  \item{result}{
The name of an object of class \code{unified_agreement}, which can be get from function \code{unified.agreement}.
}
  \item{dir}{
Directory to store the output file. The default is the current working directory.
}
  \item{file}{
File name.
}
  \item{CSS}{
Name of CSS file used. The default is \code{R2HTML}.
}
  \item{\dots}{
Arguments passed to the \code{HTMLInitFile} method.
}
}
\details{
In order to use the CSS style, the CSS file has to be placed in the same directory as the report file. The default is to use the simple CSS file \code{R2HTML.CSS} included in the package \code{R2HTML}.
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
html.unified_agreement(ua);
}
