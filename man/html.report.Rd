\name{html.report}
\alias{html.report}
\title{
Generate a HTML Report for Agreement Statistics
}
\description{
Generate a HTML report for an agreement object class \code{report}. Contain the plot of data and all agreement statistics and confidence limits. Require package \code{R2HTML}.
}
\usage{
\method{html}{report}(result, dir = getwd(), file = "report", CSS = "R2HTML", \dots)
}
\arguments{
  \item{result}{
The name of an object of class \code{report}, which can be get from function \code{agreement}.
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
\code{\link{agreement}}
}

