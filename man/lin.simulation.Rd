\name{lin.simulation}
\alias{lin.simulation}
\title{ Performs a summary of the simulation }
\description{
  In order to analyze the problem we take a bivariate gaussian distribution on \emph{(X,Y)}.
  It is true under the null hypothesis (no agreement) and the alternative (yes agreement).
}
\usage{
lin.simulation(NUM_CAMP = 5000, NUM = 30, matH0, matH1, underH0 = TRUE, ALPHA_CI = 0.05, la_CP1 = 0.9)
}
\arguments{
  \item{NUM_CAMP}{ number of samples to simulate. Its default value is \emph{5000} }
  \item{NUM}{ sample size. Its default value is \emph{30} }
  \item{matH0}{ matrix of parameters under null hypothesis. It has \emph{2} rows and \emph{3} columns. See \bold{Details}. }
  \item{matH1}{ matrix of parameters under alternative hypothesis . It has exactly the same structure of \code{matH0} }
  \item{underH0}{ logical parameter to determine what condition to simulate. Its default value is \code{TRUE} (simulation under null hypothesis) }
  \item{ALPHA_CI}{ level of significance. Its default value is \emph{0.05} }
  \item{la_CP1}{ the threshold used for TDI. Its default value is \emph{0.9} }
}
\details{
  The first column of \code{table} object is the theoretical value of the indices (\bold{Th val}) while
  the second column (\bold{Thr}) represent the threshold used to determine the rejection region.
  Theoretical values of \eqn{\alpha}{alpha} and \eqn{1 - \beta}{1 - beta} are reported in the third
  column (\bold{Th prob}). The fourth column (\bold{Mean of est}) represents the inverse transformation
  of the mean estimate of the agreement measure. We expect the first and the fourth columns to be similar
  in order to consider the estimate robust. The same conclusion is made between the fifth and the sixth
  columns which represent the standard deviation of the transformation (\bold{Std of est}) and the mean
  of the standard deviation (\bold{Mean of std}) respectively. In the seventh column (\bold{Prop rej}) it
  is calculated the proportion between \code{NUM_CAMP} runs fall in the rejection region. If we simulate
  under \eqn{H_0}{H0} then we expect that this value is about \eqn{\alpha = 0.05}{alpha = 0.05} (type one
  error probability) while we expect it is about the true value \eqn{1 - \beta}{1 - beta} (power) if we
  simulate under \eqn{H_1}{H1}.\cr \cr \code{matH0} is a matrix of six elements. In the first
  row we have the variance of \emph{X}, the covariance between \emph{X} and \emph{Y} and the expected value of \emph{X} respectively.
  In the second row we have the covariance between \emph{X} and \emph{Y}, the variance of \emph{Y} and the expected value of \emph{Y} respectively.
  \code{matH1} has exactly the same structure of \code{matH0}. 
}
\value{
  \item{ table }{ it is a matrix with \emph{6} rows and \emph{7} columns. Each row represents a measure of agreement and each column a summary for the simulation. See \bold{Details}. }
  \item{ underH0 }{ see above }
  \item{ matH0 }{ see above }
  \item{ matH1 }{ see above }
  \item{ NUM\_CAMP}{ see above }
  \item{ NUM }{ see above }
  \item{ alpha }{ see above }
  \item{ rho }{ is the value for the correlation coefficient under null hypothesis (if \code{underH0 = TRUE}) or alternative hypothesis (if \code{underH0 = FALSE})}
}
\references{
D. G. Altman, J. M. Bland (1983):
Measurement in Medicine: The Analysis of Method Comparison Studies,
The Statistician, \bold{32}, 302--317

L. Lin (1989):
A Concordance Correlation Coefficient to Evaluate Reproducibility,
Biometrics, \bold{45}, 255--258

L. Lin, A. Heyadat, B. Sinha, M. Yang (2002):
Statistical methods in assessing agreement: Models, issues, and tools,
JASA, \bold{97}, 257--270
}
\author{ Fabio Frascati \email{fabio.frascati@studenti.unimi.it}}
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
%\seealso{ \code{\link{help}}, ~~~ }
\examples{
## we define the matrix of parameters under H0
sigma2x0   <- 1 / 1.15
sigma2y0   <- 1.15
covxy0     <- 0.95 * sqrt(1 / 1.15 * 1.15)
mux0       <- 0
muy0       <- 0.15
matH0      <- matrix(0,nrow = 2,ncol = 3)
matH0[1,1] <- sigma2x0
matH0[1,2] <- covxy0
matH0[1,3] <- mux0
matH0[2,1] <- covxy0
matH0[2,2] <- sigma2y0
matH0[2,3] <- muy0
## we define the matrix of parameters under H1
sigma2x1   <- 1 / 1.1
sigma2y1   <- 1.1
covxy1     <- 0.9662055 * sqrt(1 / 1.1 * 1.1)
mux1       <- 0
muy1       <- 0.1
matH1      <- matrix(0,nrow = 2,ncol = 3)
matH1[1,1] <- sigma2x1
matH1[1,2] <- covxy1
matH1[1,3] <- mux1
matH1[2,1] <- covxy1
matH1[2,2] <- sigma2y1
matH1[2,3] <- muy1
## we run lin.simulation()
lin.simulation(matH0 = matH0,matH1 = matH1,NUM = 30,NUM_CAMP = 5000,underH0 = TRUE)
}
\keyword{ models }
