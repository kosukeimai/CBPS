\name{CBIV}
\alias{CBIV}
\title{Covariate Balancing Propensity Score for Instrumental Variable Estimates (CBIV)}
\description{
  \code{CBIV} estimates propensity scores for compliance status in an instrumental variables setup such that both covariate balance and prediction of treatment assignment are maximized.   The method, therefore, avoids an iterative process between model fitting and balance checking and implements both simultaneously.}
\usage{
	  CBIV(Tr, Z, X, iterations = 1000, 
	       method = "over", twostep = TRUE,  ...)
}
\arguments{
  \item{Tr}{A binary treatment variable.}
  \item{Z}{A binary encouragement variable.}
  \item{X}{A pre-treatment covariate matrix.}
  \item{iterations}{An optional parameter for the maximum number of iterations for the optimization.  Default is 1000.}
  \item{method}{Choose "over" to fit an over-identified model that combines the propensity score and covariate balancing conditions; choose "exact" to fit a model that only contains the covariate balancing conditions.  Our simulations suggest that "over" dramatically outperforms "exact."}
  \item{twostep}{Default is \code{TRUE} for a two-step GMM estimator, which will run substantially faster than continuous-updating.  Set to \code{FALSE} to use the continuous-updating GMM estimator.}
  \item{...}{Other parameters to be passed through to \code{optim()}.}
}
\details{Fits covariate balancing propensity scores for generalizing local average treatment effect estimates obtained from instrumental variables analysis.}
\value{
  \item{coefficients}{A named matrix of coefficients, where the first column gives the complier coefficients and the second column gives the always-taker coefficients.}
  \item{fitted.values}{The fitted N x 3 compliance score matrix.  The first column gives the estimated probability of being a complier, the second column gives the estimated probability of being an always-taker, and the third column gives the estimated probability of being a never-taker.}
  \item{weights}{The optimal weights: the reciprocal of the probability of being a complier.}
  \item{deviance}{Minus twice the log-likelihood of the CBIV fit.}
  \item{converged}{Convergence value.  Returned from the call to \code{optim()}.}
  \item{J}{The J-statistic at convergence}
  \item{df}{The number of linearly independent covariates.}
  \item{bal}{The covariate balance associated with the optimal weights, calculated as the GMM loss of the covariate balance conditions.}
 }
\references{Imai, Kosuke and Marc Ratkovic.  2014. ``Covariate Balancing Propensity Score.'' Journal of the Royal Statistical Society, Series B (Statistical Methodology). \url{http://imai.princeton.edu/research/CBPS.html}}
\author{Christian Fong}


\examples{
###
### Example: propensity score matching
### (Need to fix when we have an actual example).

##Load the LaLonde data
data(LaLonde)
## Estimate CBPS
fit <- CBPS(treat ~ age + educ + re75 + re74 + 
			I(re75==0) + I(re74==0), 
			data = LaLonde, ATT = TRUE)
summary(fit)
}