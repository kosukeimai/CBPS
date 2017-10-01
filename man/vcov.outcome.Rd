\name{vcov_outcome}
\alias{vcov_outcome}
\alias{vcov_outcome.CBPSContinuous}
\title{Calculate Variance-Covariance Matrix for Outcome Model}
\description{
  \code{vcov_outcome} Returns the variance-covariance matrix of the main parameters of a fitted CBPS object.
}
\usage{\method{vcov_outcome}{CBPSContinuous}(object, Y, Z, delta, tol = 10^(-5), lambda = 0.01)}
\arguments{
  \item{object}{A fitted CBPS object.}
  \item{Y}{The outcome.}
  \item{Z}{The covariates (including the treatment and an intercept term) that predict the outcome.}
  \item{delta}{The coefficients from regressing Y on Z, weighting by the cbpsfit$weights.}
  \item{tol}{Tolerance for choosing whether to improve conditioning of the "M" matrix prior to conversion.  Equal to 1/(condition number), i.e. the smallest eigenvalue divided by the largest.}
  \item{lambda}{The amount to be added to the diagonal of M if the condition of the matrix is worse than tol.}
}
\details{This adjusts the standard errors of the weighted regression of Y on Z for uncertainty in the weights.}
\value{
	A matrix of the estimated covariances between the parameter estimates in the weighted outcome regression, adjusted for uncertainty in the weights.
}
\references{Lunceford and Davididian 2004.}

\author{Christian Fong, Chad Hazlett, and Kosuke Imai.}

\examples{
###
### Example: Variance-Covariance Matrix
###

##Load the LaLonde data
data(LaLonde)
## Estimate CBPS via logistic regression
fit <- CBPS(treat ~ age + educ + re75 + re74 + I(re75==0) + I(re74==0), 
		    data = LaLonde, ATT = TRUE)
## Get the variance-covariance matrix.
vcov(fit)
}