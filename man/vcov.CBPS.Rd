\name{vcov.CBPS}
\alias{vcov.CBPS}
\title{Calculate Variance-Covariance Matrix for a Fitted CBPS Object }
\description{
  \code{vcov.CBPS} Returns the variance-covariance matrix of the main parameters of a fitted CBPS object.
}
\usage{\method{vcov}{CBPS}(object,...)}
\arguments{
  \item{object}{An object of class \code{formula} (or one that can be coerced to that class): a symbolic description of the model to be fitted.}
  \item{...}{Additional arguments to be passed to vcov.CBPS}
}
\details{This is the CBPS implementation of the generic function vcov().}
\value{
	A matrix of the estimated covariances between the parameter estimates in the linear or non-linear predictor of the model. 
}
\references{This documentation is modeled on the documentation of the generic \link{vcov}.}

\author{Christian Fong, Marc Ratkovic, and Kosuke Imai.}

\seealso{\link{vcov}}

\examples{
\dontrun{
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
}