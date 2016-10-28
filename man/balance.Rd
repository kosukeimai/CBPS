\name{balance}
\alias{balance}
\alias{balance.npCBPS}
\alias{balance.CBPS}
\alias{balance.CBMSM}
\title{Optimal Covariate Balance}
\description{Returns the mean and standardized mean associated with each treatment group, before and after weighting.}
\usage{\method{balance}{CBPS}(object, ...)}
\arguments{
  \item{object}{A CBPS, npCBPS, or CBMSM object.}
  \item{...}{Additional arguments to be passed to balance.}
}
\details{For binary and multi-valued treatments as well as marginal structural models, each of the matrices' rows are the covariates and whose columns are the weighted mean, and standardized mean associated with each treatment group.  The standardized mean is the weighted mean divided by the standard deviation of the covariate for the whole population.  For continuous treatments, returns the absolute Pearson correlation between the treatment and each covariate.}
\value{
  Returns a list of two matrices, "original" (before weighting) and "balanced" (after weighting).  
}
\author{Christian Fong, Marc Ratkovic, and Kosuke Imai.}

\examples{

\dontrun{
###
### Example: Optimal Covariate Balance
###
# Run for 4-treatment case
set.seed(1)
# Generate random X and underlying coefficients for probability.  
# Determine probs and treatments.
X<-cbind(rep(1,1000), mvrnorm(1000,c(0,0,0), 
		 Sigma=matrix(c(5,.5,-.03,.5,1,-.27,-.03,-.27,1),3,3)))
beta<-matrix(rnorm(12),4,3)
baseline.prob<-apply(X\%*\%beta,1,function(x) (1+sum(exp(x)))^-1)
prob<-cbind(baseline.prob, exp(X\%*\%beta[,1])*baseline.prob, 
			exp(X\%*\%beta[,2])*baseline.prob, 
			exp(X\%*\%beta[,3])*baseline.prob)
treat.latent<-runif(1000)
treat<-factor(ifelse(treat.latent < prob[,1], 1, 
			  ifelse(treat.latent < (prob[,1] + prob[,2]), 2, 
			  ifelse(treat.latent < (prob[,1] + prob[,2] + prob[,3]), 
			  3, 4))))
X<-X[,-1]
fit4<-CBPS(treat ~ X, ATT = FALSE)
balance(fit4)
}
}