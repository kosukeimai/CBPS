\name{CBPS}
\alias{CBPS}
\alias{CBPS.fit}
\alias{print.CBPS}
\title{Covariate Balancing Propensity Score (CBPS) Estimation}
\description{
  \code{CBPS} estimates propensity scores such that both covariate balance and prediction of treatment assignment are maximized.   The method, therefore, avoids an iterative process between model fitting and balance checking and implements both simultaneously.
For cross-sectional data, the method can take continuous treatments and treatments with a control (baseline) condition and either 1, 2, or 3 distinct treatment conditions.}
\usage{
	  CBPS(formula, data, na.action, ATT = 1, iterations = 1000, 
	       standardize = TRUE, method = "over", twostep = TRUE, 
		   baseline.formula = NULL, diff.formula = NULL, ...)
	  CBPS.fit(treat, X, baselineX, diffX, ATT, method, 
			   iterations, standardize, twostep, ...)
}
\arguments{
  \item{formula}{An object of class \code{formula} (or one that can be coerced to that class): a symbolic description of the model to be fitted.}
  \item{data}{An optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from \code{environment(formula)}, typically the environment from which \code{CBPS} is called.}
  \item{na.action}{A function which indicates what should happen when the data contain NAs. The default is set by the na.action setting of options, and is na.fail if that is unset.}
  \item{ATT}{Default is 1, which finds the average treatment effect on the treated interpreting the second level of the treatment factor as the treatment.  Set to 2 to find the ATT interpreting the first level of the treatment factor as the treatment.  Set to 0 to find the average treatment effect. For non-binary treatments, only the ATE is available.}
  \item{iterations}{An optional parameter for the maximum number of iterations for the optimization.  Default is 1000.}
  \item{standardize}{Default is \code{TRUE}, which normalizes weights to sum to 1 within each treatment group.  For continuous treatments, normalizes weights to sum up to 1 for the entire sample.  Set to \code{FALSE} to return Horvitz-Thompson weights.}
  \item{method}{Choose "over" to fit an over-identified model that combines the propensity score and covariate balancing conditions; choose "exact" to fit a model that only contains the covariate balancing conditions.}
  \item{twostep}{Default is \code{TRUE} for a two-step estimator, which will run substantially faster than continuous-updating.  Set to \code{FALSE} to use the continuous-updating estimator described by Imai and Ratkovic (2014).}
  \item{treat}{A vector of treatment assignments.  Binary or multi-valued treatments should be factors.  Continuous treatments should be numeric.}
  \item{X}{A covariate matrix.} 
  \item{baseline.formula}{Used only to fit iCBPS (see Fan et al). Currently only works with binary treatments.  A formula specifying the balancing covariates in the baseline outcome model, i.e., E(Y(0)|X).}
  \item{diff.formula}{Used only to fit iCBPS (see Fan et al).  Currently only works with binary treatments.  A formula specifying the balancing covariates in the difference between the treatment and baseline outcome model, i.e., E(Y(1)-Y(0)|X).}
  \item{baselineX}{Similar to \code{baseline.formula}, but in matrix form.}
  \item{diffX}{Similar to \code{diff.formula}, but in matrix form.}
  \item{...}{Other parameters to be passed through to \code{optim()}.}
}
\details{Fits covariate balancing propensity scores.}
\value{
  \item{fitted.values}{The fitted propensity score}
  \item{deviance}{Minus twice the log-likelihood of the CBPS fit}
  \item{weights}{The optimal weights.  Let \eqn{\pi_i = f(T_i | X_i)}{\pi_i = f(T_i | X_i)}.  For binary ATE, these are given by \eqn{\frac{T_i}{\pi_i} + \frac{(1 - T_i)}{(1 - \pi_i)}}{T_i/\pi_i + (1 - T_i)/(1 - \pi_i)}.  For binary ATT, these are given by \eqn{\frac{n}{n_t} * \frac{T_i - \pi_i}{1 - \pi_i}}{n/n_t * (T_i - \pi_i)/(1 - \pi_i)}.  For multi_valued treatments, these are given by \eqn{\sum_{j=0}^{J-1} T_{i,j} / \pi_{i,j}}{\sum_{j=0}^{J-1} T_i,j / \pi_i,j}.  For continuous treatments, these are given by \eqn{\frac{f(T_i)}{f(T_i | X_i)}}{f(T_i) / f(T_i | X_i) }.  These expressions for weights are all before standardization (i.e. with standardize=\code{FALSE}).  Standardization will make weights sum to 1 within each treatment group.  For continuous treatment, standardization will make all weights sum to 1.}
  \item{y}{The treatment vector used}
  \item{x}{The covariate matrix}
  \item{model}{The model frame}
  \item{converged}{Convergence value.  Returned from the call to \code{optim()}.}
  \item{call}{The matched call}
  \item{formula}{The formula supplied}
  \item{data}{The data argument}
  \item{coefficients}{A named vector of coefficients}
  \item{sigmasq}{The sigma-squared value, for continuous treatments only}
  \item{J}{The J-statistic at convergence}
  \item{mle.J}{The J-statistic for the parameters from maximum likelihood estimation}
  \item{var}{The covariance matrix, evaluated numerically from \code{optim()}.}
 }
\references{Imai, Kosuke and Marc Ratkovic.  2014. ``Covariate Balancing Propensity Score.'' Journal of the Royal Statistical Society, Series B (Statistical Methodology). \url{http://imai.princeton.edu/research/CBPS.html}			\cr
			Fong, Christian, Chad Hazlett, and Kosuke Imai.  ``Parametric and Nonparametric Covariate Balancing Propensity Score for General Treatment Regimes.''  Unpublished Manuscript. \url{http://imai.princeton.edu/research/files/CBGPS.pdf}	\cr
			Fan, Jianqing and Imai, Kosuke and Liu, Han and Ning, Yang and Yang, Xiaolin. ``Improving Covariate Balancing Propensity Score: A Doubly Robust and Efficient Approach.'' Unpublished Manuscript. \url{http://imai.princeton.edu/research/CBPStheory.html}}
\author{Christian Fong, Marc Ratkovic, Kosuke Imai, and Xiaolin Yang; The CBPS function is based on the code for version 2.15.0 of the glm function implemented in the stats package, originally written by Simon Davies.  This documentation is likewise modeled on the documentation for glm and borrows its language where the arguments and values are the same.}


\seealso{\link{summary.CBPS}}
\examples{
###
### Example: propensity score matching
###

##Load the LaLonde data
data(LaLonde)
## Estimate CBPS
fit <- CBPS(treat ~ age + educ + re75 + re74 + 
			I(re75==0) + I(re74==0), 
			data = LaLonde, ATT = TRUE)
summary(fit)
\dontrun{
## matching via MatchIt: one to one nearest neighbor with replacement
library(MatchIt)
m.out <- matchit(treat ~ fitted(fit), method = "nearest", 
				 data = LaLonde, replace = TRUE)

### Example: propensity score weighting 
###
## Simulation from Kang and Shafer (2007).
set.seed(123456)
n <- 500
X <- mvrnorm(n, mu = rep(0, 4), Sigma = diag(4))
prop <- 1 / (1 + exp(X[,1] - 0.5 * X[,2] + 
			 0.25*X[,3] + 0.1 * X[,4]))
treat <- rbinom(n, 1, prop)
y <- 210 + 27.4*X[,1] + 13.7*X[,2] + 13.7*X[,3] + 13.7*X[,4] + rnorm(n)

##Estimate CBPS with a misspecified model
X.mis <- cbind(exp(X[,1]/2), X[,2]*(1+exp(X[,1]))^(-1)+10, 
			  (X[,1]*X[,3]/25+.6)^3, (X[,2]+X[,4]+20)^2)
fit1 <- CBPS(treat ~ X.mis, ATT = 0)
summary(fit1)
	
## Horwitz-Thompson estimate
mean(treat*y/fit1$fitted.values)
## Inverse propensity score weighting
sum(treat*y/fit1$fitted.values)/sum(treat/fit1$fitted.values)

rm(list=c("y","X","prop","treat","n","X.mis","fit1"))

### Example: Continuous Treatment
set.seed(123456)
n <- 1000
X <- mvrnorm(n, mu = rep(0,2), Sigma = diag(2))
beta <- rnorm(ncol(X)+1, sd = 1)
treat <- cbind(1,X)\%*\%beta + rnorm(n, sd = 5)

treat.effect <- 1
effect.beta <- rnorm(ncol(X))
y <- rbinom(n, 1, (1 + exp(-treat.effect*treat - 
				   X\%*\%effect.beta))^-1)

fit2 <- CBPS(treat ~ X)
summary(fit2)
summary(glm(y ~ treat + X, weights = fit2$weights, 
			family = "quasibinomial"))

rm(list=c("n", "X", "beta", "treat", "treat.effect",
		  "effect.beta", "y", "fit2"))

### Simulation example: Improved CBPS (or iCBPS) from Fan et al
set.seed(123456)
n <- 500
X <- mvrnorm(n, mu = rep(0, 4), Sigma = diag(4))
prop <- 1 / (1 + exp(X[,1] - 0.5 * X[,2] + 0.25*X[,3] + 0.1 * X[,4]))
treat <- rbinom(n, 1, prop)
y1 <- 210 + 27.4*X[,1] + 13.7*X[,2] + 13.7*X[,3] + 13.7*X[,4] + rnorm(n)
y0 <- 210 + 13.7*X[,2] + 13.7*X[,3] + 13.7*X[,4] + rnorm(n)
##Estimate iCBPS with a misspecificied model
X.mis <- cbind(exp(X[,1]/2), X[,2]*(1+exp(X[,1]))^(-1)+10, 
			   (X[,1]*X[,3]/25+.6)^3, (X[,2]+X[,4]+20)^2)
fit1 <- CBPS(treat ~ X.mis, baseline.formula=~X.mis[,2:4], 
			 diff.formula=~X.mis[,1], ATT = FALSE)
summary(fit1)
}
}
