\name{npCBPS}
\alias{npCBPS}
\alias{npCBPS.fit}
\title{Non-Parametric Covariate Balancing Propensity Score (npCBPS) Estimation}
\description{
  \code{npCBPS} is a method to estimate weights interpretable as (stabilized) inverse generlized propensity score weights, w_i = f(T_i)/f(T_i|X), without actually estimating a model for the treatment to arrive at f(T|X) estimates. In brief, this works by maximizing the empirical likelihood of observing the values of treatment and covariates that were observed, while constraining the weights to be those that (a) ensure balance on the covariates, and (b) maintain the original means of the treatment and covariates. In the continuous treatment context, this balance on covariates means zero correlation of each covariate with the treatment. Furthermore, we apply a Bayesian variant that allows the correlation of each covariate with the treatment to be slightly non-zero, as might be expected in a a given finite sample.}
\usage{
	  npCBPS(formula, data, na.action, corprior = 0.1, print.level = 0, ...)
	  npCBPS.fit(treat, X, corprior, print.level, ...)
}
\arguments{
  \item{formula}{An object of class \code{formula} (or one that can be coerced to that class): a symbolic description of the model to be fitted.}
  \item{data}{An optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from \code{environment(formula)}, typically the environment from which \code{CBPS} is called.}
  \item{na.action}{A function which indicates what should happen when the data contain NAs. The default is set by the na.action setting of options, and is na.fail if that is unset.}
  \item{corprior}{Prior hyperparameter controlling the expected amount of correlation between each covariate and the treatment. Specifically, the amount of correlation between the k-dimensional covariates, X, and the treatment T after weighting is assumed to have prior distribution MVN(0,sigma^2 I_k). We conceptualize sigma^2 as a tuning parameter to be used pragmatically. It's default of 0.1 ensures that the balance constraints are not too harsh, and that a solution is likely to exist. Once the algorithm works at such a high value of sigma^2, the user may wish to attempt values closer to 0 to get finer balance.}
  \item{print.level}{Controls verbosity of output to the screen while npCBPS runs. At the default of print.level=0, little output is produced. It print.level>0, it outputs diagnostics including  the log posterior (log_post), the log empirical likelihood associated with the weights (log_el), and the log prior probability of the (weighted) correlation of treatment with the covariates.}
  \item{treat}{A vector of treatment assignments.  Binary or multi-valued treatments should be factors.  Continuous treatments should be numeric.}
  \item{X}{A covariate matrix.}
  \item{...}{Other parameters to be passed.}
}
\details{Estimates non-parametric covariate balancing propensity score weights.}
\value{
  \item{weights}{The optimal weights}
  \item{y}{The treatment vector used}
  \item{x}{The covariate matrix}
  \item{model}{The model frame}
  \item{call}{The matched call}
  \item{formula}{The formula supplied}
  \item{data}{The data argument}
  \item{log.p.eta}{The log density for the (weighted) correlation of the covariates with the treatment, given the choice of prior (\code{corprior})}
  \item{log.el}{The log empirical likelihood of the observed data at the chosen set of IPW weights.}
  \item{eta}{A vector describing the correlation between the treatment and each covariate on the weighted data at the solution.}
  \item{sumw0}{The sum of weights, provided as a check on convergence. This is always 1 when convergence occurs unproblematically. If it differs from 1 substantially, no solution perfectly satisfying the conditions was found, and the user may consider a larger value of \code{corprior}.}
}
\references{Fong, Christian, Chad Hazlett, and Kosuke Imai.  ``Parametric and Nonparametric Covariate Balancing Propensity Score for General Treatment Regimes.''  Unpublished Manuscript. \url{http://imai.princeton.edu/research/files/CBGPS.pdf}}
\author{Christian Fong, Chad Hazlett, and Kosuke Imai}

\examples{
\dontrun{
##Generate data
data(LaLonde)

##Fitting npCBPS will take a couple of minutes
fit <- npCBPS(treat ~ age + educ + re75 + re74 + I(re75==0) + I(re74==0),
			  data = LaLonde, , corprior=.1/nrow(LaLonde))
summary(lm(re78 ~ treat + age + educ + re75 + re74 + I(re75==0) + I(re74==0),
		   data = LaLonde, weights = fit$weights))
}
}
