\name{plot.CBMSM}
\alias{plot.CBMSM}
\title{Plotting CBPS Estimation for Marginal Structural Models}
\description{Plots the absolute difference in standardized means before and after weighting.}
\usage{\method{plot}{CBMSM}(x, covars = NULL, silent = TRUE, boxplot = FALSE, ...)}
\arguments{
  \item{x}{an object of class \dQuote{CBMSM}.}
  \item{covars}{Indices of the covariates to be plotted (excluding the intercept).  For example, if only the first two covariates from \code{balance} are desired, set \code{covars} to 1:2.  The default is \code{NULL}, which plots all covariates.}
  \item{silent}{If set to \code{FALSE}, returns the absolute imbalance for each treatment history pair before and after weighting.  This helps the user to create his or her own customized plot. Default is \code{TRUE}, which returns nothing.}
  \item{boxplot}{If set to \code{TRUE}, returns a boxplot summarizing the imbalance on the covariates instead of a point for each covariate.  Useful if there are many covariates.}
  \item{...}{Additional arguments to be passed to plot.}
}
\details{Covariate balance is improved if the plot's points are below the plotted line of y=x.}
\value{
The x-axis gives the imbalance for each covariate-treatment history pair without any weighting, and the y-axis gives the imbalance for each covariate-treatment history pair after CBMSM weighting.  Imbalance is measured as the absolute difference in standardized means for the two treatment histories.  Means are standardized by the standard deviation of the covariate in the full sample.
}

\author{Marc Ratkovic and Christian Fong}

\seealso{\link{CBMSM}, \link{plot}}