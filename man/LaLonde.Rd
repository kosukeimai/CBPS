\name{LaLonde}
\docType{data}
\alias{LaLonde}
\title{LaLonde Data for Covariate Balancing Propensity Score}
\description{This data set gives the outcomes a well as treatment assignments and covariates for the econometric evaluation of training programs in LaLonde (1986).}
\usage{LaLonde}
\format{A data frame consisting of 12 columns (including a treatment assignment vector) and 3212 observations.  }
\source{Data from the National Supported Work Study.  A benchmark matching dataset.  Columns consist of an indicator for whether the observed unit was in the experimental subset; an indicator for whether the individual received the treatment; age in years; schooling in years; indicators for black and Hispanic; an indicator for marriage status, one of married; an indicator for no high school degree; reported earnings in 1974, 1975, and 1978; and whether the 1974 earnings variable is missing.  Data not missing 1974 earnings are the Dehejia-Wahba subsample of the LaLonde data.  Missing values for 1974 earnings set to zero. 1974 and 1975 earnings are pre-treatment.  1978 earnings is taken as the outcome variable.}
\references{LaLonde, R.J. (1986). Evaluating the econometric evaluations of training programs with experimental data. American Economic Review 76, 4, 604-620.}
\keyword{datasets}