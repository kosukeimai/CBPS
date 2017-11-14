

#' Blackwell Data for Covariate Balancing Propensity Score
#' 
#' This data set gives the outcomes a well as treatment assignments and
#' covariates for the example from Blackwell (2013).
#' 
#' 
#' @name Blackwell
#' @docType data
#' @format A data frame consisting of 13 columns (including treatment
#' assignment, time, and identifier vectors) and 570 observations.
#' @references Blackwell, Matthew. (2013). A framework for dynamic causal
#' inference in political science. American Journal of Political Science 57, 2,
#' 504-619.
#' @source d.gone.neg is the treatment. d.gone.neg.l1, d.gone.neg.l2, and
#' d.gone.neg.l3 are lagged treatment variables. camp.length, deminc,
#' base.poll, base.und, and office covariates. year is the year of the
#' particular race, and time goes from the first measurement (time = 1) to the
#' election (time = 5). demName is the identifier, and demprcnt is the outcome.
#' @keywords datasets
NULL





#' LaLonde Data for Covariate Balancing Propensity Score
#' 
#' This data set gives the outcomes a well as treatment assignments and
#' covariates for the econometric evaluation of training programs in LaLonde
#' (1986).
#' 
#' 
#' @name LaLonde
#' @docType data
#' @format A data frame consisting of 12 columns (including a treatment
#' assignment vector) and 3212 observations.
#' @references LaLonde, R.J. (1986). Evaluating the econometric evaluations of
#' training programs with experimental data. American Economic Review 76, 4,
#' 604-620.
#' @source Data from the National Supported Work Study.  A benchmark matching
#' dataset.  Columns consist of an indicator for whether the observed unit was
#' in the experimental subset; an indicator for whether the individual received
#' the treatment; age in years; schooling in years; indicators for black and
#' Hispanic; an indicator for marriage status, one of married; an indicator for
#' no high school degree; reported earnings in 1974, 1975, and 1978; and
#' whether the 1974 earnings variable is missing.  Data not missing 1974
#' earnings are the Dehejia-Wahba subsample of the LaLonde data.  Missing
#' values for 1974 earnings set to zero. 1974 and 1975 earnings are
#' pre-treatment.  1978 earnings is taken as the outcome variable.
#' @keywords datasets
NULL



