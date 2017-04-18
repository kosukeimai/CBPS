# CBPS: Covariate Balancing Propensity Score [![Build Status](https://travis-ci.org/kosukeimai/CBPS.svg?branch=master)](https://travis-ci.org/kosukeimai/CBPS) [![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/CBPS)](https://cran.r-project.org/package=CBPS)


Implements the covariate balancing propensity score (CBPS) proposed by
[Imai and Ratkovic (2014)](https://doi.org/10.1111/rssb.12027). The propensity
score is estimated such that it maximizes the resulting covariate
balance as well as the prediction of treatment assignment. The method,
therefore, avoids an iteration between model fitting and balance
checking.  The package also implements several extensions of the CBPS
beyond the cross-sectional, binary treatment setting.  The current
version implements the CBPS for longitudinal settings so that it can
be used in conjunction with marginal structural models from [Imai and
Ratkovic (2015)](https://doi.org/10.1080/01621459.2014.956872), treatments with
three- and four-valued treatment variables, continuous-valued
treatments from [Fong, Hazlett, and Imai (2015)]
(http://imai.princeton.edu/research/files/CBGPS.pdf), and the
situation with multiple distinct binary treatments administered
simultaneously. In the future it will be extended to other settings
including the generalization of experimental and instrumental variable
estimates. Recently add the optimal CBPS which chooses the optimal
balancing function and results in doubly robust and efficient
estimator for the treatment effect.
