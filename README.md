# CBPS: Covariate Balancing Propensity Score [![Build Status](https://travis-ci.org/kosukeimai/CBPS.svg?branch=master)](https://travis-ci.org/kosukeimai/CBPS) [![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/CBPS)](https://cran.r-project.org/package=CBPS) ![CRAN downloads](http://cranlogs.r-pkg.org/badges/grand-total/CBPS)

Implements the covariate balancing propensity score (CBPS) proposed by [Imai and Ratkovic (2014)](https://doi.org/10.1111/rssb.12027). The propensity score is estimated such that it maximizes the resulting covariate balance as well as the prediction of treatment assignment. The method, therefore, avoids an iteration between model fitting and balance checking.  The package also implements optimal CBPS from [Fan et al. (in-press)](https://doi.org/10.1080/07350015.2021.2002159),   several extensions of the CBPS beyond the cross-sectional, binary treatment setting.  They include the CBPS for longitudinal settings so that it can be used in conjunction with marginal structural models from [Imai and Ratkovic (2015)](https://doi.org/10.1080/01621459.2014.956872), treatments with three- and four-valued treatment variables, continuous-valued treatments from [Fong, Hazlett, and Imai (2018)](https://doi.org/10.1214/17-AOAS1101), propensity score estimation with a large number of covariates from [Ning, Peng, and Imai (2020)](https://doi.org/10.1093/biomet/asaa020), and the situation with multiple distinct binary treatments administered simultaneously. In the future it will be extended to other settings including the generalization of experimental and instrumental variable estimates.

## Installation
```
## from CRAN
install.packages("CBPS")

## from github
library(devtools)
install_github("kosukeimai/CBPS", dependencies = TRUE)
```
