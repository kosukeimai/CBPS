#' @title Covariate Balancing Propensity Score (CBPS) Estimation
#' 
#' @description
#' \code{CBPS} estimates propensity scores such that both covariate balance and
#' prediction of treatment assignment are maximized.  The method, therefore,
#' avoids an iterative process between model fitting and balance checking and
#' implements both simultaneously. For cross-sectional data, the method can
#' take continuous treatments and treatments with a control (baseline)
#' condition and either 1, 2, or 3 distinct treatment conditions.
#' 
#' Fits covariate balancing propensity scores.
#' 
#' ### @aliases CBPS CBPS.fit print.CBPS
#' 
#' @importFrom MASS mvrnorm ginv
#' @importFrom nnet multinom
#' @importFrom numDeriv jacobian
#' @importFrom MatchIt matchit
#' @importFrom glmnet cv.glmnet
#' @importFrom graphics abline axis layout mtext par plot points
#' @importFrom stats .getXlevels as.formula binomial coef cor dnorm glm is.empty.model lm model.frame 
#' @importFrom stats model.matrix model.response naprint optim optimize pnorm predict sd symnum var terms
#' @importFrom utils packageDescription
#' 
#' @param formula An object of class \code{formula} (or one that can be coerced
#' to that class): a symbolic description of the model to be fitted.
#' @param data An optional data frame, list or environment (or object coercible
#' by as.data.frame to a data frame) containing the variables in the model. If
#' not found in data, the variables are taken from \code{environment(formula)},
#' typically the environment from which \code{CBPS} is called.
#' @param na.action A function which indicates what should happen when the data
#' contain NAs. The default is set by the na.action setting of options, and is
#' na.fail if that is unset.
#' @param ATT Default is 1, which finds the average treatment effect on the
#' treated interpreting the second level of the treatment factor as the
#' treatment.  Set to 2 to find the ATT interpreting the first level of the
#' treatment factor as the treatment.  Set to 0 to find the average treatment
#' effect. For non-binary treatments, only the ATE is available.
#' @param iterations An optional parameter for the maximum number of iterations
#' for the optimization.  Default is 1000.
#' @param standardize Default is \code{TRUE}, which normalizes weights to sum
#' to 1 within each treatment group.  For continuous treatments, normalizes
#' weights to sum up to 1 for the entire sample.  Set to \code{FALSE} to return
#' Horvitz-Thompson weights.
#' @param method Choose "over" to fit an over-identified model that combines
#' the propensity score and covariate balancing conditions; choose "exact" to
#' fit a model that only contains the covariate balancing conditions.
#' @param twostep Default is \code{TRUE} for a two-step estimator, which will
#' run substantially faster than continuous-updating.  Set to \code{FALSE} to
#' use the continuous-updating estimator described by Imai and Ratkovic (2014).
#' @param sample.weights Survey sampling weights for the observations, if
#' applicable.  When left NULL, defaults to a sampling weight of 1 for each
#' observation.
#' @param baseline.formula Used only to fit iCBPS (see Fan et al). Currently
#' only works with binary treatments.  A formula specifying the balancing
#' covariates in the baseline outcome model, i.e., E(Y(0)|X).
#' @param diff.formula Used only to fit iCBPS (see Fan et al).  Currently only
#' works with binary treatments.  A formula specifying the balancing covariates
#' in the difference between the treatment and baseline outcome model, i.e.,
#' E(Y(1)-Y(0)|X).
#' @param ... Other parameters to be passed through to \code{optim()}.
#' 
#' @return \item{fitted.values}{The fitted propensity score}
#' \item{linear.predictor}{X * beta}
#' 
#' \item{deviance}{Minus twice the log-likelihood of the CBPS fit}
#' \item{weights}{The optimal weights.  Let \eqn{\pi_i = f(T_i | X_i)}{\pi_i =
#' f(T_i | X_i)}.  For binary ATE, these are given by \eqn{\frac{T_i}{\pi_i} +
#' \frac{(1 - T_i)}{(1 - \pi_i)}}{T_i/\pi_i + (1 - T_i)/(1 - \pi_i)}.  For
#' binary ATT, these are given by \eqn{\frac{n}{n_t} * \frac{T_i - \pi_i}{1 -
#' \pi_i}}{n/n_t * (T_i - \pi_i)/(1 - \pi_i)}.  For multi_valued treatments,
#' these are given by \eqn{\sum_{j=0}^{J-1} T_{i,j} /
#' \pi_{i,j}}{\sum_{j=0}^{J-1} T_i,j / \pi_i,j}.  For continuous treatments,
#' these are given by \eqn{\frac{f(T_i)}{f(T_i | X_i)}}{f(T_i) / f(T_i | X_i)
#' }.  These expressions for weights are all before standardization (i.e. with
#' standardize=\code{FALSE}).  Standardization will make weights sum to 1
#' within each treatment group.  For continuous treatment, standardization will
#' make all weights sum to 1.  If sampling weights are used, the weight for
#' each observation is multiplied by the survey sampling weight.} \item{y}{The
#' treatment vector used} \item{x}{The covariate matrix} \item{model}{The model
#' frame} \item{converged}{Convergence value.  Returned from the call to
#' \code{optim()}.} \item{call}{The matched call} \item{formula}{The formula
#' supplied} \item{data}{The data argument} \item{coefficients}{A named vector
#' of coefficients} \item{sigmasq}{The sigma-squared value, for continuous
#' treatments only} \item{J}{The J-statistic at convergence} \item{mle.J}{The
#' J-statistic for the parameters from maximum likelihood estimation}
#' \item{var}{The covariance matrix for the coefficients.} \item{Ttilde}{For
#' internal use only.} \item{Xtilde}{For internal use only.}
#' \item{beta.tilde}{For internal use only.} \item{simgasq.tilde}{For internal
#' use only.}
#' @author Christian Fong, Marc Ratkovic, Kosuke Imai, and Xiaolin Yang; The
#' CBPS function is based on the code for version 2.15.0 of the glm function
#' implemented in the stats package, originally written by Simon Davies.  This
#' documentation is likewise modeled on the documentation for glm and borrows
#' its language where the arguments and values are the same.
#' @seealso \link{summary.CBPS}
#' @references Imai, Kosuke and Marc Ratkovic.  2014. ``Covariate Balancing
#' Propensity Score.'' Journal of the Royal Statistical Society, Series B
#' (Statistical Methodology).
#' \url{http://imai.princeton.edu/research/CBPS.html} \cr Fong, Christian, Chad
#' Hazlett, and Kosuke Imai.  ``Parametric and Nonparametric Covariate
#' Balancing Propensity Score for General Treatment Regimes.'' Unpublished
#' Manuscript. \url{http://imai.princeton.edu/research/files/CBGPS.pdf} \cr
#' Fan, Jianqing and Imai, Kosuke and Liu, Han and Ning, Yang and Yang,
#' Xiaolin. ``Improving Covariate Balancing Propensity Score: A Doubly Robust
#' and Efficient Approach.'' Unpublished Manuscript.
#' \url{http://imai.princeton.edu/research/CBPStheory.html}
#' @examples
#' 
#' ###
#' ### Example: propensity score matching
#' ###
#' 
#' ##Load the LaLonde data
#' data(LaLonde)
#' ## Estimate CBPS
#' fit <- CBPS(treat ~ age + educ + re75 + re74 + 
#' 			I(re75==0) + I(re74==0), 
#' 			data = LaLonde, ATT = TRUE)
#' summary(fit)
#' \dontrun{
#' ## matching via MatchIt: one to one nearest neighbor with replacement
#' library(MatchIt)
#' m.out <- matchit(treat ~ fitted(fit), method = "nearest", 
#' 				 data = LaLonde, replace = TRUE)
#' 
#' ### Example: propensity score weighting 
#' ###
#' ## Simulation from Kang and Shafer (2007).
#' set.seed(123456)
#' n <- 500
#' X <- mvrnorm(n, mu = rep(0, 4), Sigma = diag(4))
#' prop <- 1 / (1 + exp(X[,1] - 0.5 * X[,2] + 
#' 			 0.25*X[,3] + 0.1 * X[,4]))
#' treat <- rbinom(n, 1, prop)
#' y <- 210 + 27.4*X[,1] + 13.7*X[,2] + 13.7*X[,3] + 13.7*X[,4] + rnorm(n)
#' 
#' ##Estimate CBPS with a misspecified model
#' X.mis <- cbind(exp(X[,1]/2), X[,2]*(1+exp(X[,1]))^(-1)+10, 
#' 			  (X[,1]*X[,3]/25+.6)^3, (X[,2]+X[,4]+20)^2)
#' fit1 <- CBPS(treat ~ X.mis, ATT = 0)
#' summary(fit1)
#' 	
#' ## Horwitz-Thompson estimate
#' mean(treat*y/fit1$fitted.values)
#' ## Inverse propensity score weighting
#' sum(treat*y/fit1$fitted.values)/sum(treat/fit1$fitted.values)
#' 
#' rm(list=c("y","X","prop","treat","n","X.mis","fit1"))
#' 
#' ### Example: Continuous Treatment
#' set.seed(123456)
#' n <- 1000
#' X <- mvrnorm(n, mu = rep(0,2), Sigma = diag(2))
#' beta <- rnorm(ncol(X)+1, sd = 1)
#' treat <- cbind(1,X)%*%beta + rnorm(n, sd = 5)
#' 
#' treat.effect <- 1
#' effect.beta <- rnorm(ncol(X))
#' y <- rbinom(n, 1, (1 + exp(-treat.effect*treat - 
#' 				   X%*%effect.beta))^-1)
#' 
#' fit2 <- CBPS(treat ~ X)
#' summary(fit2)
#' summary(glm(y ~ treat + X, weights = fit2$weights, 
#' 			family = "quasibinomial"))
#' 
#' rm(list=c("n", "X", "beta", "treat", "treat.effect",
#' 		  "effect.beta", "y", "fit2"))
#' 
#' ### Simulation example: Improved CBPS (or iCBPS) from Fan et al
#' set.seed(123456)
#' n <- 500
#' X <- mvrnorm(n, mu = rep(0, 4), Sigma = diag(4))
#' prop <- 1 / (1 + exp(X[,1] - 0.5 * X[,2] + 0.25*X[,3] + 0.1 * X[,4]))
#' treat <- rbinom(n, 1, prop)
#' y1 <- 210 + 27.4*X[,1] + 13.7*X[,2] + 13.7*X[,3] + 13.7*X[,4] + rnorm(n)
#' y0 <- 210 + 13.7*X[,2] + 13.7*X[,3] + 13.7*X[,4] + rnorm(n)
#' ##Estimate iCBPS with a misspecificied model
#' X.mis <- cbind(exp(X[,1]/2), X[,2]*(1+exp(X[,1]))^(-1)+10, 
#' 			   (X[,1]*X[,3]/25+.6)^3, (X[,2]+X[,4]+20)^2)
#' fit1 <- CBPS(treat ~ X.mis, baseline.formula=~X.mis[,2:4], 
#' 			 diff.formula=~X.mis[,1], ATT = FALSE)
#' summary(fit1)
#' }
#' 
#' @export CBPS
#' 
CBPS <- function(formula, data, na.action, ATT=1, iterations=1000, standardize=TRUE, method="over", twostep=TRUE,
                 sample.weights=NULL, baseline.formula=NULL, diff.formula=NULL,...) {
  if (missing(data)) 
    data <- environment(formula)
  call <- match.call()
  family <- binomial()
  
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  Y <- model.response(mf, "any")
  if (length(dim(Y)) == 1L) {
    nm <- rownames(Y)
    dim(Y) <- NULL
    if (!is.null(nm)) 
      names(Y) <- nm
  }
  
  X <- if (!is.empty.model(mt)) model.matrix(mt, mf)#[,-2]
  else matrix(NA, NROW(Y), 0L)
  
  X<-cbind(1,X[,apply(X,2,sd)>0])
  
  #Handle sample weights
  if(is.null(sample.weights)) sample.weights<-rep(1,nrow(X))
  
  # Parse formulae 2 and 3, if they are necessary
  if (xor(is.null(baseline.formula), is.null(diff.formula))){
    stop("Either baseline.formula or diff.formula not specified.  Both must be specified to use CBPSOptimal.  Otherwise, leave both NULL.")
  }
  if(!is.null(baseline.formula))
  {	
    baselineX<-model.matrix(terms(baseline.formula))
    baselineX<-baselineX[,apply(baselineX,2,sd)>0]
    diffX<- model.matrix(terms(diff.formula))
    diffX<-diffX[,apply(as.matrix(diffX),2,sd)>0]
  }
  else
  {
    baselineX <- NULL
    diffX <- NULL
  }
  
  fit <- eval(call("CBPS.fit", X = X, treat = Y, ATT=ATT, 
                   intercept = attr(mt, "intercept") > 0L, method=method, iterations=iterations, 
                   standardize = standardize, twostep = twostep, 
                   baselineX = baselineX, diffX = diffX,sample.weights=sample.weights))	
  
  fit$na.action <- attr(mf, "na.action")
  xlevels <- .getXlevels(mt, mf)
  fit$data<-data
  fit$call <- call
  fit$formula <- formula
  fit$terms<-mt
  fit
}

#' CBPS.fit determines the proper routine (what kind of treatment) and calls the 
#' approporiate function.  It also pre- and post-processes the data
#' 
#'
#' @param ATT Default is 1, which finds the average treatment effect on the
#' treated interpreting the second level of the treatment factor as the
#' treatment.  Set to 2 to find the ATT interpreting the first level of the
#' treatment factor as the treatment.  Set to 0 to find the average treatment
#' effect. For non-binary treatments, only the ATE is available.
#' @param iterations An optional parameter for the maximum number of iterations
#' for the optimization.  Default is 1000.
#' @param standardize Default is \code{TRUE}, which normalizes weights to sum
#' to 1 within each treatment group.  For continuous treatments, normalizes
#' weights to sum up to 1 for the entire sample.  Set to \code{FALSE} to return
#' Horvitz-Thompson weights.
#' @param method Choose "over" to fit an over-identified model that combines
#' the propensity score and covariate balancing conditions; choose "exact" to
#' fit a model that only contains the covariate balancing conditions.
#' @param twostep Default is \code{TRUE} for a two-step estimator, which will
#' run substantially faster than continuous-updating.  Set to \code{FALSE} to
#' use the continuous-updating estimator described by Imai and Ratkovic (2014).
#' @param treat A vector of treatment assignments.  Binary or multi-valued
#' treatments should be factors.  Continuous treatments should be numeric.
#' @param X A covariate matrix.
#' @param sample.weights Survey sampling weights for the observations, if
#' applicable.  When left NULL, defaults to a sampling weight of 1 for each
#' observation.
#' @param baselineX Similar to \code{baseline.formula}, but in matrix form.
#' @param diffX Similar to \code{diff.formula}, but in matrix form.
#' @param ... Other parameters to be passed through to \code{optim()}.
#' 
#' @return CBPS.fit object
#'
#' @export
#' 
CBPS.fit<-function(treat, X, baselineX, diffX, ATT, method, iterations, standardize, twostep, sample.weights=sample.weights,...){
  # Special clause interprets T = 1 or 0 as a binary treatment, even if it is numeric
  if ((levels(factor(treat))[1] %in% c("FALSE","0",0)) & (levels(factor(treat))[2] %in% c("TRUE","1",1))
      & (length(levels(factor(treat))) == 2))
  {
    treat<-factor(treat)
  }
  
  # Declare some constants and orthogonalize Xdf.
  
  k=0
  if(method=="over") bal.only=FALSE
  if(method=="exact") bal.only=TRUE

  names.X<-colnames(X)
  names.X[apply(X,2,sd)==0]<-"(Intercept)"

  # Only preprocess if not doing CBPS Optimal
  X.orig<-X
  if(is.null(baselineX)){
    x.sd<-apply(as.matrix(X[,-1]),2,sd)
    Dx.inv<-diag(c(1,x.sd))
    x.mean<-apply(as.matrix(X[,-1]),2,mean)
    X[,-1]<-apply(as.matrix(X[,-1]),2,FUN=function(x) (x-mean(x))/sd(x))
    svd1<-svd(X)
    X<-svd1$u   
  }
  k<-qr(X)$rank
  if (k < ncol(X)) stop("X is not full rank")
  
  # When you take the svd, this is the identity matrix.  Perhaps
  # we forgot to work this in somewhere
  XprimeX.inv<-ginv(t(sample.weights^.5*X)%*%(sample.weights^.5*X))
  
  # Determine the number of treatments
  if (is.factor(treat)) {
    no.treats<-length(levels(treat))
    if (no.treats > 4) stop("Parametric CBPS is not implemented for more than 4 treatment values.  Consider using a continuous value.")
    if (no.treats < 2) stop("Treatment must take more than one value")
    
    if (no.treats == 2)
    {
      if (!is.null(baselineX) && !is.null(diffX))
      {
        if(ATT==1)
        {
          message("CBPSOptimal does not support ATT=1 for now. Try ATT=0.")
        }
        output<-CBPSOptimal.2Treat(treat, X, baselineX, diffX, iterations, ATT=0, standardize = standardize)
      }
      else
      {
        output<-CBPS.2Treat(treat, X, method, k, XprimeX.inv, bal.only, iterations, ATT, standardize = standardize, twostep = twostep, sample.weights=sample.weights)
      }
    }
    
    if (no.treats == 3)
    {
      output<-CBPS.3Treat(treat, X, method, k, XprimeX.inv, bal.only, iterations, standardize = standardize, twostep = twostep, sample.weights=sample.weights)
    }
    
    if (no.treats == 4)
    {
      output<-CBPS.4Treat(treat, X, method, k, XprimeX.inv, bal.only, iterations, standardize = standardize, twostep = twostep, sample.weights=sample.weights)
    }
    
    # Reverse the svd, centering and scaling
    if (is.null(baselineX)){
      d.inv<- svd1$d
      d.inv[d.inv> 1e-5]<-1/d.inv[d.inv> 1e-5]
      d.inv[d.inv<= 1e-5]<-0      
      beta.opt<-svd1$v%*%diag(d.inv)%*%coef(output)
      beta.opt[-1,]<-beta.opt[-1,]/x.sd
      beta.opt[1,]<-beta.opt[1,]-matrix(x.mean%*%beta.opt[-1,])
    }
    else{ #added by Xiaolin to deal with cases when baselineX is not null. (12/26/2016)
      beta.opt<-as.matrix(coef(output))
    }
    output$coefficients<-beta.opt
    output$x<-X.orig

    rownames(output$coefficients)<-names.X
    
    # Calculate the variance
    variance<-output$var
    
    if (no.treats == 2){
      colnames(output$coefficients)<-c("Treated")
      if (is.null(baselineX)){
         output$var<-Dx.inv%*%ginv(t(X.orig)%*%X.orig)%*%t(X.orig)%*%X%*%svd1$v%*%ginv(diag(svd1$d))%*%variance%*%ginv(diag(svd1$d))%*%t(svd1$v)%*%t(X)%*%X.orig%*%ginv(t(X.orig)%*%X.orig)%*%Dx.inv
      }
      else{
        output$var<-variance
      }
      colnames(output$var)<-names.X
      rownames(output$var)<-colnames(output$var)
    }
    
    if (no.treats == 3){
      colnames(output$coefficients)<-levels(as.factor(treat))[c(2,3)]
      var.1.1<-variance[1:k,1:k]
      var.1.2<-variance[1:k,(k+1):(2*k)]
      var.2.1<-variance[(k+1):(2*k),1:k]
      var.2.2<-variance[(k+1):(2*k),(k+1):(2*k)]
      trans.var.1.1<-Dx.inv%*%ginv(t(X.orig)%*%X.orig)%*%t(X.orig)%*%X%*%svd1$v%*%ginv(diag(svd1$d))%*%var.1.1%*%ginv(diag(svd1$d))%*%t(svd1$v)%*%t(X)%*%X.orig%*%ginv(t(X.orig)%*%X.orig)%*%Dx.inv
      trans.var.1.2<-Dx.inv%*%ginv(t(X.orig)%*%X.orig)%*%t(X.orig)%*%X%*%svd1$v%*%ginv(diag(svd1$d))%*%var.1.2%*%ginv(diag(svd1$d))%*%t(svd1$v)%*%t(X)%*%X.orig%*%ginv(t(X.orig)%*%X.orig)%*%Dx.inv
      trans.var.2.1<-Dx.inv%*%ginv(t(X.orig)%*%X.orig)%*%t(X.orig)%*%X%*%svd1$v%*%ginv(diag(svd1$d))%*%var.2.1%*%ginv(diag(svd1$d))%*%t(svd1$v)%*%t(X)%*%X.orig%*%ginv(t(X.orig)%*%X.orig)%*%Dx.inv
      trans.var.2.2<-Dx.inv%*%ginv(t(X.orig)%*%X.orig)%*%t(X.orig)%*%X%*%svd1$v%*%ginv(diag(svd1$d))%*%var.2.2%*%ginv(diag(svd1$d))%*%t(svd1$v)%*%t(X)%*%X.orig%*%ginv(t(X.orig)%*%X.orig)%*%Dx.inv
      output$var<-rbind(cbind(trans.var.1.1,trans.var.1.2),cbind(trans.var.2.1,trans.var.2.2))
      colnames(output$var)<-c(paste0(levels(as.factor(treat))[2],": ", names.X),paste0(levels(as.factor(treat))[3], ": ", names.X))
      rownames(output$var)<-colnames(output$var)
    }
    
    if (no.treats == 4)
    {
      colnames(output$coefficients)<-levels(as.factor(treat))[c(2,3,4)]
      var.1.1<-variance[1:k,1:k]
      var.1.2<-variance[1:k,(k+1):(2*k)]
      var.1.3<-variance[1:k,(2*k+1):(3*k)]
      var.2.1<-variance[(k+1):(2*k),1:k]
      var.2.2<-variance[(k+1):(2*k),(k+1):(2*k)]
      var.2.3<-variance[(k+1):(2*k),(2*k+1):(3*k)]
      var.3.1<-variance[(2*k+1):(3*k),1:k]
      var.3.2<-variance[(2*k+1):(3*k),(k+1):(2*k)]
      var.3.3<-variance[(2*k+1):(3*k),(2*k+1):(3*k)]
      trans.var.1.1<-Dx.inv%*%ginv(t(X.orig)%*%X.orig)%*%t(X.orig)%*%X%*%svd1$v%*%ginv(diag(svd1$d))%*%var.1.1%*%ginv(diag(svd1$d))%*%t(svd1$v)%*%t(X)%*%X.orig%*%ginv(t(X.orig)%*%X.orig)%*%Dx.inv
      trans.var.1.2<-Dx.inv%*%ginv(t(X.orig)%*%X.orig)%*%t(X.orig)%*%X%*%svd1$v%*%ginv(diag(svd1$d))%*%var.1.2%*%ginv(diag(svd1$d))%*%t(svd1$v)%*%t(X)%*%X.orig%*%ginv(t(X.orig)%*%X.orig)%*%Dx.inv
      trans.var.1.3<-Dx.inv%*%ginv(t(X.orig)%*%X.orig)%*%t(X.orig)%*%X%*%svd1$v%*%ginv(diag(svd1$d))%*%var.1.3%*%ginv(diag(svd1$d))%*%t(svd1$v)%*%t(X)%*%X.orig%*%ginv(t(X.orig)%*%X.orig)%*%Dx.inv
      trans.var.2.1<-Dx.inv%*%ginv(t(X.orig)%*%X.orig)%*%t(X.orig)%*%X%*%svd1$v%*%ginv(diag(svd1$d))%*%var.2.1%*%ginv(diag(svd1$d))%*%t(svd1$v)%*%t(X)%*%X.orig%*%ginv(t(X.orig)%*%X.orig)%*%Dx.inv
      trans.var.2.2<-Dx.inv%*%ginv(t(X.orig)%*%X.orig)%*%t(X.orig)%*%X%*%svd1$v%*%ginv(diag(svd1$d))%*%var.2.2%*%ginv(diag(svd1$d))%*%t(svd1$v)%*%t(X)%*%X.orig%*%ginv(t(X.orig)%*%X.orig)%*%Dx.inv
      trans.var.2.3<-Dx.inv%*%ginv(t(X.orig)%*%X.orig)%*%t(X.orig)%*%X%*%svd1$v%*%ginv(diag(svd1$d))%*%var.2.3%*%ginv(diag(svd1$d))%*%t(svd1$v)%*%t(X)%*%X.orig%*%ginv(t(X.orig)%*%X.orig)%*%Dx.inv
      trans.var.3.1<-Dx.inv%*%ginv(t(X.orig)%*%X.orig)%*%t(X.orig)%*%X%*%svd1$v%*%ginv(diag(svd1$d))%*%var.3.1%*%ginv(diag(svd1$d))%*%t(svd1$v)%*%t(X)%*%X.orig%*%ginv(t(X.orig)%*%X.orig)%*%Dx.inv
      trans.var.3.2<-Dx.inv%*%ginv(t(X.orig)%*%X.orig)%*%t(X.orig)%*%X%*%svd1$v%*%ginv(diag(svd1$d))%*%var.3.2%*%ginv(diag(svd1$d))%*%t(svd1$v)%*%t(X)%*%X.orig%*%ginv(t(X.orig)%*%X.orig)%*%Dx.inv
      trans.var.3.3<-Dx.inv%*%ginv(t(X.orig)%*%X.orig)%*%t(X.orig)%*%X%*%svd1$v%*%ginv(diag(svd1$d))%*%var.3.3%*%ginv(diag(svd1$d))%*%t(svd1$v)%*%t(X)%*%X.orig%*%ginv(t(X.orig)%*%X.orig)%*%Dx.inv
      output$var<-rbind(cbind(trans.var.1.1,trans.var.1.2,trans.var.1.3),cbind(trans.var.2.1,trans.var.2.2,trans.var.2.3),cbind(trans.var.3.1,trans.var.3.2,trans.var.3.3))
        
      colnames(output$var)<-c(paste0(levels(as.factor(treat))[2],": ", names.X),paste0(levels(as.factor(treat))[3], ": ", names.X),paste0(levels(as.factor(treat))[4], ": ", names.X))
      rownames(output$var)<-colnames(output$var)
    }
  } else if (is.numeric(treat)) {
    # Warn if it seems like the user meant to input a categorical treatment
    if (length(unique(treat)) <= 4) warning("Treatment vector is numeric.  Interpreting as a continuous treatment.  To solve for a binary or multi-valued treatment, make treat a factor.")
    output<-CBPS.Continuous(treat, X, method, k, XprimeX.inv, bal.only, iterations, standardize = standardize, twostep = twostep, sample.weights=sample.weights)
    
    # Reverse svd, centering, and scaling
    d.inv<- svd1$d
    d.inv[d.inv> 1e-5]<-1/d.inv[d.inv> 1e-5]
    d.inv[d.inv<= 1e-5]<-0
    beta.opt<-svd1$v%*%diag(d.inv)%*%coef(output)
    beta.opt[-1,]<-beta.opt[-1,]/x.sd
      
    beta.opt[1,]<-beta.opt[1,]-matrix(x.mean%*%beta.opt[-1,])
    output$coefficients<-as.matrix(beta.opt)
    output$x<-X.orig
      
    rownames(output$coefficients)<-c(names.X)
    
    # Calculate variance
    var.1<-output$var
    output$var<-Dx.inv%*%ginv(t(X.orig)%*%X.orig)%*%t(X.orig)%*%X%*%svd1$v%*%ginv(diag(svd1$d))%*%var.1%*%ginv(diag(svd1$d))%*%t(svd1$v)%*%t(X)%*%X.orig%*%ginv(t(X.orig)%*%X.orig)%*%Dx.inv
    rownames(output$var)<-names.X
    colnames(output$var)<-rownames(output$var)
  } else {
    stop("Treatment must be either a factor or numeric")
  }
  
  output$method<-method
  
  output
}

#' Print coefficients and model fit statistics
#' @param x an object of class \dQuote{CBPS} or \dQuote{npCBPS}, usually, a result of a call to \code{CBPS} or \code{npCBPS}.
#' @param digits the number of digits to keep for the numerical quantities.
#' @param ... Additional arguments to be passed to summary.
#'
#' @export
#' 
print.CBPS <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  if (length(coef(x))) {
    cat("Coefficients:\n")
    print.default(format(x$coefficients, digits = digits), 
                  print.gap = 2, quote = FALSE)
  }
  else cat("No coefficients\n\n")
  if (max(class(x) == "CBPScontinuous"))
    cat("\nSigma-Squared: ",x$sigmasq)
  if (nzchar(mess <- naprint(x$na.action))) 
    cat("  (", mess, ")\n", sep = "")
  cat("Residual Deviance:\t", format(signif(x$deviance, 
                                            digits)), "\n")
  cat("J-Statistic:\t		", format(signif(x$J)),"\n")
  cat("Log-Likelihood:\t ",-0.5*x$deviance, "\n")
  invisible(x)
}

# Expands on print by including uncertainty for coefficient estimates


#' Summarizing Covariate Balancing Propensity Score Estimation
#' 
#' Prints a summary of a fitted CBPS object.
#' 
#' Prints a summmary of a CBPS object, in a format similar to glm.  The
#' variance matrix is calculated from the numerical Hessian at convergence of
#' CBPS.
#' 
#' @param object an object of class \dQuote{CBPS}, usually, a result of a call
#' to CBPS.
#' @param ... Additional arguments to be passed to summary.
#' @return \item{call}{The matched call.} \item{deviance.residuals}{The five
#' number summary and the mean of the deviance residuals.}
#' \item{coefficients}{A table including the estimate for the each coefficient
#' and the standard error, z-value, and two-sided p-value for these estimates.}
#' \item{J}{Hansen's J-Statistic for the fitted model.}
#' \item{Log-Likelihood}{The log-likelihood of the fitted model.}
#' @author Christian Fong, Marc Ratkovic, and Kosuke Imai.
#' @seealso \link{CBPS}, \link{summary}
#'
#' @export
#' 
summary.CBPS<-function(object, ...){
  ##x <- summary.glm(object, dispersion = dispersion, correlation = correlation, symbolic.cor = symbolic.cor, ...)	  
  x<-NULL
  names.X<-as.vector(names(object$coefficients))
  sd.coef <- diag(object$var)^.5
  coef.table<-(cbind(as.vector(object$coefficients),as.vector(sd.coef),as.vector(object$coefficients/sd.coef),as.vector(2-2*pnorm(abs(object$coefficients/sd.coef)))))
  colnames(coef.table)<-c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  if (ncol(coef(object)) == 1)
  {
    rownames(coef.table)<-rownames(object$coefficients)#names.X
  }
  if (ncol(coef(object)) > 1)
  {
    rnames<-array()
    for (i in 1:ncol(coef(object)))
    {
      rnames[((i-1)*nrow(coef(object))+1):(i*nrow(coef(object)))]<-paste0(levels(as.factor(object$y))[i],": ",rownames(coef(object)))
    }
    rownames(coef.table)<-rnames
  }
  
  pval <- coef.table[,4]
  symp <- symnum(pval, corr=FALSE,
                 cutpoints = c(0,  .001,.01,.05, .1, 1),
                 symbols = c("***","**","*","."," "))
  coef.print<-cbind(signif(coef.table,3),as.vector(symp))
  coef.print[coef.print=="0"]<-"0.000"
  
  cat("\nCall:	\n", paste(deparse(object$call), sep = "\n", collapse = "\n"), 
      "\n", sep = "")
  
  
  cat("\nCoefficients:\n")
  
  print(noquote(coef.print))
  cat("---\n")
  cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 \n")
  #cat("\n	Null J:	 ",object$J)
  if(max(class(object)=="CBPScontinuous")){
    cat("\nSigma-Squared: ",object$sigmasq)
  }
  cat("\nJ - statistic:	 ",object$J)
  cat("\nLog-Likelihood: ",-0.5*object$deviance, "\n")
  
  out<-list("call"=object$call,"coefficients"=coef.table,"J"=object$J)
  invisible(out)
}




#' Calculate Variance-Covariance Matrix for a Fitted CBPS Object
#' 
#' \code{vcov.CBPS} Returns the variance-covariance matrix of the main
#' parameters of a fitted CBPS object.
#' 
#' This is the CBPS implementation of the generic function vcov().
#' 
#' @param object An object of class \code{formula} (or one that can be coerced
#' to that class): a symbolic description of the model to be fitted.
#' @param ... Additional arguments to be passed to vcov.CBPS
#' @return A matrix of the estimated covariances between the parameter
#' estimates in the linear or non-linear predictor of the model.
#' @author Christian Fong, Marc Ratkovic, and Kosuke Imai.
#' @seealso \link{vcov}
#' @references This documentation is modeled on the documentation of the
#' generic \link{vcov}.
#' @examples
#' 
#' ###
#' ### Example: Variance-Covariance Matrix
#' ###
#' 
#' ##Load the LaLonde data
#' data(LaLonde)
#' ## Estimate CBPS via logistic regression
#' fit <- CBPS(treat ~ age + educ + re75 + re74 + I(re75==0) + I(re74==0), 
#' 		    data = LaLonde, ATT = TRUE)
#' ## Get the variance-covariance matrix.
#' vcov(fit)
#' 
#'
#' @export
#' 
vcov.CBPS<-function(object,...){
  return(object$var)
}

# Plot binary and multi-valued CBPS.  Plots the standardized difference in means for each contrast
# before and after weighting.  Defined for an arbitrary number of discrete treatments.


#' Plotting Covariate Balancing Propensity Score Estimation
#' 
#' 
#' This function plots the absolute difference in standardized means before and after
#' weighting.  To access more sophisticated graphics for assessing covariate balance,
#' consider using Noah Greifer's \code{cobalt} package.
#' 
#' The "Before Weighting" plot gives the balance before weighting, and the
#' "After Weighting" plot gives the balance after weighting.
#' 
#' ### @aliases plot.CBPS plot.npCBPS
#' @param x an object of class \dQuote{CBPS} or \dQuote{npCBPS}, usually, a
#' result of a call to \code{CBPS} or \code{npCBPS}.
#' @param covars Indices of the covariates to be plotted (excluding the
#' intercept).  For example, if only the first two covariates from
#' \code{balance} are desired, set \code{covars} to 1:2.  The default is
#' \code{NULL}, which plots all covariates.
#' @param silent If set to \code{FALSE}, returns the imbalances used to
#' construct the plot.  Default is \code{TRUE}, which returns nothing.
#' @param boxplot If set to \code{TRUE}, returns a boxplot summarizing the
#' imbalance on the covariates instead of a point for each covariate.  Useful
#' if there are many covariates.
#' @param ... Additional arguments to be passed to plot.
#' @return For binary and multi-valued treatments, plots the absolute
#' difference in standardized means by contrast for all covariates before and
#' after weighting.  This quantity for a single covariate and a given pair of
#' treatment conditions is given by \eqn{\frac{\sum_{i=1}^{n} w_i * (T_i == 1)
#' * X_i}{\sum_{i=1}^{n} (T_i == 1) * w_i} - \frac{\sum_{i=1}^{n} w_i * (T_i ==
#' 0) * X_i}{\sum_{i=1}^{n} (T_i == 0) * w_i}}{[\sum w_i * (T_i == 1) *
#' X_i]/[\sum w_i * (T_i == 1)] - [\sum w_i * (T_i == 0) * X_i]/[\sum w_i *
#' (T_i == 0)]}.  For continuous treatments, plots the weighted absolute
#' Pearson correlation between the treatment and each covariate.  See
#' \url{https://en.wikipedia.org/wiki/Pearson_product-moment_correlation_coefficient#Weighted_correlation_coefficient.
#' }
#' @author Christian Fong, Marc Ratkovic, and Kosuke Imai.
#' @seealso \link{CBPS}, \link{plot}
#'
#' @export
#' 
plot.CBPS<-function(x, covars = NULL, silent = TRUE, boxplot = FALSE, ...){ 
  bal.x<-balance(x)
  if(is.null(covars))
  {
    covars<-1:nrow(bal.x[["balanced"]])
  }
  
  no.treats<-length(levels(as.factor(x$y)))	
  balanced.std.mean<-bal.x[["balanced"]][covars,]
  original.std.mean<-bal.x[["original"]][covars,]
  no.contrasts<-ifelse(no.treats == 2, 1, ifelse(no.treats == 3, 3, 6))
  
  abs.mean.ori.contrasts<-matrix(rep(0,no.contrasts*length(covars)),length(covars),no.contrasts)
  abs.mean.bal.contrasts<-matrix(rep(0,no.contrasts*length(covars)),length(covars),no.contrasts)
  contrast.names<-array()
  true.contrast.names<-array()
  contrasts<-c()
  covarlist<-c()
  ctr<-1
  for (i in 1:(no.treats-1))
  {
    for (j in (i+1):no.treats)
    {
      abs.mean.ori.contrasts[,ctr]<-abs(original.std.mean[covars,i+no.treats]-original.std.mean[covars,j+no.treats])
      abs.mean.bal.contrasts[,ctr]<-abs(balanced.std.mean[covars,i+no.treats]-balanced.std.mean[covars,j+no.treats])
      contrast.names[ctr]<-paste0(i,":",j)
      true.contrast.names[ctr]<-paste0(levels(as.factor(x$y))[i],":",levels(as.factor(x$y))[j])
      
      contrasts<-c(contrasts, rep(true.contrast.names[ctr],length(covars)))
      covarlist<-c(covarlist, rownames(balanced.std.mean))
      ctr<-ctr+1
    }
  } 
  
  max.abs.contrast<-max(max(abs.mean.ori.contrasts),max(abs.mean.bal.contrasts))
  m <- matrix(c(1,1,1,2,2,2,3,3,3),nrow = 3,ncol = 3,byrow = TRUE)
  layout(mat = m,heights = c(0.4,0.4,0.3))
  
  par(mfrow=c(2,1))
  if (!boxplot){
    plot(1, type="n", xlim=c(0,max.abs.contrast), ylim=c(1,no.contrasts),xlab="",ylab="",main="",yaxt='n', ...)
    axis(side=2, at=seq(1,no.contrasts),contrast.names)
    mtext("Absolute Difference of Standardized Means",side=1,line=2.25)
    mtext("Contrasts",side=2,line=2)
    mtext("Before Weighting",side=3,line=0.5,font=2)
    for (i in 1:no.contrasts)
    {
      for (j in 1:length(covars))
      {
        points(abs.mean.ori.contrasts[covars[j],i],i, ...)
      }
    }
    plot(1, type="n", xlim=c(0,max.abs.contrast), ylim=c(1,no.contrasts), xlab="", ylab="", main="", yaxt='n', ...)
    axis(side=2, at=seq(1,no.contrasts),contrast.names)
    mtext("Absolute Difference of Standardized Means",side=1,line=2.25)
    mtext("Contrasts",side=2,line=2)
    mtext("After Weighting",side=3,line=0.5,font=2)
    
    for (i in 1:no.contrasts)
    {
      for (j in 1:length(covars))
      {
        points(abs.mean.bal.contrasts[covars[j],i],i, ...)
      }
    }
  }
  else{
    boxplot(abs.mean.ori.contrasts, horizontal = TRUE, ylim = c(0,max.abs.contrast), xlim=c(1-0.5,no.contrasts+0.5), 
            xlab="", ylab="", main="", yaxt='n', ...)
    axis(side=2, at=seq(1,no.contrasts),contrast.names)
    mtext("Absolute Difference of Standardized Means",side=1,line=2.25)
    mtext("Contrasts",side=2,line=2)
    mtext("Before Weighting",side=3,line=0.5,font=2)
    
    boxplot(abs.mean.bal.contrasts, horizontal = TRUE, ylim = c(0,max.abs.contrast), xlim=c(1-0.5,no.contrasts+0.5), 
            xlab="", ylab="", main="", yaxt='n', ...)
    axis(side=2, at=seq(1,no.contrasts),contrast.names)
    mtext("Absolute Difference of Standardized Means",side=1,line=2.25)
    mtext("Contrasts",side=2,line=2)
    mtext("After Weighting",side=3,line=0.5,font=2)
  }
  
  par(mfrow=c(1,1))
  
  if(is.null(rownames(balanced.std.mean))) rownames(balanced.std.mean)<-paste0("X",covars)
  if(!silent) return(data.frame("contrast" = contrasts, "covariate" = covarlist, 
                                "balanced"=abs.mean.bal.contrasts,  
                                "original"=abs.mean.ori.contrasts))
}

#' Plot the pre-and-post weighting correlations between X and T
#' @param x an object of class \dQuote{CBPS} or \dQuote{npCBPS}, usually, a
#' result of a call to \code{CBPS} or \code{npCBPS}.
#' @param covars Indices of the covariates to be plotted (excluding the intercept).  For example, 
#' if only the first two covariates from \code{balance} are desired, set \code{covars} to 1:2.  
#' The default is \code{NULL}, which plots all covariates.
#' @param silent If set to \code{FALSE}, returns the imbalances used to
#' construct the plot.  Default is \code{TRUE}, which returns nothing.
#' @param boxplot If set to \code{TRUE}, returns a boxplot summarizing the
#' imbalance on the covariates instead of a point for each covariate.  Useful
#' if there are many covariates.
#' @param ... Additional arguments to be passed to balance.
#'
#' @export
#' 
plot.CBPSContinuous<-function(x, covars = NULL, silent = TRUE, boxplot = FALSE, ...){     
  bal.x<-balance(x)
  if (is.null(covars))
  {
    covars<-1:nrow(bal.x[["balanced"]])
  }
  balanced.abs.cor<-abs(bal.x[["balanced"]][covars])
  original.abs.cor<-abs(bal.x[["unweighted"]][covars])
  
  max.abs.cor<-max(max(original.abs.cor),max(balanced.abs.cor))
  
  if(!boxplot){
    plot(1, type="n", xlim=c(0,max.abs.cor), ylim=c(1.5,3.5), xlab = "Absolute Pearson Correlation", ylab = "", yaxt = "n", ...)
    axis(side=2, at=seq(2,3),c("CBPS Weighted", "Unweighted"))
    points(x=original.abs.cor, y=rep(3, length(covars)), pch=19)
    points(x=balanced.abs.cor, y=rep(2, length(covars)), pch=19)
  }
  else{
    boxplot(balanced.abs.cor, original.abs.cor, horizontal = TRUE, yaxt = 'n', xlab = "Absolute Pearson Correlation", ...)
    axis(side=2, at=c(1,2),c("CBPS Weighted", "Unweighted"))
  }
  
  if(!silent) return(data.frame("covariate"=rownames(bal.x[["balanced"]]),"balanced"=balanced.abs.cor,
                                "original"=original.abs.cor))
}



#' Optimal Covariate Balance
#' 
#' Returns the mean and standardized mean associated with each treatment group,
#' before and after weighting.  To access more comprehensive diagnotistics for 
#' assessing covariate balance, consider using Noah Greifer's \code{cobalt} package.
#' 
#' For binary and multi-valued treatments as well as marginal structural
#' models, each of the matrices' rows are the covariates and whose columns are
#' the weighted mean, and standardized mean associated with each treatment
#' group.  The standardized mean is the weighted mean divided by the standard
#' deviation of the covariate for the whole population.  For continuous
#' treatments, returns the absolute Pearson correlation between the treatment
#' and each covariate.
#' 
#' ### @aliases balance balance.npCBPS balance.CBPS balance.CBMSM
#' @param object A CBPS, npCBPS, or CBMSM object.
#' @param ... Additional arguments to be passed to balance.
#' @return Returns a list of two matrices, "original" (before weighting) and
#' "balanced" (after weighting).
#' @author Christian Fong, Marc Ratkovic, and Kosuke Imai.
#' @examples
#' 
#' ###
#' ### Example: Assess Covariate Balance
#' ###
#' data(LaLonde)
#' ## Estimate CBPS
#' fit <- CBPS(treat ~ age + educ + re75 + re74 + 
#' 			I(re75==0) + I(re74==0), 
#' 			data = LaLonde, ATT = TRUE)
#' balance(fit)
#'
#' @export
#' 
balance<-function(object, ...)
{
  UseMethod("balance")
}

#' Calculates the pre- and post-weighting difference in standardized means for covariate within each contrast
#' @param object A CBPS, npCBPS, or CBMSM object.
#' @param ... Additional arguments to be passed to balance.
#'
#' @export
#' 
balance.CBPS<-function(object, ...){
  treats<-as.factor(object$y)
  treat.names<-levels(treats)
  X<-object$x
  bal<-matrix(rep(0,(ncol(X)-1)*2*length(treat.names)),ncol(X)-1,2*length(treat.names))
  baseline<-matrix(rep(0,(ncol(X)-1)*2*length(treat.names)),ncol(X)-1,2*length(treat.names))
  w<-object$weights
  cnames<-array()
  
  jinit<-ifelse(class(object)[1] == "npCBPS", 1, 2)
  for (i in 1:length(treat.names))
  {
    for (j in jinit:ncol(X))
    {
      bal[j-1,i]<-sum((treats==treat.names[i])*X[,j]*w)/sum(w*(treats==treat.names[i]))
      bal[j-1,i+length(treat.names)]<-bal[j-1,i]/sd(X[,j])
      baseline[j-1,i]<-mean(X[which(treats==treat.names[i]),j])
      baseline[j-1,i+length(treat.names)]<-baseline[j-1,i]/sd(X[,j])
    }
    cnames[i]<-paste0(treat.names[i],".mean")
    cnames[length(treat.names)+i]<-paste0(treat.names[i],".std.mean")
  }
  colnames(bal)<-cnames
  rownames(bal)<-colnames(X)[-1]
  colnames(baseline)<-cnames
  rownames(baseline)<-colnames(X)[-1]
  out<-list(balanced=bal,original=baseline)
  out
}

#' Calculates the pre- and post-weighting correlations between each covariate and the T
#' @param object A CBPS, npCBPS, or CBMSM object.
#' @param ... Additional arguments to be passed to balance.
#'
#' @export
#' 
balance.CBPSContinuous<-function(object, ...){
  treat<-object$y
  X<-object$x
  w<-object$weights
  n<-length(w)
  cnames<-array()
  
  if ("npCBPS" %in% class(object)){
    jinit<-1
    bal<-matrix(rep(0,ncol(X)),ncol(X),1)
    baseline<-matrix(rep(0,ncol(X)),ncol(X),1)
    for (j in 1:ncol(X))
    {
      bal[j,1]<-(mean(w*X[,j]*treat) - mean(w*X[,j])*mean(w*treat)*n/sum(w))/(sqrt(mean(w*X[,j]^2) - mean(w*X[,j])^2*n/sum(w))*sqrt(mean(w*treat^2) - mean(w*treat)^2*n/sum(w)))
      baseline[j,1]<-cor(treat, X[,j], method = "pearson")
    }
    rownames(bal)<-colnames(X)
    rownames(baseline)<-colnames(X)   
  }
  else{
    bal<-matrix(rep(0,(ncol(X)-1)),ncol(X)-1,1)
    baseline<-matrix(rep(0,(ncol(X)-1)),ncol(X)-1,1)
    for (j in 2:ncol(X))
    {
      bal[j-1,1]<-(mean(w*X[,j]*treat) - mean(w*X[,j])*mean(w*treat)*n/sum(w))/(sqrt(mean(w*X[,j]^2) - mean(w*X[,j])^2*n/sum(w))*sqrt(mean(w*treat^2) - mean(w*treat)^2*n/sum(w)))
      baseline[j-1,1]<-cor(treat, X[,j], method = "pearson")
    }
    rownames(bal)<-colnames(X)[-1]
    rownames(baseline)<-colnames(X)[-1]      
  }
  
  colnames(bal)<-"Pearson Correlation"
  colnames(baseline)<-"Pearson Correlation"
  out<-list(balanced=bal,unweighted=baseline)
  out
}	
#########################################
