#npCBPS parses the formula object and passes the result to npCBPS.fit


#' Non-Parametric Covariate Balancing Propensity Score (npCBPS) Estimation
#' 
#' \code{npCBPS} is a method to estimate weights interpretable as (stabilized)
#' inverse generlized propensity score weights, w_i = f(T_i)/f(T_i|X), without
#' actually estimating a model for the treatment to arrive at f(T|X) estimates.
#' In brief, this works by maximizing the empirical likelihood of observing the
#' values of treatment and covariates that were observed, while constraining
#' the weights to be those that (a) ensure balance on the covariates, and (b)
#' maintain the original means of the treatment and covariates.
#' 
#' In the continuous treatment context, this balance on covariates means zero
#' correlation of each covariate with the treatment. In binary or categorical
#' treatment contexts, balance on covariates implies equal means on the
#' covariates for observations at each level of the treatment.  When given a
#' numeric treatment the software handles it continuously. To handle the
#' treatment as binary or categorical is must be given as a factor.
#' 
#' Furthermore, we apply a Bayesian variant that allows the correlation of each
#' covariate with the treatment to be slightly non-zero, as might be expected
#' in a a given finite sample.
#' 
#' Estimates non-parametric covariate balancing propensity score weights.
#' 
#' @aliases npCBPS npCBPS.fit
#' @param formula An object of class \code{formula} (or one that can be coerced
#' to that class): a symbolic description of the model to be fitted.
#' @param data An optional data frame, list or environment (or object coercible
#' by as.data.frame to a data frame) containing the variables in the model. If
#' not found in data, the variables are taken from \code{environment(formula)},
#' typically the environment from which \code{CBPS} is called.
#' @param na.action A function which indicates what should happen when the data
#' contain NAs. The default is set by the na.action setting of options, and is
#' na.fail if that is unset.
#' @param corprior Prior hyperparameter controlling the expected amount of
#' correlation between each covariate and the treatment. Specifically, the
#' amount of correlation between the k-dimensional covariates, X, and the
#' treatment T after weighting is assumed to have prior distribution
#' MVN(0,sigma^2 I_k). We conceptualize sigma^2 as a tuning parameter to be
#' used pragmatically. It's default of 0.1 ensures that the balance constraints
#' are not too harsh, and that a solution is likely to exist. Once the
#' algorithm works at such a high value of sigma^2, the user may wish to
#' attempt values closer to 0 to get finer balance.
#' @param print.level Controls verbosity of output to the screen while npCBPS
#' runs. At the default of print.level=0, little output is produced. It
#' print.level>0, it outputs diagnostics including the log posterior
#' (log_post), the log empirical likelihood associated with the weights
#' (log_el), and the log prior probability of the (weighted) correlation of
#' treatment with the covariates.
#' @param treat A vector of treatment assignments.  Binary or multi-valued
#' treatments should be factors.  Continuous treatments should be numeric.
#' @param X A covariate matrix.
#' @param ... Other parameters to be passed.
#' @return \item{weights}{The optimal weights} \item{y}{The treatment vector
#' used} \item{x}{The covariate matrix} \item{model}{The model frame}
#' \item{call}{The matched call} \item{formula}{The formula supplied}
#' \item{data}{The data argument} \item{log.p.eta}{The log density for the
#' (weighted) correlation of the covariates with the treatment, given the
#' choice of prior (\code{corprior})} \item{log.el}{The log empirical
#' likelihood of the observed data at the chosen set of IPW weights.}
#' \item{eta}{A vector describing the correlation between the treatment and
#' each covariate on the weighted data at the solution.} \item{sumw0}{The sum
#' of weights, provided as a check on convergence. This is always 1 when
#' convergence occurs unproblematically. If it differs from 1 substantially, no
#' solution perfectly satisfying the conditions was found, and the user may
#' consider a larger value of \code{corprior}.}
#' @author Christian Fong, Chad Hazlett, and Kosuke Imai
#' @references Fong, Christian, Chad Hazlett, and Kosuke Imai.  ``Parametric
#' and Nonparametric Covariate Balancing Propensity Score for General Treatment
#' Regimes.'' Unpublished Manuscript.
#' \url{http://imai.princeton.edu/research/files/CBGPS.pdf}
#' @examples
#' 
#' ##Generate data
#' data(LaLonde)
#' 
#' ## Restricted two only two covariates so that it will run quickly.
#' ## Performance will remain good if the full LaLonde specification is used
#' fit <- npCBPS(treat ~ age + educ, data = LaLonde, corprior=.1/nrow(LaLonde))
#' plot(fit)
#' 
npCBPS <- function(formula, data, na.action, corprior=.01, print.level=0, ...) {
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
  else matrix(, NROW(Y), 0L)
  
  X<-X[,apply(X,2,sd)>0]
  
  
  fit <- eval(call("npCBPS.fit", X = X, treat = Y, corprior = corprior, 
                   print.level = print.level))	
  
  fit$na.action <- attr(mf, "na.action")
  xlevels <- .getXlevels(mt, mf)
  fit$data<-data
  fit$call <- call
  fit$formula <- formula
  fit$terms<-mt
  fit
}
npCBPS.fit=function(treat, X, corprior, print.level, ...){
  D=treat
  rescale.orig=TRUE
  orig.X=X
  #pre-processesing: 
  X=X%*%solve(chol(var(X)))
  X=scale(X,center=TRUE, scale=TRUE)
  n=nrow(X)
  eps=1/n
  
  #Constraint matrix
  if (is.numeric(D)){
    print("Estimating npCBPS as a continuous treatment.  To estimate for a binary or multi-valued treatment, use a factor.")
    #re-orient each X to have positive correlation with T
    X=X%*%diag(as.vector(sign(cor(X,D))),nrow=ncol(X))
    D=scale(D,center=TRUE, scale=TRUE)  
    z=X*as.vector(D)
    z=cbind(z,X,D)
    ncon=ncol(z)
    ncon_cor=ncol(X)*ncol(D)
  }
  
  if(is.factor(D)){
    #For factor treatments
    Td=as.matrix(model.matrix(~D-1))
    conds=dim(Td)[2]  
    dimX=dim(X)[2]
    
    #Now divide each column of Td by it's sum
    colsums=apply(Td,2,sum)
    Td=Td%*%diag(1/colsums)
    
    #Now subtract the last column from each of the others, and remove the last
    subtractMat=Td[,conds]%*%t(as.matrix(rep(1, conds)))
    Td=Td-subtractMat
    Td=Td[,1:(conds-1)]
    
    #Center and rescale Td now
    Td=scale(x=Td, center = TRUE, scale=TRUE)
    
    #form matrix z that will be needed to setup contrasts
    z=matrix(NA,nrow=n,ncol=dimX*(conds-1))
    z=t(sapply(seq(1:n),function(x) t(kronecker(Td[x,],X[x,]))))
    
    
    #Check that correlation of Td with X is very close to colMeans of z
    cor.init=as.vector(t(apply(X = X,MARGIN = 2,function(x) cor(Td,x))))
    rescale.factors=cor.init/colMeans(z)
    if (print.level>0){print(rescale.factors)}
    
    #Add aditional constraints that E[wX*]=0, if desired
    #NB: I think we need another constraint to ensure something like E[wT*]=0
    ncon_cor=dim(z)[2]  #keep track of number of constraints not including the additional mean constraint
    z=cbind(z,X)
    ncon=dim(z)[2] #num constraints including mean constraints
    #rm(Td)
  }
  
  #-----------------------------------------------
  # Functions we will need
  #-----------------------------------------------
  
  llog = function(z, eps){
    ans = z
    avoidNA = !is.na(z)
    lo = (z < eps) & avoidNA
    ans[lo] = log(eps) - 1.5 + 2 * z[lo]/eps - 0.5 * (z[lo]/eps)^2
    ans[!lo] = log(z[!lo])
    ans
  }
  
  llogp = function(z, eps){
    ans = z
    avoidNA = !is.na(z)
    lo = (z < eps) & avoidNA
    ans[lo] = 2/eps - z[lo]/eps^2
    ans[!lo] = 1/z[!lo]
    ans
  }
  
  log_elgiven_eta=function(par,eta,z,eps,ncon_cor){
    ncon=ncol(z)
    gamma=par
    eta_long=as.matrix(c(eta, rep(0,ncon-ncon_cor)))
    #matrix version of eta for vectorization purposes 
    eta_mat=eta_long%*%c(rep(1,nrow(z)))
    
    arg = (n + t(gamma)%*%(eta_mat-t(z)))  
    #used to be:  arg = (1 + t(gamma)%*%(t(z)-eta_mat))  
    
    log_el=-sum(llog(z=arg,eps=eps))
    return(log_el)
  }
  
  get.w=function(eta,z, sumw.tol=0.05, eps){
    gam.init=rep(0, ncon)  
    opt.gamma.given.eta=optim(par=gam.init, eta=eta, method="BFGS", fn=log_elgiven_eta, z=z, eps=eps, ncon_cor=ncon_cor, control=list(fnscale=1))  
    gam.opt.given.eta=opt.gamma.given.eta$par
    eta_long=as.matrix(c(eta, rep(0,ncon-ncon_cor)))
    
    #matrix version of eta for vectorization purposes 
    eta_mat=eta_long%*%c(rep(1,nrow(z)))
    arg_temp = (n + t(gam.opt.given.eta)%*%(eta_mat-t(z)))
    
    #just use 1/x instead instead of the derivative of the pseudo-log
    w=as.numeric(1/arg_temp)
    sum.w=sum(w)
    
    #scale: should sum to 1 when actually applied:  
    w_scaled=w/sum.w    
    
    if (abs(1-sum.w)<=sumw.tol){log_el=-sum(log(w_scaled))}
    if (abs(1-sum.w)>=sumw.tol){log_el=-sum(log(w_scaled))-10^4*(1+abs(1-sum.w))}
    
    R=list()
    R$w=w
    R$sumw=sum.w
    R$log_el=log_el
    R$el.gamma=gam.opt.given.eta[1:ncon_cor]
    #R$grad.gamma=w*(eta_mat-t(z)) #gradient w.r.t. gamma
    return(R)
  }
  
  #------------------------
  # Some diagnostics:
  # (a) is the eta really the cor(X,T) you end up with?
  # (b) is balance on X and T (at 0) is maintained
  #------------------------
  
  ## eta= (.1,.1,...,.1) should produce w's that produce weighted cov = (.1,.1,...)
  #test.w=get.w(eta=rep(.1,ncon_cor),z)
  ##check convergence: is sumw near 1?
  #test.w$sumw 
  ##get w
  #wtest=test.w$w
  ##check weighted covariances: are they near 0.10?
  #sapply(seq(1,5), function(x) sum(X[,x]*T*wtest))
  ##means of X and T: are they near 0?
  #sapply(seq(1,5), function(x) sum(X[,x]*wtest))
  #sum(T*wtest)

  log_post = function(par,eta.to.be.scaled,eta_prior_sd,z, eps=eps, sumw.tol=.001){ 
    #get log(p(eta))
    eta_now=par*eta.to.be.scaled
    log_p_eta=sum(-.5*log(2*pi*eta_prior_sd^2) - (eta_now^2)/(2*eta_prior_sd^2))
    
    #get best log_el for this eta  
    el.out=get.w(eta=eta_now,z=z, sumw.tol=sumw.tol, eps=eps)
    #el.gamma=el.out$el.gamma
    #put it together into log(post)
    c=1 #in case we want to rescale the log(p(eta)), as sigma/c would.
    log_post=el.out$log_el+c*log_p_eta
    
    if(print.level>0){print(c(log_post, el.out$log_el, log_p_eta))}  
    return(log_post)
  }
  
  ###-----------------------------------------------------------
  ### The main event
  ###-----------------------------------------------------------
  
  #Now the outer optimization over eta
  #setup the prior
  eta_prior_sd=rep(corprior,ncon_cor)
  
  #get original correlations
  if (is.numeric(D)){eta.init=sapply(seq(1:ncon_cor), function(x) cor(X[,x],D))}
  
  #for factor treatment, there is probably a better analog to the intialize correlation,
  if (is.factor(D)){
    eta.init=cor.init
  }
  
  #get vector of 1's long enough to be our dummy that gets rescaled to form eta if we want
  #constant etas:
  eta.const=rep(1, ncon_cor)
  
  #note that as currently implemented, these are only the non-zero elements of eta that correspond
  # to cor(X,T). For additional constraints that hold down the mean of X and T we are assuming
  # eta=0 effectively.  They get padded in within the optimization.
  
  #Determine if we want to rescale 1's or rescale the original correlations
  #rescale.orig=FALSE
  if(rescale.orig==TRUE){eta.to.be.scaled=eta.init}else{eta.to.be.scaled=eta.const}
  
  eta.optim.out=optimize(f=log_post, interval=c(-1,1), eta.to.be.scaled=eta.to.be.scaled,
                         eps=eps, sumw.tol=.001, eta_prior_sd=eta_prior_sd,z=z, maximum=TRUE)
  
  #Some useful values:
  par.opt=eta.optim.out$maximum
  eta.opt=par.opt*eta.to.be.scaled
  
  log.p.eta.opt=sum(-.5*log(2*pi*eta_prior_sd^2) - (eta.opt^2)/(2*eta_prior_sd^2))
  
  el.out.opt=get.w(eta=eta.opt,z=z, eps=eps)
  sumw0=sum(el.out.opt$w)
  w=el.out.opt$w/sumw0
  log.el.opt=el.out.opt$log_el
  
  R=list()
  R$par=par.opt
  R$log.p.eta=log.p.eta.opt
  R$log.el=log.el.opt
  R$eta=eta.opt
  R$sumw0=sumw0  #sum of original w prior to any corrective rescaling
  R$weights=w
  R$y=D
  R$x=orig.X
    
  class(R)<-"npCBPS"
  
  return(R)
}

# Calls the appropriate plot function, based on the number of treatments
plot.npCBPS<-function(x, covars = NULL, silent = TRUE, ...){
  bal.x<-balance(x)
  
  if(is.numeric(x$y)) {out<-plot.CBPSContinuous(x, covars, silent, ...)}
  else  {out<-plot.CBPS(x, covars, silent, ...)}
  if(!is.null(out)) return(out)  
}

# Calls the appropriate balance function based on the number of treatments
balance.npCBPS<-function(object, ...){
  if(is.numeric(object$y)) {out<-balance.CBPSContinuous(object, ...)}
  else  {out<-balance.CBPS(object, ...)}
  out
}
