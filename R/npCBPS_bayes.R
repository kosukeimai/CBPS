###----------------------------------------------------------------
### Test of an idea for allowing <w,z> not equal to exactly zero, 
### by actually getting <w,z>=eta, for small eta....

npCBPS = function(X,T,iterations,ATT,standardize, amat, corprior=.1, print.level=0){  
  
  #Functions which will be called later, having to do with Taylor approximations of
  #log(.) near zero.

  llog = function(z, eps){
    ans <- z
    avoidNA <- !is.na(z)
    lo <- (z < eps) & avoidNA
    ans[lo] <- log(eps) - 1.5 + 2 * z[lo]/eps - 0.5 * (z[lo]/eps)^2
    ans[!lo] <- log(z[!lo])
    ans
  }
  
  llogp = function(z, eps){
    ans <- z
    avoidNA <- !is.na(z)
    lo <- (z < eps) & avoidNA
    ans[lo] <- 2/eps - z[lo]/eps^2
    ans[!lo] <- 1/z[!lo]
    ans
  }
  
  llogpp = function(z, eps){
    ans <- z
    avoidNA <- !is.na(z)
    lo <- (z < eps) & avoidNA
    ans[lo] <- -1/eps^2
    ans[!lo] <- -1/z[!lo]^2
    ans
  }  
  
  logelr = function(x, mu, lam) 
  {
    x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)
    if (n <= p) 
      stop("Need more observations than variables in logelr.")
    mu <- as.vector(mu)
    if (length(mu) != p) 
      stop("Length of mean doesn't match number of variables in logelr.")
    z <- t(t(x) - mu)
    arg <- 1 + z %*% lam
    return(-sum(llog(arg, eps=eps)))
  }
  
  
  # Center and rescale X prior to any other manipulations
  X=scale(as.matrix(X),center=TRUE,scale=TRUE)
  orig.X<-X
  n = nrow(X)
  eps=1/n
  
  #Operate on orthogonalization of X
  #X=svd(X)$u
  X=X%*%solve(chol(var(X)))
  
  #If user supplies the contraint matrix (unusual)  
  if(!is.null(amat)) {
    z=amat
    ncon=ncol(z)  #number of constraints 
    rm(amat) 
  }  
  
  #For continuous treatment & no user-supplied constraint matrix
  if(is.null(amat) && is.numeric(T)) {
    if (length(unique(T)) <= 4)  warning("Interpreting as a continuous treatment because T is numeric. To solve for a binary or multi-valued treatment, make T a factor.")
    T=scale(T,center=TRUE,scale=TRUE)
    z=as.matrix(X*as.vector(T))
    #optional second-order moment
    #z=cbind(z, as.matrix((X^2)*as.vector(T)))
    offset=0
    
    #Add aditional constraints that E[wT*]= and E[wX*]=0, which is typically
    #necessary
    ncon_cor=ncol(z)  #get dimensionality before mean constraints
    #z=cbind(z,T,X-offset)
    z=cbind(z,T,X)
    ncon=ncol(z)
  }
  
#For factor treatments
#(and no user-supplied constraints)
if(is.null(amat) && is.factor(T)) {
  conds=length(levels(T))
  conds=conds-1  #If we're going to drop last condition
  dimX=dim(X)[2]
  n=dim(X)[1]
  ncon_cor=dimX*conds
   
  Td = as.matrix(model.matrix(~T-1,model.frame(~T-1),contrasts=FALSE)[1:n,1:conds])
  colnames(Td)=levels(T)[1:conds]
  Tdc = scale(Td,center=TRUE,scale=TRUE)
       
  z=t(sapply(seq(1:n),function(x) t(kronecker(Tdc[x,],X[x,]))))  
  
  #add X and Tdc columns to Z
  z=cbind(z,X,Tdc)
  ncon=ncol(z)
  rm(Td)
  }
     
  ###-------------------------------------------------------------------------------
  ### Priors for the distribution of eta, parameterized as correlation coefficients
  ### For now I do these only for X*T constraints, not the mean constraints, which will be assumed
  ### to hold exactly.  Though we could relax those too.
  
  eta_prior_sd=rep(corprior,ncon_cor)
  
  ###----------------------------------------------------------------------------
  
  ###---------------------------------------------------------
  ### Solve by BFGS (optim)
  ###---------------------------------------------------------
  TINY=sqrt(.Machine$double.xmin)
    
log_elgiven_eta=function(par,eta,z,eps,ncon_cor){
  ncon=ncol(z)
  gamma=par
  eta_long=as.matrix(c(eta, rep(0,ncon-ncon_cor)))
  #matrix version of eta for vectorization purposes 
  eta_mat=eta_long%*%c(rep(1,nrow(z)))
  
  arg = (n + t(gamma)%*%(eta_mat-t(z)))  
  #used to be:  arg = (1 + t(gamma)%*%(t(z)-eta_mat))  
  
  log_el=-sum(llog(z=arg,eps=1/n))
  return(log_el)
}

get.w=function(eta,z,iterations=iterations.getw){
  tol=.005   #0.01 seems too generous; .005 seems to work well
  gam.init=rep(0, ncon)  
  opt.gamma.given.eta=optim(par=gam.init, eta=eta, method="BFGS", fn=log_elgiven_eta, z=z, eps=eps, ncon_cor=ncon_cor, control=list(fnscale=1, maxit=iterations))  
  gam.opt.given.eta=opt.gamma.given.eta$par
  eta_long=as.matrix(c(eta, rep(0,ncon-ncon_cor)))
  #matrix version of eta for vectorization purposes 
  eta_mat=eta_long%*%c(rep(1,nrow(z)))
  arg_temp = (n + t(gam.opt.given.eta)%*%(eta_mat-t(z)))
  w=as.numeric(1/arg_temp)
  sum.w=sum(w)
  
  #scale: should sum to 1 when actually applied:  
  w_scaled=w/sum.w
  
  if (abs(1-sum.w)<=tol){log_el=-sum(log(w_scaled))}
  if (abs(1-sum.w)>=tol){log_el=-sum(log(w_scaled))-10^6*(1+abs(1-sum.w))}
  
  #try a differentiable objective
  #c1=10^4
  #penalty=(c1-c1*sum.w)^2
  #log_el = -sum(log(w_scaled))-penalty  
  R=list()
  R$w=w
  R$sumw=sum.w
  R$log_el=log_el
  R$grad.gamma=w*(eta_mat-t(z)) #gradient w.r.t. gamma
  return(R)
}

log_post = function(par,eta_prior_sd,z){ 
  #get log(p(eta))
  eta_now=par
  log_p_eta=sum(-.5*log(2*pi*eta_prior_sd^2) - (eta_now^2)/(2*eta_prior_sd^2))
  
  #get best log_el for this eta  
  el.out=get.w(eta_now,z)
  
  #put it together into log(post)
  c=n  #temp measure to rescale the log(p(eta))
  log_post=el.out$log_el+c*log_p_eta
  print(c(log_post, el.out$log_el, log_p_eta))  
  return(log_post)
}

  
###-----------------------------------------------------------
### The main event
###-----------------------------------------------------------

#Now the outer optimization over eta
#initialize eta with actual correlations so we can get a good solution to start
#if (is.numeric(T)){eta.init=sapply(seq(1:ncon_cor), function(x) cor(X[,x],T))}
if (is.numeric(T)){eta.init=sapply(seq(1:ncon_cor), function(x) cor(z[,x],T))}

if (is.factor(T)){eta.init=colMeans(z[,1:ncon_cor])}

eta_prior_sd=rep(corprior,ncon_cor)

#with box constraints:
iterations.getw=100  #iterations of optim for the get.w (conditional on eta) search
iterations.eta=20 #iterations for the search over eta.

eta.optim.out=optim(par=eta.init, method="L-BFGS-B", fn=log_post, eta_prior_sd=eta_prior_sd,z=z,
                    upper=1*abs(eta.init)+rep(.01,ncon_cor),lower=-1*abs(eta.init)-rep(0.01,ncon_cor), control=list(fnscale=-1,factr=1e15, maxit=iterations.eta))
eta.opt=eta.optim.out$par
el.out.opt=get.w(eta=eta.opt,z)
log.el.opt=el.out.opt$log_el

w0=el.out.opt$w
w=w0/mean(w0)

meanflag<-0
  if (abs(sum(w0)-1) > 10^-2){
    warning("Weights do not have mean 1. Convergence has failed or the solution does not exist.")
    meanflag<-1
  }
  
  # Find the appropriate ATT if binary
  if (ATT & (length(unique(T)) == 2))
  {
    print(paste0("Finding ATT with T=", levels(T)[ATT], " as the treatment.  Set ATT=", 
                 as.character(3-ATT)," to find ATT with T=", levels(T)[3-ATT], " as the treatment"))
    treat<-levels(T)[ATT]
    n.t<-sum(T == treat)
    pscore<-1/w*(T == treat) + (1-1/w)*(T != treat)
    w<-abs(n/n.t*((T == treat)-pscore)/(1-pscore))
  }
  
  #standardize the means no matter what
  w<-w/mean(w)
  
  #if (standardize){
   # if (meanflag) warning(paste0("Standardizing weights to mean 1.  Unstandardized mean of weights is ",mean(w)))
  #  w<-w/mean(w)
  #}
  
  #rename gamma for reporting out; don't rescale it.
  #gamma = gam_opt
  
  
  R=list("fitted.values"=1/(w/n),"deviance"=2*sum(log(1/w)),"y"=T,"x"=orig.X,"converged"=eta.optim.out$conv, 
         "constraints" = z, "weights"=w, "eta"=eta.opt, "log_el"=log.el.opt)
  
  class(R)<-"npCBPS"
  
  return(R)
}
