CBPS.2Treat<-function(treat, X, method, k, XprimeX.inv, bal.only, iterations, ATT, standardize, twostep, sample.weights, ...){
  
  # There is probably an X'X missing somewhere that is causing these variance problems.
  probs.min<- 1e-6
  
  treat.orig<-treat
  treat<-sapply(treat,function(x) ifelse(x==levels(factor(treat))[2],1,0))
  if(ATT == 2)  treat<-1-treat
  
  if (ATT == 1){
    print(paste0("Finding ATT with T=",as.character(levels(factor(treat.orig))[2]),
                 " as the treatment.  Set ATT=2 to find ATT with T=",
                 as.character(levels(factor(treat.orig))[1])," as the treatment"))    
  }
  if (ATT == 2){
    print(paste0("Finding ATT with T=",as.character(levels(factor(treat.orig))[1]),
                 " as the treatment.  Set ATT=1 to find ATT with T=",
                 as.character(levels(factor(treat.orig))[2])," as the treatment"))    
  }
  
  ##Note: Sample weights sum to n and n.c and n.t measured wrt sample weigths
  sample.weights<-sample.weights/mean(sample.weights)
  n<-dim(X)[1]
  n.c<-sum(sample.weights[treat==0])
  n.t<-sum(sample.weights[treat==1])
  
  ##Generates ATT weights.	Called by loss function, etc.
  ATT.wt.func<-function(beta.curr,X.wt=X){
    X<-as.matrix(X.wt)
    n.c<-sum(sample.weights[treat==0])
    n.t<-sum(sample.weights[treat==1])
    n<-n.c+n.t
    theta.curr<-as.vector(X%*%beta.curr)
    probs.curr<-(1+exp(-theta.curr))^-1
    probs.curr<-pmin(1-probs.min,probs.curr)
    probs.curr<-pmax(probs.min,probs.curr)	
    w1<-n/n.t*(treat-probs.curr)/(1-probs.curr)
    w1
  }
  
  ##The gmm objective function--given a guess of beta, constructs the GMM J statistic.  Used for vanilla binary treatments
  gmm.func<-function(beta.curr,X.gmm=X,ATT.gmm=ATT,invV=NULL,sample.weights0=sample.weights){
    sample.weights<-sample.weights0
    ##Designate a few objects in the function.
    X<-as.matrix(X.gmm)
    ATT<-ATT.gmm
    
    ##Designate sample size, number of treated and control observations,
    ##theta.curr, which are used to generate probabilities.
    ##Trim probabilities, and generate weights.
    n.c<-sum(sample.weights[treat==0])
    n.t<-sum(sample.weights[treat==1])
    n<-n.c+n.t
    theta.curr<-as.vector(X%*%beta.curr)
    probs.curr<-(1+exp(-theta.curr))^-1
    probs.curr<-pmin(1-probs.min,probs.curr)
    probs.curr<-pmax(probs.min,probs.curr)	
    probs.curr<-as.vector(probs.curr)
    if(ATT){
      w.curr<-ATT.wt.func(beta.curr)}
    else{
      w.curr<-(probs.curr-1+treat)^-1}
    
    ##Generate the vector of mean imbalance by weights.
    w.curr.del<-1/(n)*t(sample.weights*X)%*%(w.curr)
    w.curr.del<-as.vector(w.curr.del)
    w.curr<-as.vector(w.curr)
    
    ##Generate g-bar, as in the paper.
    gbar<-c( 1/n*t(sample.weights*X)%*%(treat-probs.curr),w.curr.del)
    
    ##Generate the covariance matrix used in the GMM estimate.
    ##Was for the initial version that calculates the analytic variances.
    if(is.null(invV))
    {
      if(ATT){
        X.1<-sample.weights^.5*X*((1-probs.curr)*probs.curr)^.5
        X.2<-sample.weights^.5*X*(probs.curr/(1-probs.curr))^.5
        X.1.1<-sample.weights^.5*X*(probs.curr)^.5
      }
      else{
        X.1<-sample.weights^.5*X*((1-probs.curr)*probs.curr)^.5
        X.2<-sample.weights^.5*X*(probs.curr*(1-probs.curr))^-.5		
        X.1.1<- sample.weights^.5*X
      }
      if (ATT){
        V<-rbind(1/n*cbind(t(X.1)%*%X.1,t(X.1.1)%*%X.1.1)*n/sum(treat),
                 1/n*cbind(t(X.1.1)%*%X.1.1*n/sum(treat),t(X.2)%*%X.2*n^2/sum(treat)^2))
      }
      else{
        V<-rbind(1/n*cbind(t(X.1)%*%X.1,t(X.1.1)%*%X.1.1),
                 1/n*cbind(t(X.1.1)%*%X.1.1,t(X.2)%*%X.2))
      }
      invV<-ginv(V)
    }			   
    
    ##Calculate the GMM loss.
    loss1<-as.vector(t(gbar)%*%invV%*%(gbar))		
    out1<-list("loss"=loss1, "invV"=invV)
    out1
  }
  
  gmm.loss<-function(x,...) gmm.func(x,...)$loss
  
  ##Loss function for balance constraints, returns the squared imbalance along each dimension.
  bal.loss<-function(beta.curr,sample.weights0=sample.weights){
    sample.weights<-sample.weights0
    ##Generate theta and probabilities.
    theta.curr<-as.vector(X%*%beta.curr)
    probs.curr<-(1+exp(-theta.curr))^-1
    probs.curr<-pmin(1-probs.min,probs.curr)
    probs.curr<-pmax(probs.min,probs.curr)
    ##Generate weights.
    if(ATT){
      w.curr<-1/n*ATT.wt.func(beta.curr)
    }
    else{
      w.curr<-1/n*(probs.curr-1+treat)^-1
    }
    ##Generate mean imbalance.
    Xprimew <- t(sample.weights*X)%*%(w.curr)
    loss1<-abs(t(Xprimew)%*%XprimeX.inv%*%Xprimew)
    loss1
  }
  
  ##Does not work with ATT.  Need to fix this at some point.
  gmm.gradient<-function(beta.curr, invV, ATT.gmm=ATT,sample.weights0=sample.weights)
  {
    sample.weights<-sample.weights0
    n.c<-sum(sample.weights[treat==0])
    n.t<-sum(sample.weights[treat==1])
    n<-n.c+n.t
    ATT<-ATT.gmm
    theta.curr<-as.vector(X%*%beta.curr)
    probs.curr<-(1+exp(-theta.curr))^-1
    probs.curr<-pmin(1-probs.min,probs.curr)
    probs.curr<-pmax(probs.min,probs.curr)
    
    ##Generate the vector of mean imbalance by weights.
    if (ATT){
      w.curr<-ATT.wt.func(beta.curr)
    }
    else{
      w.curr<-(probs.curr-1+treat)^-1
    }
    w.curr.del<-1/n*t(X*sample.weights)%*%(w.curr)
    w.curr.del<-as.vector(w.curr.del)
    w.curr<-as.vector(w.curr)
    
    ##Generate g-bar, as in the paper.
    gbar<-c(1/n*t(X*sample.weights)%*%(treat-probs.curr),w.curr.del)
    
    ##Calculate derivative of g-bar
    if (ATT){
      # Need to update here
      dw<- -n/n.t*probs.curr/(1 - probs.curr)
      dw[treat==1]<-0
      dgbar<-cbind(1/n*t(-X*sample.weights*probs.curr*(1-probs.curr))%*%X, 
                   1/n.t*t(X*dw*sample.weights)%*%X)
    }
    else{
      dgbar<-cbind(-1/n*t(X*sample.weights*probs.curr*(1-probs.curr))%*%X,
                   -1/n*t(X*sample.weights*(treat - probs.curr)^2/(probs.curr*(1-probs.curr)))%*%X)
    }
    out<-2*dgbar%*%invV%*%gbar
  }
  
  bal.gradient<-function(beta.curr,sample.weights0=sample.weights)
  {
    ##Generate theta and probabilities.
    sample.weights<-sample.weights0
    theta.curr<-as.vector(X%*%beta.curr)
    probs.curr<-(1+exp(-theta.curr))^-1
    probs.curr<-pmin(1-probs.min,probs.curr)
    probs.curr<-pmax(probs.min,probs.curr)
    ##Generate weights.
    if(ATT) w.curr<-1/n*ATT.wt.func(beta.curr)
    else w.curr<-1/n*(probs.curr-1+treat)^-1
    
    if (ATT){
      dw2<- -n/n.t*probs.curr/(1 - probs.curr)
      dw2[treat==1]<-0
      dw<-1/n*t(X*dw2)
    }
    else{
      dw<-1/n*t(-X*(treat-probs.curr)^2/(probs.curr*(1-probs.curr)))
    }
    ##Generate mean imbalance.
    Xprimew <- t(X)%*%(w.curr*sample.weights)
    loss1<-t(Xprimew)%*%XprimeX.inv%*%Xprimew
    out<-sapply(2*dw%*%X%*%XprimeX.inv%*%Xprimew, function (x) ifelse((x > 0 & loss1 > 0) | (x < 0 & loss1 < 0), abs(x), -abs(x))) 
    out
  }
  
  n<-length(treat)
  n.t<-sum(treat==1)
  
  ##GLM estimation
  glm1<-suppressWarnings(glm(treat~X-1,family=binomial))
  glm1$coef[is.na(glm1$coef)]<-0
  probs.glm<-glm1$fit
  glm1$fit<-probs.glm<-pmin(1-probs.min,probs.glm)
  glm1$fit<-probs.glm<-pmax(probs.min,probs.glm)	
  beta.curr<-glm1$coef
  beta.curr[is.na(beta.curr)]<-0
  
  alpha.func<-function(alpha) gmm.loss(beta.curr*alpha)
  beta.curr<-beta.curr*optimize(alpha.func,interval=c(.8,1.1))$min
  
  ##Generate estimates for balance and CBPSE
  gmm.init<-beta.curr
  this.invV<-gmm.func(gmm.init)$invV
  
  if (twostep)
  {
    opt.bal<-optim(gmm.init, bal.loss, control=list("maxit"=iterations), method="BFGS", gr = bal.gradient, hessian=TRUE)
  }
  else
  {
    opt.bal<-optim(gmm.init, bal.loss, control=list("maxit"=iterations), method="BFGS", hessian=TRUE)
  }
  beta.bal<-opt.bal$par
  
  if(bal.only) opt1<-opt.bal
  
  if(!bal.only)
  {
    if (twostep)
    {
      gmm.glm.init<-optim(gmm.init, gmm.loss, control=list("maxit"=iterations), method="BFGS", hessian=TRUE, gr = gmm.gradient, invV = this.invV)
      gmm.bal.init<-optim(beta.bal, gmm.loss, control=list("maxit"=iterations), method="BFGS", hessian=TRUE, gr = gmm.gradient, invV = this.invV)
    }
    else
    {
      gmm.glm.init<-optim(gmm.init, gmm.loss, control=list("maxit"=iterations), method="BFGS", hessian=TRUE)
      gmm.bal.init<-optim(beta.bal, gmm.loss, control=list("maxit"=iterations), method="BFGS", hessian=TRUE)
    }
    if(gmm.glm.init$val<gmm.bal.init$val) opt1<-gmm.glm.init else opt1<-gmm.bal.init
  }
  
  ##Generate probabilities
  beta.opt<-opt1$par
  theta.opt<-as.vector(X%*%beta.opt)
  probs.opt<-(1+exp(-theta.opt))^-1
  probs.opt<-pmin(1-probs.min,probs.opt)
  probs.opt<-pmax(probs.min,probs.opt)
  
  ##Generate weights
  if(ATT){
    w.opt<-abs(ATT.wt.func(beta.opt)) 
  }else{
    w.opt<-abs((probs.opt-1+treat)^-1)
  }
  
  if (standardize)
  {
    if (ATT)
    {
      norm1<-sum(treat*sample.weights*n/n.t)
      norm2<-sum((1-treat)*sample.weights*n/n.t*(treat-probs.opt)/(1-probs.opt))
    }
    else
    {
      norm1<-sum(treat*sample.weights/probs.opt)
      norm2<-sum((1-treat)*sample.weights/(1-probs.opt))
    }
  }
  else {
    norm1 <- 1
    norm2 <- 1
  }
                
  if (ATT)
  {
    w.opt<-(treat == 1)*n/sum(treat == 1)/norm1 + abs((treat == 0)*n/sum(treat == 1)*((treat - probs.opt)/(1-probs.opt))/norm2)
  }
  else
  {		
    w.opt<-(treat == 1)/probs.opt/norm1 + (treat == 0)/(1-probs.opt)/norm2
  }
  w.opt<-w.opt*sample.weights
  
  J.opt<-ifelse(twostep, gmm.func(beta.opt, invV = this.invV)$loss, gmm.loss(beta.opt))
  
  residuals<-treat-probs.opt
  deviance <- -2*c(sum(treat*sample.weights*log(probs.opt)+(1-treat)*sample.weights*log(1-probs.opt)))
  nulldeviance <- -2*c(sum(treat*sample.weights*log(mean(treat))+(1-treat)*sample.weights*log(1-mean(treat))))
  
  XG.1<- -X*probs.opt*(1-probs.opt)*sample.weights
  XW.1<- X*(treat-probs.opt)*sample.weights^.5
  if(ATT){
    XW.2<-X*ATT.wt.func(beta.opt)*sample.weights
    dw2<- -n/n.t*probs.opt/(1 - probs.opt)
    dw2[treat==1]<-0
    XG.2 <- X*dw2*sample.weights
  } 
  else{
    XW.2 <- X*(probs.opt-1+treat)^-1*sample.weights^.5
    XG.2 <- -X*(treat - probs.opt)^2/(probs.opt*(1-probs.opt))*sample.weights
  }
  
  if (bal.only){
    G<-cbind(t(XG.2)%*%X)/n
    W1<-rbind(t(XW.2))	  
    W<-XprimeX.inv
  }
  else{
    G<-cbind(t(XG.1)%*%X,t(XG.2)%*%X)/n
    W1<-rbind(t(XW.1),t(XW.2))
    if (twostep){
      W <- this.invV
    }
    else{
      W <- gmm.func(beta.opt)$invV
    }
  }
  
  Omega<-(W1%*%t(W1)/n)
  GWGinvGW <- ginv(G%*%W%*%t(G))%*%G%*%W 
  vcov<- GWGinvGW%*%Omega%*%t(GWGinvGW)
  
  output<-list("coefficients"=matrix(beta.opt, ncol=1),"fitted.values"=probs.opt, "linear.predictor" = X%*%beta.opt,
               "deviance"=deviance,"weights"=w.opt,
               "y"=treat,"x"=X,"converged"=opt1$conv,"J"=J.opt,"var"=vcov, 
               "mle.J"=ifelse(twostep, gmm.func(glm1$coef, invV = this.invV)$loss, gmm.loss(glm1$coef)))
  
  class(output)<- c("CBPS","glm","lm")
  
  output
}
