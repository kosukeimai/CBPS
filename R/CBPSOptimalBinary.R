CBPSOptimal.2Treat<-function(treat, X, baselineX, diffX, iterations, ATT, standardize = standardize)
{
  #will add ATT=1 later.
  probs.min<- 1e-6
  n<-dim(X)[1]
  X1=baselineX
  X1new=cbind(X[,1],X1)
  treat.orig<-treat
  treat<-sapply(treat,function(x) ifelse(x==levels(factor(treat))[2],1,0))
  
  baselineX=as.matrix(baselineX)
  diffX=as.matrix(diffX)
  if(dim(baselineX)[2]+dim(diffX)[2]+1 > dim(X)[2])
  {
    bal.only=3
    xcov=NULL
  }else if(dim(baselineX)[2]+dim(diffX)[2]+1 == dim(X)[2]){
    xcov = rep(1,dim(baselineX)[2]+dim(diffX)[2]+1)
    xcov = diag(xcov)
    bal.only=1
  }else{
    stop("Invalid baseline and diff models.")
  }	
  ATT.wt.func1<-function(beta.curr,X.wt=X){
    X<-as.matrix(X.wt)
    n<-dim(X)[1]
    n.c<-sum(treat==0)
    n.t<-sum(treat==1)
    theta.curr<-as.vector(X%*%beta.curr)
    probs.curr<-(1+exp(-theta.curr))^-1
    probs.curr<-pmin(1-probs.min,probs.curr)
    probs.curr<-pmax(probs.min,probs.curr)	
    w1<-(n/n.t*(treat-probs.curr)/(1-probs.curr))
    w1[treat==1] <- n/n.t
    w1
  }
  
  gmm.func1<-function(beta.curr,invV=NULL,option=NULL){
    probs.min<- 1e-6
    n<-dim(X)[1]
    
    if(ATT == 1)
    {
      n.t<-sum(treat==1)
      w.curr.del1.att = 1/n.t*t(X1new)%*%ATT.wt.func1(beta.curr,X)
      w.curr.del3.att = 1/n.t*t(addX)%*%(n/n.t*treat-1)
      gbar = c(w.curr.del1.att,w.curr.del3.att)
      
    }else{
      theta.curr<-as.vector(X%*%beta.curr)
      probs.curr<-(1+exp(-theta.curr))^-1
      probs.curr<-pmin(1-probs.min,probs.curr)
      probs.curr<-pmax(probs.min,probs.curr)	
      probs.curr<-as.vector(probs.curr)		
      w.curr<-treat/probs.curr-(1-treat)/(1-probs.curr)    #(probs.curr-1+treat)^-1 #need to check				
      w.curr<-as.vector(w.curr)
      
      w.curr.del<-1/n * t(X)%*%(w.curr)
      w.curr.del<-as.vector(w.curr.del)
      
      ##Generate the vector of mean imbalance by weights.
      w.curr.del1<-1/n * t(X1new)%*%(w.curr)
      w.curr.del1<-as.vector(w.curr.del1)
      
      ### Generate the vector of mean imbalance by weights for h2()
      addX=diffX #as.vector(rep(1,length(X[,1])))
      w.curr3 = treat/probs.curr - 1 #(1-probs.curr)/(probs.curr-1+treat)  #need to check
      w.curr.del3 = 1/n*t(addX)%*%(w.curr3)
      w.curr.del3<-as.vector(w.curr.del3)
      w.curr3 = as.vector(w.curr3)
      ##Generate g-bar, as in the paper.
      
      if(is.null(option))
      {
        gbar<-c(w.curr.del1,w.curr.del3)
      }else if(option == "CBPS")
      {
        gbar <- w.curr.del	
      }
      
    }			 				 								
    ##Generate the covariance matrix used in the GMM estimate.
    if(is.null(invV))
    {
      if(ATT==1)
      {
        #need to fill in this part
      }else{		
        X.1<-X1new*((1-probs.curr)*probs.curr)^-.5
        X.2<-addX*(1/probs.curr-1)^.5		
        X.1.1<- X1new*probs.curr^-.5
        X.1.2 <-addX*probs.curr^-.5
        V<-rbind(1/n*cbind(t(X.1)%*%X.1,t(X.1.1)%*%X.1.2),
                 1/n*cbind(t(X.1.2)%*%X.1.1,t(X.2)%*%X.2))
      }
      invV.g<-ginv(V)		
    }else{
      invV.g <-invV
    }
    
    ##Calculate the GMM loss.
    loss1<-as.vector(t(gbar)%*%invV.g%*%(gbar))		
    out1<-list("loss"=loss1, "invV"=invV.g)
    out1
  }
  gmm.loss1<-function(x,...) gmm.func1(x,...)$loss
  
  #initial value of beta by fitting glm
  glm1<-suppressWarnings(glm(treat~X-1,family=binomial))
  glm1$coef[is.na(glm1$coef)]<-0
  probs.glm<-glm1$fit
  glm1$fit<-probs.glm<-pmin(1-probs.min,probs.glm)
  glm1$fit<-probs.glm<-pmax(probs.min,probs.glm)	
  beta.curr<-glm1$coef
  beta.curr[is.na(beta.curr)]<-0
  glm.beta.curr<-glm(treat~X-1,family=binomial)$coefficients
  
  #initial value of beta by fitting CBPS	
  invV2<-ginv(t(X)%*%X)	
  cbps.beta.curr=optim(glm.beta.curr,gmm.loss1,invV=invV2,option="CBPS")$par
  
  gmm.init = glm.beta.curr
  if(bal.only ==1) #exact constrained
  {
    opt.bal<-optim(gmm.init, gmm.loss1,method="BFGS",invV=xcov)
    #beta.bal <- opt.bal$par	
    opt1<-opt.bal
  }	  
  if(bal.only==3) #solve gmm
  {
    gmm.glm.init<-optim(glm.beta.curr, gmm.loss1, control=list("maxit"=iterations), method="BFGS", hessian=TRUE)
    gmm.cbps.init<-optim(cbps.beta.curr, gmm.loss1, control=list("maxit"=iterations), method="BFGS", hessian=TRUE)
    if(gmm.glm.init$val<gmm.cbps.init$val) opt1<-gmm.glm.init else opt1<-gmm.cbps.init
  }
  
  ##Generate probabilities
  beta.opt<-opt1$par
  theta.opt<-as.vector(X%*%beta.opt)
  probs.opt<-(1+exp(-theta.opt))^-1
  probs.opt<-pmin(1-probs.min,probs.opt)
  probs.opt<-pmax(probs.min,probs.opt)
  
  beta.opt<-opt1$par
  ##Generate weights
  
  w.opt<-abs((probs.opt-1+treat)^-1)
  norm1<-norm2<-1
  if (standardize)
  {	
    norm1<-sum(treat/probs.opt)
    norm2<-sum((1-treat)/(1-probs.opt))
  }		
  w.opt<-(treat == 1)/probs.opt/norm1 + (treat == 0)/(1-probs.opt)/norm2	
  
  J.opt<-gmm.func1(beta.opt,invV=NULL)$loss
  residuals<-treat-probs.opt
  deviance <- -2*c(sum(treat*log(probs.opt)+(1-treat)*log(1-probs.opt)))
  nulldeviance <- -2*c(sum(treat*log(mean(treat))+(1-treat)*log(1-mean(treat))))
  
  #get vcov
  #G
  if(ATT==1)
  {
    #fill in this part later
  }else{
    XG.1 <- -sqrt(abs(treat-probs.opt)/(probs.opt*(1-probs.opt)))*X	
    XG.12 <- -sqrt(abs(treat-probs.opt)/(probs.opt*(1-probs.opt)))*X1new
    XW.1 <- X1new*(probs.opt-1+treat)^-1
    
    XG.2 <- -sqrt(treat*(1-probs.opt)/probs.opt)*X
    XG.22 <- -sqrt(treat*(1-probs.opt)/probs.opt)*diffX	
    XW.2 <- diffX*(treat/probs.opt-1)
  }
  W1<-rbind(t(XW.1),t(XW.2))
  G<-cbind(t(XG.1)%*%XG.12,t(XG.2)%*%XG.22)/n
  
  Omega <- (W1%*%t(W1)/n)
  #g	
  W<-gmm.func1(beta.opt,invV=NULL)$invV
  vcov<-ginv(G%*%W%*%t(G))%*%G%*%W%*%Omega%*%W%*%t(G)%*%ginv(G%*%W%*%t(G))
  
  output<-list("coefficients"=beta.opt,"fitted.values"=probs.opt,"deviance"=deviance,"weights"=w.opt,
               "y"=treat,"x"=X,"converged"=opt1$conv,"J"=J.opt,"var"=vcov, 
               "mle.J"=gmm.loss1(glm1$coef))
  
  class(output)<- c("CBPS","glm","lm")    
  output	
}