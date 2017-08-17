CBIV <- function(Tr, Z, X, iterations=1000, method="over", twostep = TRUE, twosided = TRUE, ...) {
  probs.min<-10^-6
  pZ <- mean(Z)
  k<-0
  
  score.only<-bal.only<-FALSE
  if(method=="mle") score.only<-TRUE
  if(method=="exact") bal.only<-TRUE
  
  X<-as.matrix(X)
  X<-cbind(1,X[,apply(X,2,sd)>0])
  names.X<-colnames(X)
  names.X[apply(X,2,sd)==0]<-"(Intercept)"
  
  #######Declare some constants and orthogonalize Xdf.
  X.orig<-X
  x.sd<-apply(as.matrix(X[,-1]),2,sd)
  Dx.inv<-diag(c(1,x.sd))
  diag(Dx.inv)<-1
  x.mean<-apply(as.matrix(X[,-1]),2,mean)
  X[,-1]<-apply(as.matrix(X[,-1]),2,FUN=function(x) (x-mean(x))/sd(x))
  if(k==0) k<-sum(diag(t(X)%*%X%*%ginv(t(X)%*%X)))
  k<-floor(k+.1)
  XprimeX.inv<-ginv(t(X)%*%X)
  
  n<-length(Tr)
  
  gmm.func <- function(beta.curr, invV = NULL, twosided)
  {
    if (twosided){
      beta.curr.c<-beta.curr[1:k]
      beta.curr.a<-beta.curr[k+(1:k)]
      
      baseline.prob<-(1 + exp(X%*%beta.curr.c) + exp(X%*%beta.curr.a))^-1
      
      probs.curr.c<-pmin(pmax(exp(X%*%beta.curr.c)*baseline.prob,probs.min),1-probs.min)
      probs.curr.a<-pmin(pmax(exp(X%*%beta.curr.a)*baseline.prob,probs.min),1-probs.min)
      probs.curr.n<-pmin(pmax(baseline.prob,probs.min),1-probs.min)
      
      sums<-probs.curr.c+probs.curr.a+probs.curr.n
      probs.curr.c<-probs.curr.c/sums
      probs.curr.a<-probs.curr.a/sums
      probs.curr.n<-probs.curr.n/sums
      
      w.curr<-cbind(Z*Tr/(pZ*(probs.curr.c + probs.curr.a)) - 1, 
                    (1-Z)*Tr/((1-pZ)*probs.curr.a) - 1,
                    Z*(1-Tr)/(pZ*probs.curr.n) - 1, 
                    (1-Z)*(1-Tr)/((1-pZ)*(probs.curr.c + probs.curr.n)) - 1)
      
      w.curr.del<-1/n*t(X)%*%w.curr
      w.curr.del<-as.matrix(w.curr.del)
      w.curr<-as.matrix(w.curr)
      
      gbar<-c(1/n*t(X)%*%((Z*Tr/(1-probs.curr.n) + (1-Z)*(1-Tr)/(1-probs.curr.a) - 1)*probs.curr.c),
              1/n*t(X)%*%((Z*Tr/(1-probs.curr.n) + (1-Z)*Tr/probs.curr.a - 1)*probs.curr.a),
              w.curr.del)
      
      if (is.null(invV))
      {
        X.1.1<-X*as.vector((pZ/(1 - probs.curr.n) + (1 - pZ)/(1 - probs.curr.a) - 1)*probs.curr.c^2)
        X.1.2<-X*as.vector((pZ/(1 - probs.curr.n) - 1)*probs.curr.a*probs.curr.c)
        X.1.3<-X*as.vector(probs.curr.c*((1 - probs.curr.n)^-1 - 1))
        X.1.4<-X*as.vector(-probs.curr.c)
        X.1.5<-X*as.vector(-probs.curr.c)
        X.1.6<-X*as.vector(probs.curr.c*((probs.curr.c + probs.curr.n)^-1 - 1))
        X.2.2<-X*as.vector((pZ/(1-probs.curr.n) + (1 - pZ)/probs.curr.a - 1)*probs.curr.a^2)
        X.2.3<-X*as.vector(probs.curr.a*((probs.curr.c+probs.curr.a)^-1 - 1))
        X.2.4<-X*as.vector(probs.curr.a*((probs.curr.a)^-1 - 1))
        X.2.5<-X*as.vector(-probs.curr.a)
        X.2.6<-X*as.vector(-probs.curr.a)
        X.3.3<-X*as.vector((pZ*(probs.curr.c + probs.curr.a))^-1 - 1)
        X.3.4<- -X
        X.3.5<- -X
        X.3.6<- -X
        X.4.4<-X*as.vector(((1-pZ)*probs.curr.a)^-1 - 1)
        X.4.5<- -X
        X.4.6<- -X
        X.5.5<-X*as.vector((pZ*probs.curr.n)^-1 - 1)
        X.5.6<- -X
        X.6.6<-X*as.vector(((1-pZ)*(probs.curr.c + probs.curr.n))^-1 - 1)
        
        V<-1/n*rbind(cbind(t(X.1.1)%*%X, t(X.1.2)%*%X, t(X.1.3)%*%X, t(X.1.4)%*%X, t(X.1.5)%*%X, t(X.1.6)%*%X),
                     cbind(t(X.1.2)%*%X, t(X.2.2)%*%X, t(X.2.3)%*%X, t(X.2.4)%*%X, t(X.2.5)%*%X, t(X.2.6)%*%X),
                     cbind(t(X.1.3)%*%X, t(X.2.3)%*%X, t(X.3.3)%*%X, t(X.3.4)%*%X, t(X.3.5)%*%X, t(X.3.6)%*%X),
                     cbind(t(X.1.4)%*%X, t(X.2.4)%*%X, t(X.3.4)%*%X, t(X.4.4)%*%X, t(X.4.5)%*%X, t(X.4.6)%*%X),
                     cbind(t(X.1.5)%*%X, t(X.2.5)%*%X, t(X.3.5)%*%X, t(X.4.5)%*%X, t(X.5.5)%*%X, t(X.5.6)%*%X),
                     cbind(t(X.1.6)%*%X, t(X.2.6)%*%X, t(X.3.6)%*%X, t(X.4.6)%*%X, t(X.5.6)%*%X, t(X.6.6)%*%X))
        invV<-ginv(V)
      }
    }
    else{
      probs.curr <- pmin(pmax(exp(X%*%beta.curr)/(1 + exp(X%*%beta.curr)), probs.min), 1-probs.min)
      
      w.curr <- cbind(Z*Tr/(pZ * probs.curr) - 1, Z*(1-Tr)/(pZ *(1 - probs.curr)) - 1)
      w.curr.del<-1/n*t(X)%*%w.curr
      w.curr.del<-as.matrix(w.curr.del)
      w.curr<-as.matrix(w.curr)
      
      gbar<-c(1/n*t(X)%*%(Tr*Z*(1-probs.curr) - Z*(1-Tr)*probs.curr),
              w.curr.del)
      
      if (is.null(invV))
      {
        X.1.1<-X*as.vector(pZ * probs.curr * (1 - probs.curr))
        X.1.2<-X*as.vector(1 - probs.curr)
        X.1.3<-X*as.vector(-probs.curr)
        X.2.2<-X*as.vector((pZ * probs.curr)^-1 - 1)
        X.2.3<- -X
        X.3.3<-X*as.vector((pZ * (1 - probs.curr))^-1 - 1)
        
        V<-1/n*rbind(cbind(t(X.1.1)%*%X, t(X.1.2)%*%X, t(X.1.3)%*%X),
                     cbind(t(X.1.2)%*%X, t(X.2.2)%*%X, t(X.2.3)%*%X),
                     cbind(t(X.1.3)%*%X, t(X.2.3)%*%X, t(X.3.3)%*%X))
        invV<-ginv(V)
      }
      
    }
    
    loss1<-as.vector(t(gbar)%*%invV%*%(gbar))
    out1<-list("loss"=loss1, "invV"=invV)
    out1
  }
  
  gmm.loss <- function(beta.curr, invV = NULL, twosided) gmm.func(beta.curr, invV, twosided = twosided)$loss
  
  gmm.gradient <- function(beta.curr, invV, twosided)
  {
    if(twosided){
      beta.curr.c<-beta.curr[1:k]
      beta.curr.a<-beta.curr[k+(1:k)]
      
      baseline.prob<-(1 + exp(X%*%beta.curr.c) + exp(X%*%beta.curr.a))^-1
      
      probs.curr.c<-pmin(pmax(exp(X%*%beta.curr.c)*baseline.prob,probs.min),1-probs.min)
      probs.curr.a<-pmin(pmax(exp(X%*%beta.curr.a)*baseline.prob,probs.min),1-probs.min)
      probs.curr.n<-pmin(pmax(baseline.prob,probs.min),1-probs.min)
      
      sums<-probs.curr.c+probs.curr.a+probs.curr.n
      probs.curr.c<-probs.curr.c/sums
      probs.curr.a<-probs.curr.a/sums
      probs.curr.n<-probs.curr.n/sums
      
      w.curr<-cbind(Z*Tr/(pZ*(probs.curr.c + probs.curr.a)) - 1, 
                    (1-Z)*Tr/((1-pZ)*probs.curr.a) - 1,
                    Z*(1-Tr)/(pZ*probs.curr.n) - 1, 
                    (1-Z)*(1-Tr)/((1-pZ)*(probs.curr.c + probs.curr.n)) - 1)
      
      w.curr.del<-1/n*t(X)%*%w.curr
      w.curr.del<-as.matrix(w.curr.del)
      w.curr<-as.matrix(w.curr)
      
      gbar<-c(1/n*t(X)%*%((Z*Tr/(1-probs.curr.n) + (1-Z)*(1-Tr)/(1-probs.curr.a) - 1)*probs.curr.c),
              1/n*t(X)%*%((Z*Tr/(1-probs.curr.n) + (1-Z)*Tr/probs.curr.a - 1)*probs.curr.a),
              w.curr.del)
      
      Ac<- -probs.curr.c*probs.curr.n/(probs.curr.c + probs.curr.a)^2
      Bc<- probs.curr.c/probs.curr.a
      Cc<- -probs.curr.c*probs.curr.a/(1-probs.curr.a)^2
      Dc<- probs.curr.c/probs.curr.n
      Aa<- -probs.curr.a*probs.curr.n/(probs.curr.c + probs.curr.a)^2
      Ba<- -(1-probs.curr.a)/probs.curr.a
      Ca<- probs.curr.a/(1 - probs.curr.a)
      Da<- probs.curr.a/probs.curr.n
      
      dgbar<-rbind(cbind(t(X*as.vector(probs.curr.c*(Z*Tr*Ac + (1-Z)*(1-Tr)*Cc + (Z*Tr/(probs.curr.c + probs.curr.a) + 
                                                                                    (1-Z)*(1-Tr)/(1-probs.curr.a) - 1)*(1 - probs.curr.c))))%*%X,
                         t(X*as.vector(probs.curr.a*(Z*Tr*Ac + (1-Z)*Tr*Bc - (Z*Tr/(probs.curr.c + probs.curr.a) + 
                                                                                (1-Z)*Tr/probs.curr.a - 1)*probs.curr.c)))%*%X,
                         t(X*as.vector(Z*Tr/pZ*Ac))%*%X,
                         t(X*as.vector((1-Z)*Tr/(1-pZ)*Bc))%*%X,
                         t(X*as.vector(Z*(1-Tr)/pZ*Dc))%*%X,
                         t(X*as.vector((1-Z)*(1-Tr)/(1-pZ)*Cc))%*%X),
                   cbind(t(X*as.vector(probs.curr.c*(Z*Tr*Aa + (1-Z)*(1-Tr)*Ca - (Z*Tr/(probs.curr.c + probs.curr.a) + 
                                                                                    (1-Z)*(1-Tr)/(1-probs.curr.a) - 1)*probs.curr.a)))%*%X,
                         t(X*as.vector(probs.curr.a*(Z*Tr*Aa + (1-Z)*Tr*Ba + (Z*Tr/(probs.curr.c + probs.curr.a) + 
                                                                                (1-Z)*Tr/probs.curr.a - 1)*(1-probs.curr.a))))%*%X,
                         t(X*as.vector(Z*Tr/pZ*Aa))%*%X,
                         t(X*as.vector((1-Z)*Tr/(1-pZ)*Ba))%*%X,
                         t(X*as.vector(Z*(1-Tr)/pZ*Da))%*%X,
                         t(X*as.vector((1-Z)*(1-Tr)/(1-pZ)*Ca))%*%X))/n
    }
    else{
      probs.curr <- pmin(pmax(exp(X%*%beta.curr)/(1 + exp(X%*%beta.curr)), probs.min), 1-probs.min)
      
      w.curr <- cbind(Z*Tr/(pZ * probs.curr) - 1, Z*(1-Tr)/(pZ *(1 - probs.curr)) - 1)
      w.curr.del<-1/n*t(X)%*%w.curr
      w.curr.del<-as.matrix(w.curr.del)
      w.curr<-as.matrix(w.curr)
      
      gbar<-c(1/n*t(X)%*%(Tr*Z*(1-probs.curr) - Z*(1-Tr)*probs.curr),
              w.curr.del)
      
      dgbar<-cbind(t(X*as.vector((-Z*Tr - Z*(1-Tr))*probs.curr*(1-probs.curr)))%*%X,
                   t(X*as.vector(-Z*Tr*(1-probs.curr)/(pZ*probs.curr)))%*%X,
                   t(X*as.vector(Z*(1-Tr)*probs.curr/(pZ*(1-probs.curr))))%*%X)/n
    }
    
    out<-2*dgbar%*%invV%*%gbar
    out
  }
  
  mle.loss <- function(beta.curr, twosided)
  {
    if (twosided){
      beta.curr.c<-beta.curr[1:k]
      beta.curr.a<-beta.curr[k+(1:k)]
      
      baseline.prob<-(1 + exp(X%*%beta.curr.c) + exp(X%*%beta.curr.a))^-1
      
      probs.curr.c<-pmin(pmax(exp(X%*%beta.curr.c)*baseline.prob,probs.min),1-probs.min)
      probs.curr.a<-pmin(pmax(exp(X%*%beta.curr.a)*baseline.prob,probs.min),1-probs.min)
      probs.curr.n<-pmin(pmax(1-probs.curr.c-probs.curr.a,probs.min),1-probs.min)
      
      sums<-probs.curr.c+probs.curr.a+probs.curr.n
      probs.curr.c<-probs.curr.c/sums
      probs.curr.a<-probs.curr.a/sums
      probs.curr.n<-probs.curr.n/sums
      
      # Take the negative because we are minimizing.   Want to minimize negative log-likelihood.
      loss<- -sum(Z*Tr*log(probs.curr.c+probs.curr.a) + Z*(1-Tr)*log(probs.curr.n) + (1-Z)*Tr*log(probs.curr.a) + (1-Z)*(1-Tr)*log(1-probs.curr.a))
    }
    else{
      probs.curr <- pmin(pmax(exp(X%*%beta.curr)/(1 + exp(X%*%beta.curr)),probs.min),1-probs.min)
      loss <- -sum(Z*Tr*log(probs.curr) + Z*(1 - Tr)*log(1 - probs.curr))
    }
    loss
  }
  
  mle.gradient <- function(beta.curr, twosided)
  {
    if (twosided){
      beta.curr.c<-beta.curr[1:k]
      beta.curr.a<-beta.curr[k+(1:k)]
      
      baseline.prob<-(1 + exp(X%*%beta.curr.c) + exp(X%*%beta.curr.a))^-1
      
      probs.curr.c<-pmin(pmax(exp(X%*%beta.curr.c)*baseline.prob,probs.min),1-probs.min)
      probs.curr.a<-pmin(pmax(exp(X%*%beta.curr.a)*baseline.prob,probs.min),1-probs.min)
      probs.curr.n<-pmin(pmax(baseline.prob,probs.min),1-probs.min)
      
      sums<-probs.curr.c+probs.curr.a+probs.curr.n
      probs.curr.c<-probs.curr.c/sums
      probs.curr.a<-probs.curr.a/sums
      probs.curr.n<-probs.curr.n/sums
      
      ds<- -c(t(X)%*%((Z*Tr/(probs.curr.c + probs.curr.a) + (1-Z)*(1-Tr)/(1 - probs.curr.a) - 1)*probs.curr.c),
              t(X)%*%((Z*Tr/(probs.curr.c + probs.curr.a) + (1-Z)*Tr/probs.curr.a - 1)*probs.curr.a))
    }
    else{
      probs.curr<-pmin(pmax(exp(X%*%beta.curr)/(1 + exp(X%*%beta.curr)),probs.min),1-probs.min)
      ds <- -t(X)%*%((Z*Tr/probs.curr - Z*(1-Tr)/(1-probs.curr))*probs.curr*(1-probs.curr))
    }
    ds
  }
  
  bal.loss <- function(beta.curr, invV, twosided)
  {
    if (twosided){
      beta.curr.c<-beta.curr[1:k]
      beta.curr.a<-beta.curr[k+(1:k)]
      
      baseline.prob<-(1 + exp(X%*%beta.curr.c) + exp(X%*%beta.curr.a))^-1
      
      probs.curr.c<-pmin(pmax(exp(X%*%beta.curr.c)*baseline.prob,probs.min),1-probs.min)
      probs.curr.a<-pmin(pmax(exp(X%*%beta.curr.a)*baseline.prob,probs.min),1-probs.min)
      probs.curr.n<-pmin(pmax(baseline.prob,probs.min),1-probs.min)
      
      sums<-probs.curr.c+probs.curr.a+probs.curr.n
      probs.curr.c<-probs.curr.c/sums
      probs.curr.a<-probs.curr.a/sums
      probs.curr.n<-probs.curr.n/sums
      
      w.curr<-cbind(Z*Tr/(pZ*(probs.curr.c + probs.curr.a)) - 1, 
                    (1-Z)*Tr/((1-pZ)*probs.curr.a) - 1,
                    Z*(1-Tr)/(pZ*probs.curr.n) - 1, 
                    (1-Z)*(1-Tr)/((1-pZ)*(probs.curr.c + probs.curr.n)) - 1)
      
      invV <- invV[(2*ncol(X)+1):ncol(invV),(2*ncol(X)+1):ncol(invV)]
    }
    else{
      probs.curr<-pmin(pmax(exp(X%*%beta.curr)/(1 + exp(X%*%beta.curr)),probs.min),1-probs.min)
      
      w.curr<-cbind(Z*Tr/(pZ * probs.curr) - 1, Z*(1-Tr)/(pZ *(1 - probs.curr)) - 1)
      invV <- invV[(ncol(X)+1):ncol(invV),(ncol(X)+1):ncol(invV)]
    }
    
    w.curr.del<-as.matrix(1/n*t(X)%*%w.curr)
    wbar <- c(w.curr.del)
    
    loss1<-sum(diag(t(wbar)%*%invV%*%wbar))
    loss1
  }
  
  bal.gradient <- function(beta.curr, invV, twosided)
  {
    if(twosided){
      beta.curr.c<-beta.curr[1:k]
      beta.curr.a<-beta.curr[k+(1:k)]
      
      baseline.prob<-(1 + exp(X%*%beta.curr.c) + exp(X%*%beta.curr.a))^-1
      
      probs.curr.c<-pmin(pmax(exp(X%*%beta.curr.c)*baseline.prob,probs.min),1-probs.min)
      probs.curr.a<-pmin(pmax(exp(X%*%beta.curr.a)*baseline.prob,probs.min),1-probs.min)
      probs.curr.n<-pmin(pmax(baseline.prob,probs.min),1-probs.min)
      
      sums<-probs.curr.c+probs.curr.a+probs.curr.n
      probs.curr.c<-probs.curr.c/sums
      probs.curr.a<-probs.curr.a/sums
      probs.curr.n<-probs.curr.n/sums
      
      Ac<- -probs.curr.c*probs.curr.n/(probs.curr.c + probs.curr.a)^2
      Bc<- probs.curr.c/probs.curr.a
      Cc<- -probs.curr.c*probs.curr.a/(1-probs.curr.a)^2
      Dc<- probs.curr.c/probs.curr.n
      Aa<- -probs.curr.a*probs.curr.n/(probs.curr.c + probs.curr.a)^2
      Ba<- -(1-probs.curr.a)/probs.curr.a
      Ca<- probs.curr.a/(1 - probs.curr.a)
      Da<- probs.curr.a/probs.curr.n
      
      w.curr<-cbind(Z*Tr/(pZ*(probs.curr.c + probs.curr.a)) - 1, 
                    (1-Z)*Tr/((1-pZ)*probs.curr.a) - 1,
                    Z*(1-Tr)/(pZ*probs.curr.n) - 1, 
                    (1-Z)*(1-Tr)/((1-pZ)*(probs.curr.c + probs.curr.n)) - 1)
      
      w.curr.del<-as.matrix(1/n*t(X)%*%w.curr)
      wbar <- c(w.curr.del)
      
      dw.beta.c<-1/n*cbind(t(X*as.vector(Z*Tr/pZ*Ac))%*%X,
                           t(X*as.vector((1-Z)*Tr/(1-pZ)*Bc))%*%X,
                           t(X*as.vector(Z*(1-Tr)/pZ*Dc))%*%X,
                           t(X*as.vector((1-Z)*(1-Tr)/(1-pZ)*Cc))%*%X)
      
      dw.beta.a<-1/n*cbind(t(X*as.vector(Z*Tr/pZ*Aa))%*%X,
                           t(X*as.vector((1-Z)*Tr/(1-pZ)*Ba))%*%X,
                           t(X*as.vector(Z*(1-Tr)/pZ*Da))%*%X,
                           t(X*as.vector((1-Z)*(1-Tr)/(1-pZ)*Ca))%*%X)
      
      invV <- invV[(2*ncol(X)+1):ncol(invV),(2*ncol(X)+1):ncol(invV)]
      
      out.1<-2*dw.beta.c%*%invV%*%wbar
      out.2<-2*dw.beta.a%*%invV%*%wbar
      
      out<-c(out.1, out.2)
    }
    else{
      probs.curr<-pmin(pmax(exp(X%*%beta.curr)/(1 + exp(X%*%beta.curr)),probs.min),1-probs.min)
      
      w.curr<-cbind(Z*Tr/(pZ * probs.curr) - 1, Z*(1-Tr)/(pZ *(1 - probs.curr)) - 1)
      w.curr.del<-as.matrix(1/n*t(X)%*%w.curr)
      wbar <- c(w.curr.del)
      
      invV <- invV[(ncol(X)+1):ncol(invV),(ncol(X)+1):ncol(invV)]
      
      dw.beta <- 1/n*cbind(t(X*as.vector(-Z*Tr*(1-probs.curr)/(pZ*probs.curr)))%*%X,
                           t(X*as.vector(Z*(1-Tr)*probs.curr/(pZ*(1-probs.curr))))%*%X)
      
      out<-2*dw.beta%*%invV%*%wbar
    }
    out
  }
  
  # Get starting point for optim
  # This block needs to be separated for one-sided and two-sided
  if (twosided){
    beta.n0 <- coef(glm(I(1-Tr) ~ - 1 + X, family = binomial(logit), subset = which(Z == 1)))
    beta.a0 <- coef(glm(Tr ~ -1 + X, family = binomial(logit), subset = which(Z == 0)))

    p.hat.a0 <- pmin(pmax(exp(X%*%beta.a0)/(1 + exp(X%*%beta.a0) + exp(X%*%beta.n0)), probs.min), 1-probs.min)
    p.hat.n0 <- pmin(pmax(exp(X%*%beta.n0)/(1 + exp(X%*%beta.a0) + exp(X%*%beta.n0)), probs.min), 1-probs.min)
    p.hat.c0 <- pmin(pmax(1/(1 + exp(X%*%beta.a0) + exp(X%*%beta.n0)), probs.min), 1-probs.min)
    
    sums <- p.hat.c0 + p.hat.a0 + p.hat.n0
    p.hat.c0 <- p.hat.c0/sums
    p.hat.a0 <- p.hat.a0/sums
    p.hat.n0 <- p.hat.n0/sums
    
    beta.init <- c(coef(lm(log(p.hat.c0/(p.hat.n0)) ~ -1 + X)), coef(lm(log(p.hat.a0/(p.hat.n0)) ~ -1 + X)))
  }
  else{
    beta.init <- coef(glm(Tr ~ -1 + X, subset = which(Z == 1)))
  }

  # All optimization functions need a one-sided or two-sided option
  mle.opt<-optim(beta.init, mle.loss, control=list("maxit"=iterations), method = "BFGS", gr = mle.gradient, twosided = twosided)
  beta.mle<-mle.opt$par

  this.invV<-gmm.func(beta.mle, twosided = twosided)$invV
  
  if (score.only)   gmm.opt<-mle.opt
  else {
    bal.init.opt<-optim(beta.init, bal.loss, control=list("maxit"=iterations), method = "BFGS", invV = this.invV, gr = bal.gradient, twosided = twosided)
    bal.mle.opt<-optim(beta.mle, bal.loss, control=list("maxit"=iterations), method = "BFGS", invV = this.invV, gr = bal.gradient, twosided = twosided)
    if(bal.init.opt$value > bal.mle.opt$value){
      bal.opt <- bal.mle.opt
    }
    else{
      bal.opt <- bal.init.opt
    }
    
    beta.bal<-bal.opt$par
    
    if (bal.only) gmm.opt<-bal.opt
    else {
      gmm.mle.opt<-optim(beta.mle, gmm.loss, control=list("maxit"=iterations), method = "BFGS", invV = this.invV, gr = gmm.gradient, twosided = twosided)      
      gmm.bal.opt<-optim(beta.bal, gmm.loss, control=list("maxit"=iterations), method = "BFGS", invV = this.invV, gr = gmm.gradient, twosided = twosided)      
      if (gmm.mle.opt$value > gmm.bal.opt$value)
      {
        gmm.opt<-gmm.bal.opt
      }
      else
      {
        gmm.opt<-gmm.mle.opt
      }    
    }
  }
  
  beta.opt<-matrix(gmm.opt$par,nrow=k)
  J.opt<-gmm.loss(beta.opt, this.invV, twosided)
  bal.loss.opt <- bal.loss(beta.opt, invV = this.invV, twosided)

  class(beta.opt)<-"coef"
  
  if (twosided){
    pi.c.opt<-as.vector(exp(X%*%beta.opt[,1]))
    pi.a.opt<-as.vector(exp(X%*%beta.opt[,2]))
    pi.n.opt<-1
    sums<-pi.c.opt+pi.a.opt+pi.n.opt
    
    #Normalize, then trim
    pi.c.opt<-pmax(pmin(pi.c.opt/sums, 1-probs.min),probs.min)
    pi.a.opt<-pmax(pmin(pi.a.opt/sums, 1-probs.min),probs.min)
    pi.n.opt<-pmax(pmin(pi.n.opt/sums, 1-probs.min),probs.min)
    
    #Renormalize, so that they add up to 1
    sums<-pi.c.opt+pi.a.opt+pi.n.opt
    fitted.values <- cbind(pi.c.opt/sums, pi.a.opt/sums, pi.n.opt/sums)
    colnames(fitted.values)<-c("Compliers","Always","Never")  
    
    beta.opt[-1,]<-beta.opt[-1,]/x.sd
    
    deviance<- -2*sum(Z*Tr*log(fitted.values[,1]+fitted.values[,2]) + Z*(1-Tr)*log(fitted.values[,3]) + (1-Z)*Tr*log(fitted.values[,2]) + (1-Z)*(1-Tr)*log(1-fitted.values[,2]))
    
    if (k > 2)
    {
      beta.opt[1,]<-beta.opt[1,]-matrix(x.mean%*%beta.opt[-1,])
    }
    else
    {
      beta.opt[1,]<-beta.opt[1,]-x.mean*beta.opt[-1,]
    }
  }
  else{
    fitted.values<-pmax(pmin(exp(X%*%beta.opt)/(1+exp(X%*%beta.opt)),1-probs.min),probs.min)
    deviance<- -2*sum(Z*Tr*log(fitted.values) + Z*(1-Tr)*log(1-fitted.values))
    
    beta.opt[-1]<-beta.opt[-1]/x.sd
    beta.opt[1]<-beta.opt[1]-sum(x.mean*beta.opt[-1])
  }
  
  if (method == "exact" & twosided == TRUE){
    this.invV <- this.invV[(2*ncol(X)+1):ncol(this.invV),(2*ncol(X)+1):ncol(this.invV)]
  }
  
  output<-list("coefficients"=beta.opt,"fitted.values"=fitted.values,"weights"=1/fitted.values[,1],
               "deviance"=deviance,"converged"=gmm.opt$conv,"J"=J.opt,"df"=k,
               "bal"=bal.loss.opt, "x" = X.orig, "z" = Z, "y" = Tr, "invSigma" = this.invV)
  class(output)<-"CBIV"
  output
}

vcov_outcome.CBIV <- function(object, Y, Ttilde, Ztilde, delta){
  X <- object$x

  N <- length(Y)
  K <- ncol(X)
  P <- ncol(Ttilde)
  
  w <- object$weights
  beta <- object$coefficients
  Z <- object$z
  Tr <- object$y
  invSigma <- object$invSigma
  pZ <- mean(object$z)
  
  pi.min<-10^-6
  
  baseline.prob <- (1 + exp(X%*%beta[,1]) + exp(X%*%beta[,2]))^-1
  pi.c<-pmin(pmax(exp(X%*%beta[,1])*baseline.prob,pi.min),1-pi.min)
  pi.a<-pmin(pmax(exp(X%*%beta[,2])*baseline.prob,pi.min),1-pi.min)
  pi.n<-pmin(pmax(baseline.prob,pi.min),1-pi.min)
  
  sums<-pi.c+pi.a+pi.n
  pi.c<-pi.c/sums
  pi.a<-pi.a/sums
  pi.n<-pi.n/sums
  
  mtheta <- t(as.vector(w*(Y - Ttilde%*%delta))*Ztilde)
  
  gbeta <- rbind(t(as.vector((Z*Tr/(1-pi.n) + (1-Z)*(1-Tr)/(1-pi.a) - 1)*pi.c)*X),
                 t(as.vector((Z*Tr/(1-pi.n) + (1-Z)*Tr/pi.a - 1)*pi.a)*X),
                 t(as.vector(Z*Tr/(pZ*(pi.c + pi.a)) - 1)*X), 
                 t(as.vector((1-Z)*Tr/((1-pZ)*pi.a) - 1)*X),
                 t(as.vector(Z*(1-Tr)/(pZ*pi.n) - 1)*X), 
                 t(as.vector((1-Z)*(1-Tr)/((1-pZ)*(pi.c + pi.n)) - 1)*X))
  
  gtheta <- rbind(mtheta, gbeta)
  
  W <- rbind(cbind(diag(P), matrix(0, nrow = P, ncol = ncol(invSigma))), 
             cbind(matrix(0, nrow = nrow(invSigma), ncol = P), invSigma))

  Mdelta <- t(Ztilde)%*%(-w*Ttilde)/N
  
  dw <- cbind(as.vector(-(1-pi.c)/pi.c)*X,
              as.vector(pi.a/pi.c)*X)
  
  Mbeta <- t(Ztilde*as.vector(Y - Ttilde%*%delta))%*%dw/N
  
  G <- rbind(cbind(t(as.vector((Z*Tr*pi.a/(1-pi.n)^2 + (1-Z)*(1-Tr)*pi.n/(1-pi.a)^2 + pi.c - 1)*pi.c)*X)%*%X/N,
                   t(as.vector((1 - Z*Tr/(1-pi.n)^2)*pi.c*pi.a)*X)%*%X/N),
             cbind(t(as.vector((1 - Z*Tr/(1-pi.n)^2)*pi.c*pi.a)*X)%*%X/N,
                   t(as.vector((Z*Tr*pi.c/(1-pi.n)^2 + pi.a - 1)*pi.a)*X)%*%X/N),
             cbind(t(as.vector(-Z*Tr*pi.c*pi.n/(pZ*(1-pi.n)^2))*X)%*%X/N,
                   t(as.vector(-Z*Tr*pi.a*pi.n/(pZ*(1-pi.n)^2))*X)%*%X/N),
             cbind(t(as.vector((1-Z)*Tr*pi.c/((1-pZ)*pi.a))*X)%*%X/N,
                   t(as.vector(-(1-Z)*Tr*(1-pi.a)/((1-pZ)*pi.a))*X)%*%X/N),
             cbind(t(as.vector(Z*(1-Tr)*pi.c/(pZ*pi.n))*X)%*%X/N,
                   t(as.vector(Z*(1-Tr)*pi.a/(pZ*pi.n))*X)%*%X/N),
             cbind(t(as.vector(-(1-Z)*(1-Tr)*pi.c*pi.a/((1-pZ)*(1-pi.a)^2))*X)%*%X/N,
                   t(as.vector((1-Z)*(1-Tr)*pi.a/((1-pZ)*(1-pi.a)))*X)%*%X/N))
  
  Gtilde <- rbind(cbind(Mdelta, Mbeta), cbind(matrix(0, nrow = nrow(G), ncol = P), G))
  
  outer.inv <- t(Gtilde)%*%W%*%Gtilde
  
  #cond.num <- svd(outer.inv)$d[1]/svd(outer.inv)$d[nrow(outer.inv)]
  #if (cond.num>(1/tol)){outer.inv <- outer.inv+lambda*diag(rep(1,nrow(outer.inv)))}
  
  outer <- ginv(outer.inv)
  inner <- gtheta%*%t(gtheta)/N
  
  V <- outer %*% t(Gtilde) %*% W %*% inner %*% W %*% Gtilde %*% outer

  return(V/N)
}

wols_vcov_outcome.CBIV <- function(object, Y, Ztilde, delta){
  X <- object$x
  
  N <- length(Y)
  K <- ncol(X)
  P <- ncol(Ztilde)
  
  w <- object$weights
  beta <- object$coefficients
  Z <- object$z
  Tr <- object$y
  invSigma <- object$invSigma
  pZ <- mean(object$z)
  
  Zind <- which(apply(Ztilde, 2, function(u) all(u == Z)))
  Intind <- which(apply(Ztilde, 2, function(u) all(u == 1)))
  Xtilde <- Ztilde[,-c(Zind,Intind),drop=FALSE]
  deltaXtilde <- delta[-c(Zind,Intind)]
  if (length(Xtilde) > 0){
    Ztilde[,-c(Zind,Intind)] <- Xtilde*w
  }
  
  pi.min<-10^-6
  
  baseline.prob <- (1 + exp(X%*%beta[,1]) + exp(X%*%beta[,2]))^-1
  pi.c<-pmin(pmax(exp(X%*%beta[,1])*baseline.prob,pi.min),1-pi.min)
  pi.a<-pmin(pmax(exp(X%*%beta[,2])*baseline.prob,pi.min),1-pi.min)
  pi.n<-pmin(pmax(baseline.prob,pi.min),1-pi.min)
  
  sums<-pi.c+pi.a+pi.n
  pi.c<-pi.c/sums
  pi.a<-pi.a/sums
  pi.n<-pi.n/sums
  
  mtheta <- t(as.vector((w*Y - Ztilde%*%delta))*Ztilde)
  
  gbeta <- rbind(t(as.vector((Z*Tr/(1-pi.n) + (1-Z)*(1-Tr)/(1-pi.a) - 1)*pi.c)*X),
                 t(as.vector((Z*Tr/(1-pi.n) + (1-Z)*Tr/pi.a - 1)*pi.a)*X),
                 t(as.vector(Z*Tr/(pZ*(pi.c + pi.a)) - 1)*X), 
                 t(as.vector((1-Z)*Tr/((1-pZ)*pi.a) - 1)*X),
                 t(as.vector(Z*(1-Tr)/(pZ*pi.n) - 1)*X), 
                 t(as.vector((1-Z)*(1-Tr)/((1-pZ)*(pi.c + pi.n)) - 1)*X))
  
  gtheta <- rbind(mtheta, gbeta)
  
  W <- rbind(cbind(diag(P), matrix(0, nrow = P, ncol = ncol(invSigma))), 
             cbind(matrix(0, nrow = nrow(invSigma), ncol = P), invSigma))
  
  Mdelta <- t(-Ztilde)%*%Ztilde/N
  
  dw <- cbind(as.vector(-(1-pi.c)/pi.c)*X,
              as.vector(pi.a/pi.c)*X)

  Mbeta <- matrix(0, nrow = P, ncol = 2*K)
  Mbeta[Intind,] <- t(as.vector(Y - Xtilde%*%deltaXtilde))%*%dw/N
  Mbeta[Zind,] <- t(Z*(Y - Xtilde%*%deltaXtilde))%*%dw/N
  if (length(Xtilde) > 0){
    Mbeta[-c(Intind,Zind),] <- t(Xtilde*as.vector(w*Y - Ztilde%*%delta) + 
                                 w*Xtilde*as.vector(Y - Xtilde%*%deltaXtilde))%*%dw/N   
  }

  #Mbeta <- (t(Ztilde*as.vector(Y))%*%dw)/N
  
  G <- rbind(cbind(t(as.vector((Z*Tr*pi.a/(1-pi.n)^2 + (1-Z)*(1-Tr)*pi.n/(1-pi.a)^2 + pi.c - 1)*pi.c)*X)%*%X/N,
                   t(as.vector((1 - Z*Tr/(1-pi.n)^2)*pi.c*pi.a)*X)%*%X/N),
             cbind(t(as.vector((1 - Z*Tr/(1-pi.n)^2)*pi.c*pi.a)*X)%*%X/N,
                   t(as.vector((Z*Tr*pi.c/(1-pi.n)^2 + pi.a - 1)*pi.a)*X)%*%X/N),
             cbind(t(as.vector(-Z*Tr*pi.c*pi.n/(pZ*(1-pi.n)^2))*X)%*%X/N,
                   t(as.vector(-Z*Tr*pi.a*pi.n/(pZ*(1-pi.n)^2))*X)%*%X/N),
             cbind(t(as.vector((1-Z)*Tr*pi.c/((1-pZ)*pi.a))*X)%*%X/N,
                   t(as.vector(-(1-Z)*Tr*(1-pi.a)/((1-pZ)*pi.a))*X)%*%X/N),
             cbind(t(as.vector(Z*(1-Tr)*pi.c/(pZ*pi.n))*X)%*%X/N,
                   t(as.vector(Z*(1-Tr)*pi.a/(pZ*pi.n))*X)%*%X/N),
             cbind(t(as.vector(-(1-Z)*(1-Tr)*pi.c*pi.a/((1-pZ)*(1-pi.a)^2))*X)%*%X/N,
                   t(as.vector((1-Z)*(1-Tr)*pi.a/((1-pZ)*(1-pi.a)))*X)%*%X/N))
  
  Gtilde <- rbind(cbind(Mdelta, Mbeta), cbind(matrix(0, nrow = nrow(G), ncol = P), G))
  
  outer.inv <- t(Gtilde)%*%W%*%Gtilde
  outer <- ginv(outer.inv)
  inner <- gtheta%*%t(gtheta)/N
  
  V <- outer %*% t(Gtilde) %*% W %*% inner %*% W %*% Gtilde %*% outer
  
  return(V/N)
}

ht_vcov_outcome.CBIV <- function(object, Y){
  X <- object$x
  
  N <- length(Y)
    
  beta <- object$coefficients
  Z <- object$z
  Tr <- object$y
  invSigma <- object$invSigma
  pZ <- mean(object$z)
  
  pi.min<-10^-6
  
  baseline.prob <- (1 + exp(X%*%beta[,1]) + exp(X%*%beta[,2]))^-1
  pi.c<-pmin(pmax(exp(X%*%beta[,1])*baseline.prob,pi.min),1-pi.min)
  pi.a<-pmin(pmax(exp(X%*%beta[,2])*baseline.prob,pi.min),1-pi.min)
  pi.n<-pmin(pmax(baseline.prob,pi.min),1-pi.min)
  
  sums<-pi.c+pi.a+pi.n
  pi.c<-pi.c/sums
  pi.a<-pi.a/sums
  pi.n<-pi.n/sums
  
  gbeta <- rbind(t(as.vector((Z*Tr/(1-pi.n) + (1-Z)*(1-Tr)/(1-pi.a) - 1)*pi.c)*X),
                 t(as.vector((Z*Tr/(1-pi.n) + (1-Z)*Tr/pi.a - 1)*pi.a)*X),
                 t(as.vector(Z*Tr/(pZ*(pi.c + pi.a)) - 1)*X), 
                 t(as.vector((1-Z)*Tr/((1-pZ)*pi.a) - 1)*X),
                 t(as.vector(Z*(1-Tr)/(pZ*pi.n) - 1)*X), 
                 t(as.vector((1-Z)*(1-Tr)/((1-pZ)*(pi.c + pi.n)) - 1)*X))

  G <- rbind(cbind(t(as.vector((Z*Tr*pi.a/(1-pi.n)^2 + (1-Z)*(1-Tr)*pi.n/(1-pi.a)^2 + pi.c - 1)*pi.c)*X)%*%X/N,
                   t(as.vector((1 - Z*Tr/(1-pi.n)^2)*pi.c*pi.a)*X)%*%X/N),
             cbind(t(as.vector((1 - Z*Tr/(1-pi.n)^2)*pi.c*pi.a)*X)%*%X/N,
                   t(as.vector((Z*Tr*pi.c/(1-pi.n)^2 + pi.a - 1)*pi.a)*X)%*%X/N),
             cbind(t(as.vector(-Z*Tr*pi.c*pi.n/(pZ*(1-pi.n)^2))*X)%*%X/N,
                   t(as.vector(-Z*Tr*pi.a*pi.n/(pZ*(1-pi.n)^2))*X)%*%X/N),
             cbind(t(as.vector((1-Z)*Tr*pi.c/((1-pZ)*pi.a))*X)%*%X/N,
                   t(as.vector(-(1-Z)*Tr*(1-pi.a)/((1-pZ)*pi.a))*X)%*%X/N),
             cbind(t(as.vector(Z*(1-Tr)*pi.c/(pZ*pi.n))*X)%*%X/N,
                   t(as.vector(Z*(1-Tr)*pi.a/(pZ*pi.n))*X)%*%X/N),
             cbind(t(as.vector(-(1-Z)*(1-Tr)*pi.c*pi.a/((1-pZ)*(1-pi.a)^2))*X)%*%X/N,
                   t(as.vector((1-Z)*(1-Tr)*pi.a/((1-pZ)*(1-pi.a)))*X)%*%X/N))
  
  outer.inv <- t(G)%*%invSigma%*%G
  
  outer <- ginv(outer.inv)
  inner <- gbeta%*%t(gbeta)/N
  
  U <- outer %*% t(G) %*% invSigma %*% inner %*% invSigma %*% G %*% outer
  
  dw <- cbind(as.vector(-(1-pi.c)/pi.c)*X,
              as.vector(pi.a/pi.c)*X)
  V <- apply(as.vector(Y*Z/pZ - Y*(1-Z)/(1-pZ))*dw, 2, mean)

  return(t(V)%*%U%*%V/N)
}