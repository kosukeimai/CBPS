#' hdCBPS: high dimensional CBPS method
#'
#' hdCBPS calculates ATE using CBPS method in a high dimensional setting.
#' @param formula	An object of class formula (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param data An optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which CBPS is called.
#' @param na.action	A function which indicates what should happen when the data contain NAs. The default is set by the na.action setting of options, and is na.fail if that is unset.
#' @param iterations An optional parameter for the maximum number of iterations for the optimization. Default is 1000.
#' @param method Choose among "linear", "binomial", and "possion".
#' @param y An outcome variable.
#' @return
#' \item{ATT}{Average treatment effect on the treated.}
#' \item{ATE}{Average treatment effect.}
#' \item{s}{Standard Error.}
#' \item{fitted.values}{The fitted propensity score}
#' \item{coefficients1}{Coefficients for the treated propensity score}
#' \item{coefficients0}{Coefficients for the untreated propensity score}
#' \item{model}{The model frame}
#' @author Sida Peng
#' @export


# hdCBPS parses the formula object and passes the result to hdCBPS.fit


#' hdCBPS: high dimensional CBPS method
#' 
#' hdCBPS calculates ATE using CBPS method in a high dimensional setting.
#' 
#' 
#' @aliases hdCBPS hdCBPS
#' @param formula An object of class formula (or one that can be coerced to
#' that class): a symbolic description of the model to be fitted.
#' @param data An optional data frame, list or environment (or object coercible
#' by as.data.frame to a data frame) containing the variables in the model. If
#' not found in data, the variables are taken from environment(formula),
#' typically the environment from which CBPS is called.
#' @param na.action A function which indicates what should happen when the data
#' contain NAs. The default is set by the na.action setting of options, and is
#' na.fail if that is unset.
#' @param y An outcome variable.
#' @param ATT Option to calculate ATT
#' @param iterations An optional parameter for the maximum number of iterations
#' for the optimization. Default is 1000.
#' @param method Choose among "linear", "binomial", and "possion".
#' @return \item{ATT}{Average treatment effect on the treated.}
#' \item{ATE}{Average treatment effect.} \item{s}{Standard Error.}
#' \item{fitted.values}{The fitted propensity score}
#' \item{coefficients1}{Coefficients for the treated propensity score}
#' \item{coefficients0}{Coefficients for the untreated propensity score}
#' \item{model}{The model frame}
#' @author Sida Peng
hdCBPS <- function(formula, data, na.action, y, ATT = 0, iterations=1000, method="linear") {
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
  T <- model.response(mf, "any")
  if (length(dim(T)) == 1L) {
    nm <- rownames(T)
    dim(T) <- NULL
    if (!is.null(nm))
      names(T) <- nm
  }
  
  X <- if (!is.empty.model(mt)) model.matrix(mt, mf)#[,-2]
  else matrix(, NROW(T), 0L)
  
  X<-cbind(1,X[,apply(X,2,sd)>0])
  
  fit <- eval(call("hdCBPS.fit", x = X, y = y,  treat = T, Att = ATT, iterations=iterations, methd=method))
  
  fit$na.action <- attr(mf, "na.action")
  xlevels <- .getXlevels(mt, mf)
  fit$data<-data
  fit$call <- call
  fit$formula <- formula
  fit$terms<-mt
  fit
}




hdCBPS.fit <- function(x,y,treat, Att, iterations=1000, methd="linear") {
  n = dim(x)[1]
  p = dim(x)[2]
  
  y1hat = y[treat==1]
  x1hat = x[treat==1,]
  
  y0hat = y[treat==0]
  x0hat = x[treat==0,]
  
  if (methd =="linear"){
    cov1  = cv.glmnet(x1hat,y1hat)
    cov0  = cv.glmnet(x0hat,y0hat)
  } else if (methd=="binomial") {
    cov1  = cv.glmnet(x1hat,y1hat,family = "binomial",intercept=FALSE)
    cov0  = cv.glmnet(x0hat,y0hat,family = "binomial",intercept=FALSE)
  } else if (methd=="poisson") {
    cov1  = cv.glmnet(x1hat,y1hat,family = "poisson",intercept=FALSE)
    cov0  = cv.glmnet(x0hat,y0hat,family = "poisson",intercept=FALSE)
  }
  
  covb  = cv.glmnet(x,treat,family="binomial")
  
  S1 = which(coef(cov1)!=0)
  S0 = which(coef(cov0)!=0)
  
  
  ##Generates ATE weights.	Called by loss function, etc.
  ATE.wt.func<-function(beta.curr, S, tt, X.wt, beta.ini= coef(covb)){
    x2<-as.matrix(X.wt)
    n2<-dim(x2)[1]
    X2 = cbind(rep(1,n2),x2)
    beta.all = beta.ini
    beta.all[S] = beta.curr
    theta.curr<-as.vector(X2%*%beta.all)
    probs.curr<-1-(1+exp(theta.curr))^-1
    if (tt == 0){
      W<- (1-((1-treat)/(1-probs.curr)))
    }else{
      W<- (treat/probs.curr-1)
    }
    out<-list("W"=W)
    out
  }
  
  ##Generates ATE weights nolinear.	Called by loss function, etc.
  ATE.wt.nl.func<-function(beta.curr, S, tt, X.wt, beta.ini= coef(covb)){
    x2<-as.matrix(X.wt)
    n2<-dim(x2)[1]
    X2 = cbind(rep(1,n2),x2)
    SS = c(1,S)
    beta.all = beta.ini
    beta.all[SS] = beta.curr
    #beta.all = as.vector(beta.all)
    theta.curr<-as.vector(X2%*%beta.all)
    probs.curr<-1-(1+exp(theta.curr))^-1
    if (tt == 0){
      W<- (1-((1-treat)/(1-probs.curr)))
    }else{
      W<- (treat/probs.curr-1)
    }
    out<-list("W"=W)
    out
  }
  
  ##Generates ATT weights.	Called by loss function, etc.
  ATT.wt.func<-function(beta.curr, S, X.wt, beta.ini= coef(covb)){
    x2<-as.matrix(X.wt)
    n2<-dim(x2)[1]
    X2 = cbind(rep(1,n2),x2)
    beta.all = beta.ini
    beta.all[S] = beta.curr
    #beta.all = as.vector(beta.all)
    theta.curr<-as.vector(X2%*%beta.all)
    probs.curr<-1-(1+exp(theta.curr))^-1
    
    W<- (treat-(((1-treat)*probs.curr)/(1-probs.curr)))
    
    out<-list("W"=W)
    out
  }
  
  ##Generates ATT weights nolinear.	Called by loss function, etc.
  ATT.wt.nl.func<-function(beta.curr, S, X.wt, beta.ini= coef(covb)){
    x2<-as.matrix(X.wt)
    n2<-dim(x2)[1]
    X2 = cbind(rep(1,n2),x2)
    SS = c(1,S)
    beta.all = beta.ini
    beta.all[SS] = beta.curr
    #beta.all = as.vector(beta.all)
    theta.curr<-as.vector(X2%*%beta.all)
    probs.curr<-1-(1+exp(theta.curr))^-1
    
    W<- (treat-(((1-treat)*probs.curr)/(1-probs.curr)))
    
    out<-list("W"=W)
    out
  }
  
  
  gmm.func<-function(beta.curr, S, tt, X.gmm, methd){
    ##Designate a few objects in the function.
    x1<-as.matrix(X.gmm)
    n1<-dim(x1)[1]
    
    ##Generate the vector of mean imbalance by weights.
    if (methd =="linear"){
      w.curr<-ATE.wt.func(beta.curr, S, tt, x1)
      X1 = cbind(rep(1,n1),x1)
      if (length(S) != 0){
        w.curr.del<- t(X1[,S])%*%(w.curr$W)
        w.curr.del<- as.vector(w.curr.del)
      }else{
        w.curr.del = matrix(0, nrow = n1, ncol = 0)
      }
    }else if (methd=="poisson") {
      w.curr<-ATE.wt.nl.func(beta.curr, S, tt, x1)
      X1 = cbind(rep(1,n1),x1)
      if (tt==1){
        pweight = exp(as.vector(X1%*%coef(cov1)))
      }else{
        pweight = exp(as.vector(X1%*%coef(cov0)))
      }
      ##Generate the vector of mean imbalance by weights.
      if (length(S) != 0){
        w.curr.del<- t(cbind(as.matrix(pweight), matrix(rep(pweight, length(S)),ncol=length(S))*X1[,S]))%*%(w.curr$W)
        w.curr.del<- as.vector(w.curr.del)
      }else{
        w.curr.del = t(as.matrix(pweight))%*%(w.curr$W)
        w.curr.del<- as.vector(w.curr.del)
      }
    }else if (methd=="binomial") {
      w.curr<-ATE.wt.nl.func(beta.curr, S, tt, x1)
      X1 = cbind(rep(1,n1),x1)
      if (tt==1){
        pweight1 = exp(as.vector(X1%*%coef(cov1)))/(1+exp(as.vector(X1%*%coef(cov1))))
        pweight2 = exp(as.vector(X1%*%coef(cov1)))/(1+exp(as.vector(X1%*%coef(cov1))))^2
      }else{
        pweight1 = exp(as.vector(X1%*%coef(cov0)))/(1+exp(as.vector(X1%*%coef(cov0))))
        pweight2 = exp(as.vector(X1%*%coef(cov0)))/(1+exp(as.vector(X1%*%coef(cov0))))^2
      }
      
      ##Generate the vector of mean imbalance by weights.
      if (length(S) != 0){
        w.curr.del<- t(cbind(as.matrix(pweight1), matrix(rep(pweight2, length(S)),ncol=length(S))*X1[,S]))%*%(w.curr$W)
        w.curr.del<- as.vector(w.curr.del)
      }else{
        w.curr.del = t(as.matrix(pweight1))%*%(w.curr$W)
        w.curr.del<- as.vector(w.curr.del)
      }
    }
    
    
    
    ##Generate g-bar, as in the paper.
    gbar<- w.curr.del
    
    ##Calculate the GMM loss.
    loss1<-as.vector(t(gbar)%*%(gbar))
    out1<-list("loss"=loss1)
    out1
  }
  
  ATT.gmm.func<-function(beta.curr, S, X.gmm, methd){
    ##Designate a few objects in the function.
    x1<-as.matrix(X.gmm)
    n1<-dim(x1)[1]
    
    ##Generate the vector of mean imbalance by weights.
    if (methd =="linear"){
      w.curr<-ATT.wt.func(beta.curr, S, x1)
      X1 = cbind(rep(1,n1),x1)
      if (length(S) != 0){
        w.curr.del<- t(X1[,S])%*%(w.curr$W)
        w.curr.del<- as.vector(w.curr.del)
      }else{
        w.curr.del = matrix(0, nrow = n1, ncol = 0)
      }
    }else if (methd=="poisson") {
      w.curr<-ATT.wt.nl.func(beta.curr, S, x1)
      X1 = cbind(rep(1,n1),x1)
      pweight = exp(as.vector(X1%*%coef(cov0)))
      
      ##Generate the vector of mean imbalance by weights.
      if (length(S) != 0){
        w.curr.del<- t(cbind(as.matrix(pweight), matrix(rep(pweight, length(S)),ncol=length(S))*X1[,S]))%*%(w.curr$W)
        w.curr.del<- as.vector(w.curr.del)
      }else{
        w.curr.del = t(as.matrix(pweight))%*%(w.curr$W)
        w.curr.del<- as.vector(w.curr.del)
      }
    }else if (methd=="binomial") {
      w.curr<-ATT.wt.nl.func(beta.curr, S, x1)
      X1 = cbind(rep(1,n1),x1)
      pweight1 = exp(as.vector(X1%*%coef(cov0)))/(1+exp(as.vector(X1%*%coef(cov0))))
      pweight2 = exp(as.vector(X1%*%coef(cov0)))/(1+exp(as.vector(X1%*%coef(cov0))))^2
      
      
      ##Generate the vector of mean imbalance by weights.
      if (length(S) != 0){
        w.curr.del<- t(cbind(as.matrix(pweight1), matrix(rep(pweight2, length(S)),ncol=length(S))*X1[,S]))%*%(w.curr$W)
        w.curr.del<- as.vector(w.curr.del)
      }else{
        w.curr.del = t(as.matrix(pweight1))%*%(w.curr$W)
        w.curr.del<- as.vector(w.curr.del)
      }
    }
    
    
    
    ##Generate g-bar, as in the paper.
    gbar<- w.curr.del
    
    ##Calculate the GMM loss.
    loss1<-as.vector(t(gbar)%*%(gbar))
    out1<-list("loss"=loss1)
    out1
  }
  
  
  
  
  tol = 1e-5
  kk0 =0
  kk1 =0
  ATT.kk0 =0
  diff0 = 1
  diff1 = 1
  ATT.diff0 = 1 
  gmm.loss0<-function(xx,...) gmm.func(xx, S=S0, tt=0, X.gmm = x, methd = methd)$loss
  gmm.loss1<-function(xx,...) gmm.func(xx, S=S1, tt=1, X.gmm = x, methd = methd)$loss
  if (methd =="linear"){
    
    beta0 = optim(coef(covb)[S0], gmm.loss0, method = "Nelder-Mead")
    while (diff0>tol & kk0<iterations){
      beta0 = optim(beta0$par, gmm.loss0, method = "Nelder-Mead")
      diff0 = beta0$value
      kk0 = kk0+1
    }
    
    beta1 = optim(coef(covb)[S1], gmm.loss1, method = "Nelder-Mead")
    while (diff1>tol & kk1<iterations){
      beta1 = optim(beta1$par, gmm.loss1, method = "Nelder-Mead")
      diff1 = beta1$value
      kk1 = kk1+1
    }
    w.curr1<-ATE.wt.func(beta1$par, S1, 1, x)
    w.curr0<-ATE.wt.func(beta0$par, S0, 0, x)
  } else{
    
    beta0 = optim(coef(covb)[c(1,S0)], gmm.loss0, method = "Nelder-Mead")
    while (diff0>tol & kk0<iterations){
      beta0 = optim(beta0$par, gmm.loss0, method = "Nelder-Mead")
      diff0 = beta0$value
      kk0 = kk0+1
    }
    
    beta1 = optim(coef(covb)[c(1,S1)], gmm.loss1, method = "Nelder-Mead")
    while (diff1>tol & kk1<iterations){
      beta1 = optim(beta1$par, gmm.loss1, method = "Nelder-Mead")
      diff1 = beta1$value
      kk1 = kk1+1
    }
    w.curr1<-ATE.wt.nl.func(beta1$par, S1, 1, x)
    w.curr0<-ATE.wt.nl.func(beta0$par, S0, 0, x)
  }
  
  ATE = 1/(n)*(t(y1hat)%*%(w.curr1$W[treat==1]+1)+t(y0hat)%*%(w.curr0$W[treat==0]-1))
  ATT = NULL
  w = NULL
  
  if (Att==1){
    if (methd =="linear"){
      ATT.gmm.loss<-function(xx,...) ATT.gmm.func(xx, S=S0,  X.gmm = x, methd = methd)$loss
      #beta0 = optim(coef(covb)[S0], gmm.loss0,method="BFGS")
      ATT.beta0 = optim(coef(covb)[S0], ATT.gmm.loss, method = "Nelder-Mead")
      while (ATT.diff0>tol & ATT.kk0<iterations){
        ATT.beta0 = optim(ATT.beta0$par, ATT.gmm.loss, method = "Nelder-Mead")
        ATT.diff0 = ATT.beta0$value
        ATT.kk0 = ATT.kk0+1
      }
      
      ATT.w.curr0<-ATT.wt.func(ATT.beta0$par, S0, x)
      X = cbind(rep(1,n),x)
      ATT.beta.0 = coef(covb)
      ATT.beta.0[S0] = ATT.beta0$par
      ATT.beta.0 = as.matrix(ATT.beta.0)
      ATT.theta.0<-as.vector(X%*%ATT.beta.0)
      ATT.r_yhatb0<-1-(1+exp(ATT.theta.0))^-1
      
    }else{
      ATT.gmm.loss<-function(xx,...) ATT.gmm.func(xx, S=S0,  X.gmm = x, methd = methd)$loss
      #beta0 = optim(coef(covb)[S0], gmm.loss0,method="BFGS")
      ATT.beta0 = optim(coef(covb)[c(1,S0)], ATT.gmm.loss, method = "Nelder-Mead")
      while (ATT.diff0>tol & ATT.kk0<iterations){
        ATT.beta0 = optim(ATT.beta0$par, ATT.gmm.loss, method = "Nelder-Mead")
        ATT.diff0 = ATT.beta0$value
        ATT.kk0 = ATT.kk0+1
      }
      
      ATT.w.curr0<-ATT.wt.nl.func(ATT.beta0$par, S0, x)
      X = cbind(rep(1,n),x)
      ATT.beta.0 = coef(covb)
      ATT.beta.0[c(1,S0)] = ATT.beta0$par
      ATT.beta.0 = as.matrix(ATT.beta.0)
      ATT.theta.0<-as.vector(X%*%ATT.beta.0)
      ATT.r_yhatb0<-1-(1+exp(ATT.theta.0))^-1
    }
    ATT = 1/(sum(treat))*sum(y1hat)-1/sum(ATT.w.curr0$W[treat==0])*(t(y0hat)%*%(ATT.w.curr0$W[treat==0]))
  }
  
  
  
  
  
  
  if (methd =="linear"){
    
    X = cbind(rep(1,n),x)
    beta.0 = coef(covb)
    beta.0[S0] = beta0$par
    beta.0 = as.matrix(beta.0)
    theta.0<-as.vector(X%*%beta.0)
    r_yhatb0<-1-(1+exp(theta.0))^-1
    
    beta.1 = coef(covb)
    beta.1[S1] = beta1$par
    beta.1 = as.matrix(beta.1)
    theta.1<-as.vector(X%*%beta.1)
    r_yhatb1<-1-(1+exp(theta.1))^-1
    
    r_yhat1   <- predict(cov1,newx=x1hat,s='lambda.min')
    r_yhat0   <- predict(cov0,newx=x0hat,s='lambda.min')
    
    r_yhat1full = as.vector(predict(cov1,newx=x,s='lambda.min'))
    r_yhat0full = as.vector(predict(cov0,newx=x,s='lambda.min'))
    
  } else{
    
    X = cbind(rep(1,n),x)
    beta.0 = coef(covb)
    beta.0[c(1,S0)] = beta0$par
    beta.0 = as.matrix(beta.0)
    theta.0<-as.vector(X%*%beta.0)
    r_yhatb0<-1-(1+exp(theta.0))^-1
    
    beta.1 = coef(covb)
    beta.1[c(1,S1)] = beta1$par
    beta.1 = as.matrix(beta.1)
    theta.1<-as.vector(X%*%beta.1)
    r_yhatb1<-1-(1+exp(theta.1))^-1
    
    r_yhat1   <- predict(cov1,newx=x1hat,s='lambda.min',type = "response")
    r_yhat0   <- predict(cov0,newx=x0hat,s='lambda.min',type = "response")
    
    r_yhat1full = as.vector(predict(cov1,newx=x,s='lambda.min', type = "response"))
    r_yhat0full = as.vector(predict(cov0,newx=x,s='lambda.min', type = "response"))
  }
  
  delta_K <- sum((r_yhat1full-r_yhat0full-rep(ATE,n))^2)
  sigma_1 <- sum((r_yhat1-y1hat)^2/r_yhatb1[treat==1])/n
  sigma_0 <- sum((r_yhat0-y0hat)^2/(1-r_yhatb0[treat==0]))/n
  s = sqrt((delta_K+sum(sigma_1/r_yhatb1)+sum(sigma_0/r_yhatb0))/n)/sqrt(n)
  
  if (Att==1){
    ATT.delta_K <- sum(ATT.r_yhatb0*(r_yhat1full-r_yhat0full-rep(ATT,n))^2)
    w = (n/sum(treat))*sqrt((ATT.delta_K+sum(ATT.r_yhatb0[treat==1]*(r_yhat1-y1hat)^2)+sum(ATT.r_yhatb0[treat==0]^2*(r_yhat0-y0hat)^2/(1-ATT.r_yhatb0[treat==0])))/n)/sqrt(n)
  }
  
  
  fitted.values = rep(1,n)
  fitted.values[treat==1] = 1/(w.curr1$W[treat==1]+1)
  fitted.values[treat==0] = 1-1/(1-w.curr0$W[treat==0])
  
  output =list()
  output$ATE = ATE
  output$ATT = ATT
  output$s = s
  output$w = w
  output$test1 = w.curr1$W
  output$test0 = w.curr0$W
  output$coefficients1 = beta.1
  output$coefficients0 = beta.0
  output$fitted.values = fitted.values
  output$fitted.y = y
  output$fitted.x = x
  output
}
