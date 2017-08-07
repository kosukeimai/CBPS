CBPS.2Treat<-function(treat, X, method, k, XprimeX.inv, bal.only, iterations, ATT, standardize, twostep, sample.weights, ...){

	##Initialize data and some options
	#print("Test version")
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
		
		
	##Function for generating ATT weights. Called by loss function, etc.  
	##Note: Returns unsample weighted ATT weights
	ATT.wt.func<-function(beta.curr,X.wt=X,sample.weights0=sample.weights){
		X<-as.matrix(X.wt)
		n<-dim(X)[1]
		n.c<-sum(sample.weights0[treat==0])
		n.t<-sum(sample.weights0[treat==1])
		theta.curr<-as.vector(X%*%beta.curr)
		probs.curr<-(1+exp(-theta.curr))^-1
		probs.curr<-pmin(1-probs.min,probs.curr)
		probs.curr<-pmax(probs.min,probs.curr)	
		w1<-(n/n.t*(treat-probs.curr)/(1-probs.curr))
		w1[treat==1]<-n/n.t
		w1
	}
  

		
  ##The gmm objective function--given a guess of beta, constructs the GMM J statistic.  Used for binary treatments
	gmm.func<-function(beta.curr,X.gmm=X,ATT.gmm=ATT,invV=NULL,sample.weights0=sample.weights){

		sample.weights<-sample.weights0
		##Designate a few objects in the function.
		X<-as.matrix(X.gmm)
		ATT<-ATT.gmm
		
		##Designate sample size, number of treated and control observations,
		##theta.curr, which are used to generate probabilities.
		##Trim probabilities, and generate weights.
		##Note: n.c, n.t defined wrt sample weights
		n<-dim(X)[1]
		n.c<-sum(sample.weights[treat==0])
		n.t<-sum(sample.weights[treat==1])
		n<-n.c+n.t
		theta.curr<-as.vector(X%*%beta.curr)
		probs.curr<-(1+exp(-theta.curr))^-1
		probs.curr<-pmin(1-probs.min,probs.curr)
		probs.curr<-pmax(probs.min,probs.curr)	
		probs.curr<-as.vector(probs.curr)
		if(ATT)
			w.curr<-ATT.wt.func(beta.curr)
		else
			w.curr<-(probs.curr-1+treat)^-1
		  
	  
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
		V<-rbind(1/n*cbind(t(X.1)%*%X.1,t(X.1.1)%*%X.1.1)*n/n.t,
			     1/n*cbind(t(X.1.1)%*%X.1.1*n/n.t,t(X.2)%*%X.2*n^2/n.t^2))
		}
		else{
		V<-rbind(1/n*cbind(t(X.1)%*%X.1,t(X.1.1)%*%X.1.1),
			     1/n*cbind(t(X.1.1)%*%X.1.1,t(X.2)%*%X.2))
		}
		invV<-ginv(V)
		}			   
	
		##Calculate the GMM loss.
		loss1<-as.vector(t(gbar)%*%invV%*%(gbar))		
		out1<-list("loss"=max(loss1*n,loss1*n), "invV"=invV)
		out1
	}
	
	gmm.loss<-function(x,...) gmm.func(x,...)$loss

	XprimeX.inv<-ginv(t(sample.weights^.5*X)%*%(sample.weights^.5*X))

	##Loss function for balance constraints, returns the squared imbalance along each dimension.
	##Note zeroes out sample mean imbalance wrt sample weights
	bal.loss<-function(beta.curr,sample.weights0=sample.weights){
		sample.weights<-sample.weights0
		##Generate theta and probabilities.
		theta.curr<-as.vector(X%*%beta.curr)
		probs.curr<-(1+exp(-theta.curr))^-1
		probs.curr<-pmin(1-probs.min,probs.curr)
		probs.curr<-pmax(probs.min,probs.curr)
		##Generate weights.
		if(ATT)
			w.curr<-ATT.wt.func(beta.curr)
		else
			w.curr<-(probs.curr-1+treat)^-1
		X.2<-X
		##Generate mean imbalance.
		loss1<-abs(t(w.curr)%*%(sample.weights*X)%*%XprimeX.inv%*%t(sample.weights*X)%*%(w.curr))
		loss1
	}
  
  
  	#Redundant from above; delete
	#n<-sum(sample.weights)#length(treat)
	#n.c<-sum(sample.weights[treat==0])
	#n.t<-sum(sample.weights[treat==1])
	
	x.orig<-x<-cbind(as.matrix(X))

	##GLM estimation
	glm1<-suppressWarnings(glm(treat~X-1,family=binomial,weights=sample.weights))
	glm1$coef[is.na(glm1$coef)]<-0
	probs.glm<-glm1$fit
	glm1$fit<-probs.glm<-pmin(1-probs.min,probs.glm)
	glm1$fit<-probs.glm<-pmax(probs.min,probs.glm)	
	beta.curr<-glm1$coef
	beta.curr[is.na(beta.curr)]<-0

	##Do quick line search off beta.curr
	alpha.func<-function(alpha) gmm.loss(beta.curr*alpha)
	beta.curr<-beta.curr*optimize(alpha.func,interval=c(.8,1.1))$min
	
	##Generate estimates for balance and CBPSE
	gmm.init<-beta.curr
	this.invV<-gmm.func(gmm.init)$invV
  
	##Calculate imbalance minimizing estimates.
	opt.bal<-optim(gmm.init, bal.loss, control=list("maxit"=iterations), method="BFGS", hessian=TRUE)
	beta.bal<-opt.bal$par
    
    ##If only doing balance, stop, else GMM
	if(bal.only) {
		opt1<-opt.bal
		}else {
		
		##Initialize from glm and balance; choose minimum
		gmm.glm.init<-optim(gmm.init, gmm.loss, control=list("maxit"=iterations), method="BFGS", hessian=TRUE)
		gmm.bal.init<-optim(beta.bal, gmm.loss, control=list("maxit"=iterations), method="BFGS", hessian=TRUE)

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
  
	norm1<-norm2<-1
if (standardize)
	{
		if (ATT)
		{
			norm1<-sum(treat*n*sample.weights/sum(treat==1))
			norm2<-sum((1-treat)*n*sample.weights/sum(treat==1)*(treat-probs.opt)/(1-probs.opt))
		}
		else
		{
			norm1<-sum(treat/probs.opt*sample.weights)
			norm2<-sum((1-treat)/(1-probs.opt)*sample.weights)
		}
	
	
	if (ATT)
	{
		w.opt<-(treat == 1)*n/n.t/norm1*sample.weights + abs((treat == 0)*n/n.t*((treat - probs.opt)/(1-probs.opt))/norm2)*sample.weights
	}
	else
	{		
		w.opt<-(treat == 1)/probs.opt/norm1 + (treat == 0)/(1-probs.opt)/norm2
	}
	}
  
  	twostep<-F
	J.opt<-ifelse(twostep, gmm.func(beta.opt, invV = this.invV)$loss, gmm.loss(beta.opt))
  
	residuals<-treat-probs.opt
	deviance <- -2*c(sum(treat*log(probs.opt)+(1-treat)*log(1-probs.opt)))
	nulldeviance <- -2*c(sum(treat*log(mean(treat))+(1-treat)*log(1-mean(treat))))

	XG.1<- -X*(probs.opt)^.5*(1-probs.opt)^.5*sample.weights^.5
	XW.1<- X*(treat-probs.opt)*sample.weights^.5
	if(ATT==T){
	  XW.2<-X*(treat-probs.opt)/(1-probs.opt)*n/n.t*sample.weights^.5
	  XG.2<-X*((1-treat)*probs.opt/(1-probs.opt)*n/n.t)^.5*sample.weights^.5
	} 
	else{
	  XW.2<- X*(probs.opt-1+treat)^-1*sample.weights^.5
	  XG.2<- -X*probs.opt^.5*(1-probs.opt)^.5*abs((probs.opt-1+treat)^-1)*sample.weights^.5#*(abs(probs.opt-treat)/(probs.opt*(1-probs.opt)))^.5
  }
	if (twostep){
	  W<-this.invV
	}
	else{
	    W<-gmm.func(beta.opt)$invV
	}
	W1<-rbind(t(XW.1),t(XW.2))
	Omega<-(W1%*%t(W1)/n)
	G<-cbind(t(XG.1)%*%XG.1,t(XG.2)%*%XG.2)/n
	vcov<-ginv(G%*%W%*%t(G))%*%G%*%W%*%Omega%*%W%*%t(G)%*%ginv(G%*%W%*%t(G))


	output<-list("coefficients"=matrix(beta.opt, ncol=1),"fitted.values"=probs.opt,"deviance"=deviance,"weights"=w.opt*sample.weights,
				 "y"=treat,"x"=X,"converged"=opt1$conv,"J"=J.opt,"var"=vcov, 
				 "mle.J"=ifelse(twostep, gmm.func(glm1$coef, invV = this.invV)$loss, gmm.loss(glm1$coef)))

	class(output)<- c("CBPS","glm","lm")
    
	output
}


vcov_outcome.CBPS<-function(object,Y,Z=NULL,delta=NULL,tol=10^-5,lambda=0.01){
	
#	object: A fitted CBPS object
#Y: The outcome
#Z: The outcome covariates
#delta: Coefficients from regressing Y on Z
#tol = 10^-5: The tolerance for choosing whether to condition M prior to inversion.
#lambda = 0.01: The amount to be added to the diagonal of M if the 
#condition matrix is worse than tol.  More on this below.
	##Align objects with rest of package	
	y<-Y
	obj<-object
	
	method<-cb1$method
	##Gather model components form obj
	treat<-obj$y
	if(is.null(Z)) {
		X<-obj$x
		}else{
		X<-Z
		}
	sds<-apply(X,2,sd)
	X[,-1]<-t(apply(X[,-1],1,FUN=function(x) x/sds[-1]))
	Xt<-cbind(1,treat,X[,-1])
	sds2<-c(1,1,sds[-1])
	n<-length(treat)
	n1<-sum(treat)
	probs<-	 obj$fit

	##Test for ATT or ate
	ATT<-sd(obj$weights[treat==1])/sd(obj$weights[treat==0])<1e-6

	##Generate ATT
	if(ATT){
	wts<-(treat-probs)/(1-probs)*(n/n1)	
	wts.deriv<- -X*(n/n1)*(1-treat)*probs/(1-probs)
	} else{
	wts<-treat/probs-(1-treat)/(1-probs)
	wts.deriv<- -X*(treat*(1-probs)/probs +(1-treat)*probs/(1-probs))
	}
	wts.deriv.abs<-wts.deriv*sign(wts)
	wts.sq<-wts^2

	##Point estimates
	ests<-lm(y~Xt-1,w=abs(wts))
	errs<-ests$residuals

	##Don't need
	#wts.sq<-wts^2
	#W11<-t(Xt*errs^2*wts.sq)%*%Xt
	#W22<-t(X*wts.sq)%*%X
	#W1<-cbind(solve(W11),matrix(0,nrow=nrow(W11),ncol=ncol(W22)))
	#W2<-cbind(matrix(0,ncol=nrow(W11),nrow=nrow(W22)),solve(W22))
	#W<-rbind(W1,W2)

	##Calculate Jacobian
	if(method=="exact"){
	Gtilde1<- cbind(t(-Xt*abs(wts))%*%(Xt), t(Xt*errs)%*%wts.deriv.abs)
	Gtilde2<- cbind(matrix(0,nrow=ncol(X),ncol=ncol(Xt)),t(X)%*% (wts.deriv) ) 
	G<-rbind(Gtilde1,Gtilde2)
		}else{
	deriv.score<-X*(treat-probs)
	deriv.score.abs<-X*(2)
	Gtilde1<- cbind(t(-Xt*abs(wts))%*%(Xt), t(Xt*errs)%*%wts.deriv.abs)
	Gtilde2<- cbind(matrix(0,nrow=ncol(X),ncol=ncol(Xt)),t(X)%*% (wts.deriv) ) 
	G<-rbind(Gtilde1,Gtilde2)
	}

	##Calculate sample moment conditions
	if(method=="exact"){
	gtilde<-cbind(Xt*errs*abs(wts),X*(wts))
	}else{
	gtilde<-cbind(Xt*errs*abs(wts),X*(wts))
	}

	M<-t(gtilde)%*%gtilde
	   cond.num=svd(M)$d[1]/svd(M)$d[nrow(M)]
	   if (cond.num>(1/tol)){M = M+lambda*diag(rep(1,nrow(M)))}

	#Variance estimates
	V<-solve(t(G)%*%solve(M)%*%G)
	V2<-V[1:ncol(Xt),1:ncol(Xt)]/(sds2%*%t(sds2))
	
#	sd.out<-diag(V)^.5

#	return(list("coef"=ests$coef/sds2,"se"=sd.out[1:(ncol(Xt))]/sds2))
	V2
}
