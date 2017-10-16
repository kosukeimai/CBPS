library(numDeriv)
library(MASS)



#' Covariate Balancing Propensity Score (CBPS) for Marginal Structural Models
#' 
#' \code{CBMSM} estimates propensity scores such that both covariate balance
#' and prediction of treatment assignment are maximized.  With longitudinal
#' data, the method returns marginal structural model weights that can be
#' entered directly into a linear model.  The method also handles multiple
#' binary treatments administered concurrently.
#' 
#' Fits covariate balancing propensity scores for marginal structural models.
#' 
#' @aliases CBMSM CBMSM.fit
#' @param formula A list of formulas of the form treat ~ X.  The function
#' assumes that there is one formula for each time, and they are ordered from
#' the first time to the last time.
#' @param id A vector which identifies the unit associated with each row of
#' treat and X.
#' @param time A vector which identifies the time period associated with each
#' row of treat and X.
#' @param data An optional data frame, list or environment (or object coercible
#' by as.data.frame to a data frame) containing the variables in the model. If
#' not found in data, the variables are taken from \code{environment(formula)},
#' typically the environment from which \code{CBMSM} is called.
#' @param twostep Set to \code{TRUE} to use a two-step estimator, which will
#' run substantially faster than continuous-updating.  Default is \code{FALSE},
#' which uses the continuous-updating estimator described by Imai and Ratkovic
#' (2014).
#' @param msm.variance Default is \code{FALSE}, which uses the low-rank
#' approximation of the variance described in Imai and Ratkovic (2014).  Set to
#' \code{TRUE} to use the full variance matrix.
#' @param time.vary Default is \code{FALSE}, which uses the same coefficients
#' across time period.  Set to \code{TRUE} to fit one set per time period.
#' @param treat A vector of treatment assignments.  For N observations over T
#' time periods, the length of treat should be N*T.
#' @param X A covariate matrix.  For N observations over T time periods, X
#' should have N*T rows.
#' @param type "MSM" for a marginal structural model, with multiple time
#' periods or "MultiBin" for multiple binary treatments at the same time
#' period.
#' @param MultiBin.fit A parameter for whether the multiple binary treatments
#' occur concurrently (\code{FALSE}) or over consecutive time periods
#' (\code{TRUE}) as in a marginal structural model.  Setting type = "MultiBin"
#' when calling \code{CBMSM} will set MultiBin.fit to \code{TRUE} when
#' CBMSM.fit is called.
#' @param ... Other parameters to be passed through to \code{optim()}
#' @return \item{weights}{The optimal weights.} \item{fitted.values}{The fitted
#' propensity score for each observation.} \item{y}{The treatment vector used.}
#' \item{x}{The covariate matrix.} \item{id}{The vector id used in CBMSM.fit.}
#' \item{time}{The vector time used in CBMSM.fit.} \item{model}{The model
#' frame.} \item{call}{The matched call.} \item{formula}{The formula supplied.}
#' \item{data}{The data argument.} \item{treat.hist}{A matrix of the treatment
#' history, with each observation in rows and time in columns.}
#' \item{treat.cum}{A vector of the cumulative treatment history, by
#' individual.}
#' @author Marc Ratkovic, Christian Fong, and Kosuke Imai; The CBMSM function
#' is based on the code for version 2.15.0 of the glm function implemented in
#' the stats package, originally written by Simon Davies.  This documenation is
#' likewise modeled on the documentation for glm and borrows its language where
#' the arguments and values are the same.
#' @seealso \link{plot.CBMSM}
#' @references
#' 
#' Imai, Kosuke and Marc Ratkovic.  2014. ``Covariate Balancing Propensity
#' Score.'' Journal of the Royal Statistical Society, Series B (Statistical
#' Methodology). \url{http://imai.princeton.edu/research/CBPS.html}
#' 
#' Imai, Kosuke and Marc Ratkovic.  2015.  ``Robust Estimation of Inverse
#' Probability Weights for Marginal Structural Models.'' Journal of the
#' American Statistical Association.
#' \url{http://imai.princeton.edu/research/MSM.html}
#' @examples
#' 
#' 
#' ##Load Blackwell data
#' 
#' data(Blackwell)
#' 
#' ## Quickly fit a short model to test
#' form0 <- "d.gone.neg ~ d.gone.neg.l1 + camp.length"
#' fit0<-CBMSM(formula = form0, time=Blackwell$time,id=Blackwell$demName,
#' 			data=Blackwell, type="MSM",  iterations = NULL, twostep = TRUE, 
#' 			msm.variance = "approx", time.vary = FALSE)
#' 
#' \dontrun{
#' ##Fitting the models in Imai and Ratkovic  (2014)		
#' ##Warning: may take a few mintues; setting time.vary to FALSE
#' ##Results in a quicker fit but with poorer balance
#' ##Usually, it is best to use time.vary TRUE
#' form1<-"d.gone.neg ~ d.gone.neg.l1 + d.gone.neg.l2 + d.neg.frac.l3 + 
#' 		camp.length + camp.length + deminc + base.poll + year.2002 + 
#' 		year.2004 + year.2006 + base.und + office"
#' 		
#' fit1<-CBMSM(formula = form1, time=Blackwell$time,id=Blackwell$demName,
#' 			data=Blackwell, type="MSM",  iterations = NULL, twostep = TRUE, 
#' 			msm.variance = "full", time.vary = TRUE)
#' 
#' fit2<-CBMSM(formula = form1, time=Blackwell$time,id=Blackwell$demName,
#' 			data=Blackwell, type="MSM",  iterations = NULL, twostep = TRUE, 
#' 			msm.variance = "approx", time.vary = TRUE)
#' 
#' 
#' ##Assessing balance
#' 
#' bal1<-balance.CBMSM(fit1)
#' bal2<-balance.CBMSM(fit2)
#' 
#' ##Effect estimation: Replicating Effect Estimates in 
#' ##Table 3 of Imai and Ratkovic (2014)
#' 
#' lm1<-lm(demprcnt[time==1]~fit1$treat.hist,data=Blackwell,
#' weights=fit1$glm.weights)
#' lm2<-lm(demprcnt[time==1]~fit1$treat.hist,data=Blackwell,
#' weights=fit1$weights)
#' lm3<-lm(demprcnt[time==1]~fit1$treat.hist,data=Blackwell,
#' weights=fit2$weights)
#' 
#' lm4<-lm(demprcnt[time==1]~fit1$treat.cum,data=Blackwell,
#' weights=fit1$glm.weights)
#' lm5<-lm(demprcnt[time==1]~fit1$treat.cum,data=Blackwell,
#' weights=fit1$weights)
#' lm6<-lm(demprcnt[time==1]~fit1$treat.cum,data=Blackwell,
#' weights=fit2$weights)
#' 
#' 
#' 
#' ### Example: Multiple Binary Treatments Administered at the Same Time
#' n<-200
#' k<-4
#' set.seed(1040)
#' X1<-cbind(1,matrix(rnorm(n*k),ncol=k))
#' 
#' betas.1<-betas.2<-betas.3<-c(2,4,4,-4,3)/5
#' probs.1<-probs.2<-probs.3<-(1+exp(-X1 %*% betas.1))^-1
#' 
#' treat.1<-rbinom(n=length(probs.1),size=1,probs.1)
#' treat.2<-rbinom(n=length(probs.2),size=1,probs.2)
#' treat.3<-rbinom(n=length(probs.3),size=1,probs.3)
#' treat<-c(treat.1,treat.2,treat.3)
#' X<-rbind(X1,X1,X1)
#' time<-c(rep(1,nrow(X1)),rep(2,nrow(X1)),rep(3,nrow(X1)))
#' id<-c(rep(1:nrow(X1),3))
#' y<-cbind(treat.1,treat.2,treat.3) %*% c(2,2,2) + 
#' X1 %*% c(-2,8,7,6,2) + rnorm(n,sd=5)
#' 
#' multibin1<-CBMSM(treat~X,id=id,time=time,type="MultiBin",twostep=TRUE)
#' summary(lm(y~-1+treat.1+treat.2+treat.3+X1, weights=multibin1$w))
#' }
#' 
#' @export CBMSM
#' 
CBMSM<-function(formula, id, time, data, type="MSM", twostep = TRUE, msm.variance = "approx", time.vary = FALSE, ...){
  if (missing(data)) 
    data <- environment(formula)
  call <- match.call()
  family <- binomial()
  
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)

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
  
  
  ##Format treatment matrix
     	id<-as.numeric(as.factor(id))
		unique.id<-sort(unique(id))
		treat.hist<-matrix(NA,nrow=length(unique.id),ncol=length(unique(time)))
		colnames(treat.hist)<-sort(unique(time))
		rownames(treat.hist)<-unique.id
		for(i in 1:length(unique(unique.id))) for(j in sort(unique(time)))
		{{
		treat.hist[i,j]<-Y[id==unique.id[i] & time==j]
		}}
		#treat.hist.fac<-apply(treat.hist,1,function(x) paste(x, collapse="+"))
		cm.treat<-rowSums(treat.hist)

  
  if(type=="MSM")   {
    MultiBin.fit<-FALSE
  }
  
  if(type=="MultiBin"){
    MultiBin.fit<-TRUE
    X<-cbind(1,X[,apply(X,2,sd)>0])
    names.X<-c("Intercept",colnames(X)[-1])
  }

  fit <- eval(call("CBMSM.fit", treat = Y, X = X, id = id, time=time, 
                   MultiBin.fit = MultiBin.fit, twostep = twostep, msm.variance = msm.variance,
                   time.vary = time.vary))    
  
  fit$call<-call
  fit$formula<-formula
  fit$y<-Y
  fit$x<-X
  fit$id<-id
  fit$time<-time
  fit$model<-mf
  fit$data<-data
  fit$treat.hist<-treat.hist
  fit$treat.cum<-rowSums(treat.hist)
  fit$weights<-fit$weights[time==min(time)]
  fit
}


########################
###Calls loss function
########################
CBMSM.fit<-function(treat, X, id, time, MultiBin.fit, twostep, msm.variance, time.vary, ...){
id0<-id
id<-as.numeric(as.factor(id0))

if(msm.variance=="approx") full.var<-FALSE
if(msm.variance=="full") full.var<-TRUE

X.mat<-X
X.mat<-X.mat[,apply(X.mat,2,sd)>0, drop = FALSE]

##Format design matrix, run glm
glm1<-glm(treat~X.mat,family="binomial")


##################
##Make SVD matrix of covariates
##and matrix of treatment history
##################
#if(time.vary==FALSE){
X.svd<-X.mat
#X.svd<-apply(X.svd,2,FUN=function(x) (x-mean(x))/sd(x), drop=FALSE)
X.svd<-scale(X.svd) # Edit by Christian; this was causing an error
#X.svd[,c(1,2,7)]<-X.svd[,c(1,2,7)]*10
X.svd<-svd(X.svd)$u%*%diag(svd(X.svd)$d>0.0001)
X.svd<-X.svd[,apply(X.svd,2,sd)>0,drop=FALSE]
glm1<-glm(treat~X.svd,family="binomial")
if(time.vary==TRUE){
#} else{
X.svd<-NULL
for(i in sort(unique(time))){
	X.sub<-X.mat[time==i,,drop=FALSE]
	#X.sub<-apply(X.sub,2,FUN=function(x) (x-mean(x))/sd(x))
	X.sub <- scale(X.sub) # Edit by Christian; this was causing an error
	X.sub[is.na(X.sub)]<-0
	X.sub<-svd(X.sub)$u%*%diag(svd(X.sub)$d>0.0001)
	X.sub<-X.sub[,apply(X.sub,2,sd)>0,drop=FALSE]
	X.svd<-rbind(X.svd,X.sub)
	}
##Make matrix of time-varying glm starting vals
glm.coefs<-NULL
n.time<-length(unique(time))
for(i in 1:n.time){
	glm.coefs<-cbind(glm.coefs, summary(glm(treat~X.svd, subset=(time==i)))$coef[,1])
	}
	glm.coefs[is.na(glm.coefs)]<-0
	glm1$coefficients<-as.vector(glm.coefs)
}

##################
## Start optimization
##################
#Twostep  is true
msm.loss1<-function(x,...) msm.loss.func(betas=x, X=cbind(1,X.svd), treat=treat, time=time,...)$loss
glm.fit<-msm.loss.func(glm1$coef,X=cbind(1,X.svd),time=time,treat=treat,full.var=full.var,twostep=FALSE)
#print(head(glm.fit$V))[,1:10]
##Twostep is true; full variance option is passed

if(twostep==TRUE){
 Vcov.inv<-glm.fit
  
 msm.opt<-optim(glm1$coef,msm.loss1,full.var=full.var,Vcov.inv=Vcov.inv$V,bal.only=TRUE,twostep=TRUE,method="BFGS")
 
 msm.fit<-msm.loss.func(msm.opt$par,X=cbind(1,X.svd), treat=treat, time=time, full.var=full.var,Vcov.inv=Vcov.inv$V,bal.only=TRUE,twostep=TRUE)
}

##Twostep is false; full variance option is passed

if(twostep==FALSE) {
	msm.opt<-optim(glm1$coef,msm.loss1,full.var=full.var,bal.only=TRUE,twostep=FALSE,method="BFGS")

	msm.fit<-msm.loss.func(msm.opt$par,X=cbind(1,X.svd), treat=treat, time=time, full.var=full.var,Vcov.inv=Vcov.inv$V,bal.only=TRUE,twostep=FALSE)
	}
	
##################
## Calculate unconditional probs and treatment matrix
##################

n.obs<-length(unique(id))
n.time<-length(unique(time))

treat.hist<-matrix(NA,nrow=n.obs,ncol=n.time)
name.cands<-sort(unique(id))
for(i in 1:n.obs) for(j in 1:n.time) treat.hist[i,j]<-treat[id==name.cands[i] & time==j ]

treat.hist.unique<-unique(treat.hist,MAR=1)

treat.unique<-rep(NA,n.obs)

for(i in 1:n.obs) treat.unique[i]<- which(apply(treat.hist.unique,1,FUN=function(x) sum((x-treat.hist[i,])^2) )==0)

treat.unique<-as.factor(treat.unique)

uncond.probs.cand<-rep(0,n.obs)
for(i in 1:n.obs) {for(j in 1:n.obs) {
		check<-mean(treat.hist[j,]==treat.hist[i,])==1
	if(check) uncond.probs.cand[i]<-uncond.probs.cand[i]+1
	}
	}
	
	uncond.probs.cand<-uncond.probs.cand/n.obs
	
###########
##Produce Weights
###########

wts.out<-rep(uncond.probs.cand/msm.fit$pr,n.time)[time==1]
probs.out<-msm.fit$pr
uncond.probs<-uncond.probs.cand

loss.glm<-glm.fit$loss
loss.msm<-msm.fit$loss

if(loss.glm<loss.msm){
	warning("CBMSM fails to improve covariate balance relative to MLE.  \n GLM loss:    ", glm.fit$loss, "\n CBMSM loss:  ", msm.fit$loss, "\n")
}

# I know I'm putting probs.out in the weights and wts.out in the fitted values, but that
# is what Marc said to do
out<-list("weights"=probs.out,"fitted.values"=wts.out,"id"=id0[1:n.obs],"glm.g"=glm.fit$g.all,"msm.g"=msm.fit$g.all,"glm.weights"=(uncond.probs/glm.fit$pr)[time==1])
class(out)<-c("CBMSM","list")
return(out)
}


########################
###Loss function for MSM
########################
		


msm.loss.func<-function(betas,X=X,treat=treat,time=time,bal.only=F,time.sub=0,twostep=FALSE, Vcov.inv=NULL,full.var=FALSE){

	if((length(betas)==dim(X)[2]) ) betas<-rep(betas, dim(X)[2]/length(betas))

	time<-time-min(time)+1
	unique.time<-sort(unique(time))
	n.t<-length(unique.time)
	n<-dim(X)[1]/n.t
	treat.use<-betas.use<-NULL
	X.t<-NULL
	for(i in 1:n.t){
	betas.use<-cbind(betas.use,betas[1:dim(X)[2]+(i-1)*dim(X)[2]  ])
	treat.use<-cbind(treat.use,treat[time==unique.time[i]])
	X.t<-cbind(X.t,X[time==unique.time[i],])
	}
	betas<-betas.use
	betas[is.na(betas)]<-0
	treat<-treat.use
	thetas<-NULL
	for(i in 1:n.t)
	thetas<-cbind(thetas,X[time==i,]%*%betas[,i] )

	probs.trim<-.0001
	probs<-(1+exp(-thetas))^(-1)
	probs<-pmax(probs,probs.trim)
	probs<-pmin(probs,1-probs.trim)
	probs.obs<-treat*probs+(1-treat)*(1-probs)



	w.each<-treat/probs+(1-treat)/(1-probs)#+(treat-probs)^2/(probs*(1-probs))
	w.all<-apply(w.each,1,prod)#*probs.uncond

	bin.mat<-matrix(0,nrow=(2^n.t-1),ncol=n.t)
	for(i in 1:(2^n.t-1)) bin.mat[i,(n.t-length(integer.base.b(i))+1):n.t]<-
	integer.base.b(i)

	num.valid.outer<-constr.mat.outer<-NULL
	for(i.time in 1:n.t){
		num.valid<-rep(0,dim(treat)[1])
		constr.mat.prop<-constr.mat<-matrix(0,nrow=dim(treat)[1],ncol=dim(bin.mat)[1])
		for(i in 1:dim(bin.mat)[1]){
			is.valid<-sum(bin.mat[i,(i.time):dim(bin.mat)[2]])>0
			if(is.valid){
				#for(i.wt in i.time:n.t) w.all.now<-w.all.now*1/(1+3*probs[,i.wt]*(1-probs[,i.wt]))
				constr.mat[,i]<-(w.all*(-1)^(treat%*%bin.mat[i,]))
				num.valid<-num.valid+1
				}else{
				constr.mat[,i]<-0
				}

			}	
			num.valid.outer<-c(num.valid.outer,num.valid)
			constr.mat.outer<-rbind(constr.mat.outer,constr.mat)
		}

	if(twostep==FALSE){
	if(full.var==TRUE){
	var.big<-0
	X.t.big<-matrix(NA,nrow=n.t*dim(X.t)[1],ncol=dim(X.t)[2])
	for(i in 1:n.t){X.t.big[1:dim(X.t)[1]+(i-1)*dim(X.t)[1],]<-X.t}
	for(i in 1:dim(X.t.big)[1]){
		mat1<-(X.t.big[i,])%*%t(X.t.big[i,]) 
		mat2<-constr.mat.outer[i,]%*%t(constr.mat.outer[i,])
		var.big<-var.big+mat2 %x%mat1
	}
	}
	}

	X.wt<-X.prop<-g.wt<-g.prop<-NULL
	for(i in 1:n.t){
		g.prop<-c(g.prop, 1/n*t(X[time==i,])%*%(treat[,i]-probs[,i]))
		g.wt<-rbind(g.wt,1/n*t(X[time==i,])%*%cbind(constr.mat.outer[time==i,])*(i>time.sub))
		X.prop.curr<-matrix(0,ncol=n,nrow=dim(X)[2])
		X.wt.curr<-matrix(0,ncol=n,nrow=dim(X)[2])
		X.prop<-rbind(X.prop,1/n^.5*t((X[time==i,]*(probs.obs[,i]*(1-probs.obs[,i]))^.5)))
		if(bal.only){
				X.wt<-rbind(X.wt,1/n^.5*t(X[time==i,]*unique(num.valid.outer)[i]^.5)) } else{
		X.wt<-rbind(X.wt,1/n^.5*t(X[time==i,]*w.all^.5*unique(num.valid.outer)[i]^.5))
		}
	}

	mat.prop<-matrix(0,nrow=n, ncol=dim(X.wt)[2])
	mat.prop[,1]<-1
	g.prop.all<-0*g.wt
	g.prop.all[,1]<-g.prop
	#g.prop.all<-g.prop
	if(bal.only==T) g.prop.all<-0*g.prop.all

	g.all<-rbind(g.prop.all,g.wt)
	X.all<-rbind(X.prop*(1-bal.only),X.wt)

	if(twostep==TRUE){
		var.X.inv<-Vcov.inv
	} else{
		if(full.var==FALSE){
			var.X.inv<-ginv((X.all)%*%t(X.all))}else{
			var.X.inv<-ginv(var.big/n)
			}
		}
	length.zero<-dim(g.prop.all)[2]#length(g.prop)
	#var.X[(length.zero+1):(2*length.zero),1:length.zero]<-0
	#var.X[1:length.zero,(length.zero+1):(2*length.zero)]<-0


	#print(dim(g.all))
	#print(dim(var.X.inv))
	if(full.var==TRUE) g.all<-as.vector(g.wt)
	loss<-t(g.all)%*%var.X.inv%*%g.all
	out=list("loss"=(sum(diag(loss)))*n,"Var.inv"=var.X.inv,"probs"=w.all,"g.all"=g.all)
	
	#t(g.prop)%*%ginv(X.prop%*%t(X.prop))%*%g.prop +sum(diag(t(g.wt)%*%ginv(X.wt%*%t(X.wt))%*%g.wt ))
		
	
}#closes msm.loss.func



########################
###Makes binary representation
########################
integer.base.b <-
function(x, b=2){
        xi <- as.integer(x)
        if(any(is.na(xi) | ((x-xi)!=0)))
                print(list(ERROR="x not integer", x=x))
        N <- length(x)
        xMax <- max(x)
        ndigits <- (floor(logb(xMax, base=2))+1)
        Base.b <- array(NA, dim=c(N, ndigits))
        for(i in 1:ndigits){#i <- 1
                Base.b[, ndigits-i+1] <- (x %% b)
                x <- (x %/% b)
        }
        if(N ==1) Base.b[1, ] else Base.b
} 

#' @export
balance.CBMSM<-function(object, ...)
{
  treat.hist<-matrix(NA,nrow=length(unique(object$id)),ncol=length(unique(object$time)))
  ids<-sort(unique(object$id))
  times<-sort(unique(object$time))
  for(i in 1:length(ids)) {
    for(j in 1:length(times)){
      treat.hist[i,j]<-object$y[object$id== ids[i] & object$time==j]
    }
  }
  
  treat.hist.fac<-apply(treat.hist,1,function(x) paste(x, collapse="+"))
  bal<-matrix(NA,nrow=(ncol(object$x)-1),ncol=length(unique(treat.hist.fac))*2)
  baseline<-matrix(NA,nrow=(ncol(object$x)-1),ncol=length(unique(treat.hist.fac))*2)
  cnames<-array()
    
  for (i in 1:length(unique(treat.hist.fac)))
  {
    for (j in 2:ncol(object$x))
    {
      bal[j-1,i]<-sum((treat.hist.fac==unique(treat.hist.fac)[i])*object$x[which(object$time == times[1]),j]*object$w)/sum(object$w*(treat.hist.fac == unique(treat.hist.fac)[i]))
      #bal[j-1,i]<-sum((treat.hist.fac==unique(treat.hist.fac)[i])*object$x[,j]*object$w)/sum(object$w*(treat.hist.fac == unique(treat.hist.fac)[i]))
      # print(c(j,i,bal[j-1,i]))
      bal[j-1,i+length(unique(treat.hist.fac))]<-bal[j-1,i]/sd(object$w*object$x[which(object$time == times[1]),j])
      #bal[j-1,i+length(unique(treat.hist.fac))]<-bal[j-1,i]/sd(object$w*object$x[,j])
      baseline[j-1,i]<-sum((treat.hist.fac==unique(treat.hist.fac)[i])*object$x[which(object$time == times[1]),j]*object$glm.w)/sum(object$glm.w*(treat.hist.fac == unique(treat.hist.fac)[i]))
      baseline[j-1,i+length(unique(treat.hist.fac))]<-bal[j-1,i]/sd(object$glm.w*object$x[which(object$time == times[1]),j])
      #baseline[j-1,i]<-sum((treat.hist.fac==unique(treat.hist.fac)[i])*object$x[,j]*object$glm.w)/sum(object$glm.w*(treat.hist.fac == unique(treat.hist.fac)[i]))
      #baseline[j-1,i+length(unique(treat.hist.fac))]<-bal[j-1,i]/sd(object$glm.w*object$x[,j])
    }
    bal[is.na(bal)]<-0
    baseline[is.na(baseline)]<-0
    cnames[i]<-paste0(unique(treat.hist.fac)[i],".mean")
    cnames[i+length(unique(treat.hist.fac))]<-paste0(unique(treat.hist.fac)[i],".std.mean")
  }
  colnames(bal)<-cnames
  rnames<-colnames(object$x)[-1]
  rownames(bal)<-rnames
  colnames(baseline)<-cnames
  rownames(baseline)<-rnames
  statbal<-sum((bal-bal[,1])*(bal!=0)^2)
  statloh<-sum((baseline-baseline[,1])*(baseline!=0)^2)

  list("Balanced"=bal, "Unweighted"=baseline, "StatBal")
}



#' Plotting CBPS Estimation for Marginal Structural Models
#' 
#' Plots the absolute difference in standardized means before and after
#' weighting.
#' 
#' Covariate balance is improved if the plot's points are below the plotted
#' line of y=x.
#' 
#' @param x an object of class \dQuote{CBMSM}.
#' @param covars Indices of the covariates to be plotted (excluding the
#' intercept).  For example, if only the first two covariates from
#' \code{balance} are desired, set \code{covars} to 1:2.  The default is
#' \code{NULL}, which plots all covariates.
#' @param silent If set to \code{FALSE}, returns the absolute imbalance for
#' each treatment history pair before and after weighting.  This helps the user
#' to create his or her own customized plot. Default is \code{TRUE}, which
#' returns nothing.
#' @param boxplot If set to \code{TRUE}, returns a boxplot summarizing the
#' imbalance on the covariates instead of a point for each covariate.  Useful
#' if there are many covariates.
#' @param ... Additional arguments to be passed to plot.
#' @return The x-axis gives the imbalance for each covariate-treatment history
#' pair without any weighting, and the y-axis gives the imbalance for each
#' covariate-treatment history pair after CBMSM weighting.  Imbalance is
#' measured as the absolute difference in standardized means for the two
#' treatment histories.  Means are standardized by the standard deviation of
#' the covariate in the full sample.
#' @author Marc Ratkovic and Christian Fong
#' @seealso \link{CBMSM}, \link{plot}
#' 
#' @export
#' 
plot.CBMSM<-function(x, covars = NULL, silent = TRUE, boxplot = FALSE, ...)
{
  bal.out<-balance.CBMSM(x)
  bal<-bal.out$Balanced
  baseline<-bal.out$Unweighted
  no.treats<-ncol(bal)/2
  
  if (is.null(covars))
  {
    covars<-1:nrow(bal)
  }
  
  covarlist<-c()
  contrast<-c()
  bal.std.diff<-c()
  baseline.std.diff<-c()
  
  treat.hist.names<-sapply(colnames(bal)[1:no.treats],function(s) substr(s, 1, nchar(s)-5))  
  
  for (i in covars)
  {    
    for (j in 1:(no.treats-1))
    {
      for (k in (j+1):no.treats)
      {
        covarlist<-c(covarlist, rownames(bal)[i])        
        contrast<-c(contrast, paste(treat.hist.names[j],treat.hist.names[k],sep=":",collapse=""))
        bal.std.diff<-c(bal.std.diff,abs(bal[i,no.treats+j] - bal[i,no.treats+k]))
        baseline.std.diff<-c(baseline.std.diff,abs(baseline[i,no.treats+j] - baseline[i,no.treats+k]))
      }
    }
  }
    
  range.x<-range.y<-range(c(bal.std.diff,baseline.std.diff))
  if (!boxplot){
    plot(x=baseline.std.diff,y=bal.std.diff,asp="1",xlab="Unweighted Regression Imbalance",ylab="CBMSM Imbalance", 
         xlim=range.x, ylim = range.y, main = "Difference in Standardized Means", ...)
    abline(0,1)    
  }
  else{
    boxplot(baseline.std.diff, bal.std.diff, horizontal = TRUE, yaxt = 'n', xlab = "Difference in Standardized Means", ...)
    axis(side=2, at=c(1,2),c("CBMSM Weighted", "Unweighted"))
  }

  if(!silent) return(data.frame("Covariate" = covarlist, "Contrast"=contrast, "Unweighted"=baseline.std.diff, "Balanced"=bal.std.diff))
}
