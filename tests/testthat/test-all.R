rm(list=ls())
library(CBPS)
library(testthat)
context("tests CBPS")

accuracy <- 0.0001
if (grepl('SunOS',Sys.info()['sysname'], ignore.case = TRUE)) {
  accuracy <- 0.005
}
if (grepl('Solaris',Sys.info()['sysname'], ignore.case = TRUE)) {
  accuracy <- 0.005
}
if (grepl('solaris',R.version$platform, ignore.case = TRUE)) {
  accuracy <- 0.005
}
if (grepl('solaris',R.version$os, ignore.case = TRUE)) {
  accuracy <- 0.005
}
accuracy <- ifelse(capabilities("long.double"), accuracy, 0.05)

test_that("tests CBMS on the Lalonde data", {
  # set random seed
  set.seed(12345)
  
  # Load the LaLonde data
  data(LaLonde)

  # Estimate CBPS via logistic regression
  fit <- CBPS(treat ~ age + educ + re75 + re74 + I(re75==0) + I(re74==0), data = LaLonde, ATT = TRUE)
  x <- vcov(fit)

  expect_equal(fit$coefficients[2], -0.06388547, tolerance = accuracy)
  expect_equal(fit$coefficients["educ",1], -0.08668682, tolerance = accuracy)
  expect_that(dim(x), is_equivalent_to(c(7,7)))
  
  if (R.Version()$arch != 'i386') {
    expect_equal(x["age", "age"], 0.01207999, tolerance = accuracy)
    expect_equal(x["re74", "re75"], -0.03111142, tolerance = accuracy)
  } else {
    expect_equal(x["age", "age"], 0.009398465, tolerance = accuracy)
    expect_equal(x["re74", "re75"], -0.03388868, tolerance = accuracy)
  }
})  

test_that("tests CBMS on the Blackwell data", {
  # set random seed
  set.seed(12345)

  #Load Blackwell data
  data(Blackwell)

  # Quickly fit a short model to test
  form0 <- "d.gone.neg ~ d.gone.neg.l1 + camp.length"
  x<-CBMSM(formula = form0, time=Blackwell$time,id=Blackwell$demName, data=Blackwell, type="MSM",  
           iterations = NULL, twostep = TRUE, msm.variance = "approx", time.vary = FALSE)
 
  expect_that(length(x), is_equivalent_to(15))
  expect_true("glm.weights" %in% names(x))
  expect_equal(as.numeric(x$fitted.values[100]), 0.01368907, tolerance = accuracy)
  expect_equal(x$glm.g[25,5], -1.23216, tolerance = accuracy)
  expect_equal(x$msm.g[25,5], -1.104169, tolerance = accuracy)

  donotrun <- TRUE
  if (! donotrun) {
    # Fitting the models in Imai and Ratkovic  (2014)		
    form1<-"d.gone.neg ~ d.gone.neg.l1 + d.gone.neg.l2 + d.neg.frac.l3 + camp.length + camp.length + 
            deminc + base.poll + year.2002 + year.2004 + year.2006 + base.und + office"
 		
    fit1<-CBMSM(formula = form1, time=Blackwell$time,id=Blackwell$demName,
  	  		data=Blackwell, type="MSM",  iterations = NULL, twostep = TRUE, 
  		  	msm.variance = "full", time.vary = TRUE)
  
    fit2<-CBMSM(formula = form1, time=Blackwell$time,id=Blackwell$demName,
  	  		data=Blackwell, type="MSM",  iterations = NULL, twostep = TRUE, 
  		  	msm.variance = "approx", time.vary = TRUE)

    #Assessing balance
    bal1<-balance.CBMSM(fit1)
    bal2<-balance.CBMSM(fit2)
 
    #Effect estimation: Replicating Effect Estimates in Table 3 of Imai and Ratkovic (2014)
    lm1<-lm(demprcnt[time==1]~fit1$treat.hist,data=Blackwell, weights=fit1$glm.weights)
    lm2<-lm(demprcnt[time==1]~fit1$treat.hist,data=Blackwell, weights=fit1$weights)
    lm3<-lm(demprcnt[time==1]~fit1$treat.hist,data=Blackwell, weights=fit2$weights)
    lm4<-lm(demprcnt[time==1]~fit1$treat.cum,data=Blackwell, weights=fit1$glm.weights)
    lm5<-lm(demprcnt[time==1]~fit1$treat.cum,data=Blackwell, weights=fit1$weights)
    lm6<-lm(demprcnt[time==1]~fit1$treat.cum,data=Blackwell, weights=fit2$weights)
 
    # Example: Multiple Binary Treatments Administered at the Same Time
    n<-200
    k<-4
    set.seed(1040)
    X1<-cbind(1,matrix(rnorm(n*k),ncol=k))
  
    betas.1<-betas.2<-betas.3<-c(2,4,4,-4,3)/5
    probs.1<-probs.2<-probs.3<-(1+exp(-X1 %*% betas.1))^-1
 
    treat.1<-rbinom(n=length(probs.1),size=1,probs.1)
    treat.2<-rbinom(n=length(probs.2),size=1,probs.2)
    treat.3<-rbinom(n=length(probs.3),size=1,probs.3)
    treat<-c(treat.1,treat.2,treat.3)
    X<-rbind(X1,X1,X1)
    time<-c(rep(1,nrow(X1)),rep(2,nrow(X1)),rep(3,nrow(X1)))
    id<-c(rep(1:nrow(X1),3))
    y<-cbind(treat.1,treat.2,treat.3) %*% c(2,2,2) + 
    X1 %*% c(-2,8,7,6,2) + rnorm(n,sd=5)
  
    multibin1<-CBMSM(treat~X,id=id,time=time,type="MultiBin",twostep=TRUE)
    x <- summary(lm(y~-1+treat.1+treat.2+treat.3+X1, weights=multibin1$w))
  }
})  
