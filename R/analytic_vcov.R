#' @title Calculate Variance-Covariance Matrix for Outcome Model
#' 
#' @description
#' \code{vcov_outcome} Returns the variance-covariance matrix of the main
#' parameters of a fitted CBPS object.
#' 
#' This adjusts the standard errors of the weighted regression of Y on Z for
#' uncertainty in the weights.
#' 
#' @aliases vcov_outcome vcov_outcome.CBPSContinuous
#' @param object A fitted CBPS object.
#' @param Y The outcome.
#' @param Z The covariates (including the treatment and an intercept term) that
#' predict the outcome.
#' @param delta The coefficients from regressing Y on Z, weighting by the
#' cbpsfit$weights.
#' @param tol Tolerance for choosing whether to improve conditioning of the "M"
#' matrix prior to conversion.  Equal to 1/(condition number), i.e. the
#' smallest eigenvalue divided by the largest.
#' @param lambda The amount to be added to the diagonal of M if the condition
#' of the matrix is worse than tol.
#' @return A matrix of the estimated covariances between the parameter
#' estimates in the weighted outcome regression, adjusted for uncertainty in
#' the weights.
#' @author Christian Fong, Chad Hazlett, and Kosuke Imai.
#' @references Lunceford and Davididian 2004.
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
#' @export vcov_outcome
#' 
vcov_outcome<-function(object, Y, Z, delta, tol=10^(-5), lambda=0.01)
{
  UseMethod("vcov_outcome")
}

#' vcov_outcome
#'
#' @export
#'
vcov_outcome.CBPSContinuous <- function(object, Y, Z, delta, tol=10^(-5), lambda=0.01){
  Xtilde <- object$Xtilde
  Ttilde <- object$Ttilde
  w <- object$weights
  beta.tilde <- object$beta.tilde
  sigmasq.tilde <- object$sigmasq.tilde
  
  N <- length(Y)
  K <- ncol(Xtilde)
  P <- ncol(Z)
   
  Sdelta <- matrix(0, nrow = P, ncol = P)
  Stheta <- matrix(0, nrow = P, ncol = K+1)
  
  # Precompute residuals since we use them a lot
  eps.beta <- as.vector(Ttilde - Xtilde%*%beta.tilde)
  eps.delta <- as.vector(Y - Z%*%delta)
  
  M11 <- apply(-2/sigmasq.tilde*eps.beta*Xtilde, 2, mean)
  M12 <- mean(-1/sigmasq.tilde^2*eps.beta^2)
  M22 <- apply(as.vector(1/(2*sigmasq.tilde)*w*(1 - 1/sigmasq.tilde*eps.beta^2)*Ttilde)*Xtilde, 2, mean)
  M21 <- matrix(0, nrow = K, ncol = K)

  for (i in 1:N){
    # Just added a -1 to Sdelta, I think it's correct.  Doesn't actually make a difference in V
    Sdelta <- Sdelta - w[i]*Z[i,]%*%t(Z[i,])/N
    M21 <- M21 + as.vector(-1/sigmasq.tilde*w[i]*Ttilde[i]*eps.beta[i])*Xtilde[i,]%*%t(Xtilde[i,])/N
    Stheta <- Stheta + cbind(-1/sigmasq.tilde*w[i]*eps.beta[i]*eps.delta[i]*Z[i,]%*%t(Xtilde[i,]),
                             1/(2*sigmasq.tilde)*w[i]*(1 - 1/sigmasq.tilde*eps.beta[i]^2)*eps.delta[i]*Z[i,])/N
  }
  M <- rbind(c(M11, M12), cbind(M21,M22))
  
  #Improve conditioning of M if necessary
  cond.num=svd(M)$d[1]/svd(M)$d[nrow(M)]
  if (cond.num>(1/tol)){M = M+lambda*diag(rep(1,nrow(M)))}
  
  s <- as.vector(w*eps.delta)*Z
  mtheta <- cbind(1/sigmasq.tilde*(eps.beta)^2 - 1, 
                  as.vector(w*Ttilde)*Xtilde)

  M.inv = solve(M)

  inner <- matrix(0, nrow = P, ncol = P)
  for (i in 1:N){
    inner.part <- s[i,] - Stheta%*%M.inv%*%mtheta[i,]
    inner <- inner + inner.part%*%t(inner.part)/N
  }
  
  Sdelta.inv = solve(Sdelta)
  
  V <- Sdelta.inv %*% inner %*% t(Sdelta.inv)/N
  
  return(V)
}


