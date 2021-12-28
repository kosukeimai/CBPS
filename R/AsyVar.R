#' Asymptotic Variance of the Average Treatment Effect Obtained with the oCBPS Procedure
#'
#' @param Y The vector of actual outcome values (observations).
#' @param Y_1_hat The vector of estimated outcomes according to the treatment model (user specified).
#' @param Y_0_hat The vector of estimated outcomes according to the control model (user specified).
#' @param CBPS_obj An object obtained with the CBPS function. This should include the vector of propensity scores and the estimated average treatment effect with the oCBPS method.
#' @param method The specific CBPS method to be considered. Default is "CBPS", but "oCBPS" can also be selected.
#' @param X The matrix of covariates with the rows corresponding to the observations and the columns corresponding to the variables. The left most column must be a column of 1's for the intercept.
#' @param TL The vector of treatment labels. More specifically, the label is 1 if it is in the treatment group and 0 if it is in the control group.
#' @param pi The vector of estimated propensity scores
#' @param mu The estimated average treatment effect obtained with the oCBPS procedure.
#' @param CI The specified confidence level (between 0 and 1) for calculating the confidence interval forr the average treatment effect with oCBPS.
#'
#' @import stats
#' @return The estimated asymptotic variance of the average treatment effect estimate obtained with the oCBPS procedure.
#' @return The confidence interval of the average treatment effect estimate obtained with the oCBPS procedure.
#' @export
#'
#' @examples #GENERATING THE DATA
#'n=300
#'           
#'#Initialize the X matrix.
#'X_v1 <- rnorm(n,3,sqrt(2))
#'X_v2 <- rnorm(n,0,1)
#'X_v3 <- rnorm(n,0,1)
#'X_v4 <- rnorm(n,0,1)
#'X_mat <- as.matrix(cbind(rep(1,n), X_v1, X_v2, X_v3, X_v4)) 
#'           
#'#Initialize the Y_1 and Y_0 vector using the treatment model and the control model.
#'Y_1 <- X_mat %*% matrix(c(200, 27.4, 13.7, 13.7, 13.7), 5, 1) + rnorm(n)
#'Y_0 <- X_mat %*% matrix(c(200, 0 , 13.7, 13.7, 13.7), 5, 1) + rnorm(n)
#'           
#'#True Propensity Score calculation.
#'pre_prop <- X_mat[,2:5] %*% matrix(c(0, 0.5, -0.25, -0.1), 4, 1)
#'propensity_true <- (exp(pre_prop))/(1+(exp(pre_prop)))
#'           
#'#Generate T_vec, the treatment vector, with the true propensity scores.
#'T_vec <- rbinom(n, size=1, prob=propensity_true)
#'           
#'#Now generate the actual outcome Y_outcome (accounting for treatment/control groups).
#'Y_outcome <- Y_1*T_vec + Y_0*(1-T_vec)
#'           
#'#Use oCBPS.
#'ocbps.fit <- CBPS(T_vec ~ X_mat, ATT=0, baseline.formula = ~X_mat[,c(1,3:5)], 
#'                  diff.formula = ~X_mat[,2])
#'           
#'#Use the AsyVar function to get the asymptotic variance of the 
#'#estimated average treatment effect and its confidence interval when using oCBPS.
#'AsyVar(Y=Y_outcome, CBPS_obj=ocbps.fit, method="oCBPS", CI=0.95)
#'
#'#Use CBPS.
#'cbps.fit <- CBPS(T_vec ~ X_mat, ATT=0)
#'                 
#' 
#'#Use the AsyVar function to get the asymptotic variance of the
#'#estimated average treatment effect and its confidence interval when using CBPS.
#'AsyVar(Y=Y_outcome, CBPS_obj=cbps.fit, method="CBPS", CI=0.95)
#'           
AsyVar <- function(Y, Y_1_hat=NULL, Y_0_hat=NULL, CBPS_obj, method="CBPS",
                   X=NULL, TL=NULL, pi=NULL, mu=NULL, CI=0.95){
  
  #Define X in all cases.
  if(is.null(X)==TRUE && is.null(CBPS_obj)==TRUE){
    stop('Need to specify either a matrix of covariates (X) 
           or an object obtained with the CBPS function.')
  } else {
    if(is.null(X)==TRUE && is.null(CBPS_obj)==FALSE){
      X <- CBPS_obj$x
    }
  } 
  
  #Define TL in all cases.
  if(is.null(TL)==TRUE && is.null(CBPS_obj)==TRUE){
    stop('Need to specify either a vector of treatment labels (TL) 
           or an object obtained with the CBPS function.')
  } else {
    if(is.null(TL)==TRUE && is.null(CBPS_obj)==FALSE){
      TL <- CBPS_obj$y
    }
  } 
  
  #Define pi in all cases.
  if(is.null(pi)==TRUE && is.null(CBPS_obj)==TRUE){
    stop('Need to specify either a vector of estimated propensity scores (pi)
           or an object obtained with the CBPS function.')
  } else {
    if(is.null(pi)==TRUE && is.null(CBPS_obj)==FALSE){
      pi <- CBPS_obj$fitted.values
    }
  }
  
  #Define mu in all cases.
  if(is.null(mu)==TRUE && is.null(CBPS_obj)==TRUE){
    stop('Need to specify either an estimate of the average treatment effect (mu)
           or an object obtained with the CBPS function.')
  } else {
    if(is.null(mu)==TRUE && is.null(CBPS_obj)==FALSE){
      mu <- mean(((TL)*Y/pi) - 
                   (((1-TL)*Y)/(1-pi)))
    }
  }
  
  
  p <- ncol(X)
  n <- length(Y)
  n_1 <- length(Y_1)
  n_0 <- length(Y_0)
  
  #Fit the linear regression model for TL=1 and TL=0 separately.
  if(is.null(Y_1_hat)==TRUE){
    Y_1 <- Y[which(TL==1)]
    
    #Separate the treatment covariates and the control covariates. 
    X_1 <- as.matrix(X[which(TL==1),-1])
    
    #Perform linear regression separately.
    lin_1 <- lm(Y_1 ~ X_1) 
    
    #Y_hat(1|X_i) values 
    Y_1_hat <- X%*%as.matrix(lin_1$coefficients,p,1)  
  }
  
  if(is.null(Y_0_hat)==TRUE){
    Y_0 <- Y[which(TL==0)]
    
    #Separate the treatment covariates and the control covariates. 
    X_0 <- as.matrix(X[which(TL==0),-1])
    
    #Perform linear regression separately.
    lin_0 <- lm(Y_0 ~ X_0) 
    
    #Y_hat(0|X_i) values 
    Y_0_hat <- X%*%as.matrix(lin_0$coefficients,p,1)  
  }
  
  L_hat <- Y_1_hat - Y_0_hat
  K_hat <- Y_0_hat
  
  if(method=="oCBPS"){
    
    #Calculate the Var(Y_1|X_i) and Var(Y_0|X_i). 
    sigma_hat_1_squared <- sum(((Y - Y_1_hat)*TL)^2)/(n_1-p)
    sigma_hat_0_squared <- sum(((Y - Y_0_hat)*(1-TL))^2)/(n_0-p)
    
    result <- list()
    result[[1]] <- mean((sigma_hat_1_squared)/pi + (sigma_hat_0_squared)/(1-pi) + (L_hat - mu)^2)
    
    result[[2]] <- result[[1]]/n
    
    diff_ocbps <- stats::qnorm(1-(1-CI)/2)*sqrt(result[[2]])
    lower_ocbps <- mu - diff_ocbps
    upper_ocbps <- mu + diff_ocbps  
    
    result[[3]] <- c(lower_ocbps,upper_ocbps)
    names(result) <- c("Asymptotic Variance of sqrt(n)*hat(mu)_{oCBPS}", "Asymptotic Variance of hat(mu)_{oCBPS}",  
                       paste(CI, "% Confidence Interval for hat(mu)_{oCBPS}", sep=""))
    print(result)
  } else {
    if(method=="CBPS"){
      
      #omega_hat (Look at the Var(g_beta) formula on page 7 in the paper! 
      #Because this is omega_hat according to page 35 right under (A.2).)
      new_X_2 <- array(rep(NA,p*p*n),dim=c(p,p,n))
      for(i in 1:n){
        new_X_2[,,i] <- matrix(X[i,],p,1)%*%matrix(X[i,],1,p)/(pi[i]*(1-pi[i]))
      }
      omega_hat <- apply(new_X_2,c(1,2),mean)
      
      
      #Sigma_mu_hat
      Sigma_mu_hat <- mean(((Y_1_hat)^2/pi) + ((Y_0_hat)^2/(1-pi))) - (mu)^2
      
      
      #Cov(mu_beta, g_beta)
      new_X_3 <- matrix(rep(NA,n*p),n,p)
      for(i in 1:n){
        new_X_3[i,] <- X[i,]*( (Y_1_hat[i]/pi[i]) + (Y_0_hat[i]/(1-pi[i])))
      }
      
      cov_hat <- apply(new_X_3,2,mean)
      
      
      #H_0_hat
      prop_modified <- pi/(1+exp(matrix(CBPS_obj$coefficients,1,p)%*%t(X)))
      der_mat <- matrix(rep(prop_modified,each=p),p,n)
      derivative <- matrix(rep(NA,p*n),p,n)
      for(i in 1:n){
        derivative[,i] <- matrix(der_mat[,i]*t(X[i,]),p,1) 
      }
      
      temp_sum <- matrix(rep(0,p),p,1)
      for(i in 1:n){
        temp_sum <- temp_sum + derivative[,i]*( (Y_1_hat[i]/pi[i]) + (Y_0_hat[i]/(1-pi[i])))
      }
      H_0_hat <- -temp_sum/n
      
      
      #H_f_hat
      temp_sum_2 <- matrix(rep(0,p*p),p,p)
      for(i in 1:n){
        temp_sum_2 <- temp_sum_2 + ( matrix(X[i,],p,1)%*%matrix(derivative[,i],1,p)/(pi[i]*(1-pi[i]))  )
      }
      H_f_hat <- -temp_sum_2/n
      
      
      #Put it all together
      result <- list()
      
      result[[1]] <- Sigma_mu_hat + t(H_0_hat)%*%solve(t(H_f_hat)%*%solve(omega_hat)%*%H_f_hat)%*%H_0_hat -
        2*t(H_0_hat)%*%solve(H_f_hat)%*%cov_hat
      
      result[[2]] <- result[[1]]/n
      
      diff_ocbps <- stats::qnorm(1-(1-CI)/2)*sqrt(result[[2]])
      lower_ocbps <- mu - diff_ocbps
      upper_ocbps <- mu + diff_ocbps  
      
      result[[3]] <- c(lower_ocbps, upper_ocbps)
      names(result) <- c("Asymptotic Variance of sqrt(n)*hat(mu)_{CBPS}", "Asymptotic Variance of hat(mu)_{CBPS}",  
                         paste(CI, "% Confidence Interval for hat(mu)_{CBPS}", sep=""))
    }
  }
}



