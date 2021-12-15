#' Asymptotic Variance of the Average Treatment Effect Obtained with the oCBPS Procedure
#'
#' @param Y The vector of actual outcome values (observations).
#' @param Y_1_hat The vector of estimated outcomes according to the treatment model (user specified).
#' @param Y_0_hat The vector of estimated outcomes according to the control model (user specified).
#' @param CBPS_obj An object obtained with the CBPS function. This should include the vector of propensity scores and the estimated average treatment effect with the oCBPS method.
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
#'X_1 <- rnorm(n,3,sqrt(2))
#'X_2 <- rnorm(n,0,1)
#'X_3 <- rnorm(n,0,1)
#'X_4 <- rnorm(n,0,1)
#'X_mat <- as.matrix(cbind(rep(1,n), X_1, X_2, X_3, X_4)) 
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
#'#NOW ACTUALLY STARTING THE ESTIMATION
#'#First, fit the linear regression model for T_vec=1 and T_vec=0 separately.
#'Y_outcome_1 <- Y_outcome[which(T_vec==1)]
#'Y_outcome_0 <- Y_outcome[which(T_vec==0)]
#'           
#'#Separate the treatment covariates and the control covariates. 
#'X_mat_1 <- X_mat[which(T_vec==1),c(2,3,4,5)]
#'X_mat_0 <- X_mat[which(T_vec==0),c(2,3,4,5)]
#'           
#'#Perform linear regression separately.
#'lin_1 <- lm(Y_outcome_1 ~ X_mat_1) 
#'lin_0 <- lm(Y_outcome_0 ~ X_mat_0)
#'           
#'#Y_hat(1|X_i), #Y_hat(0|X_i) values 
#'Y_1_hat <- X_mat%*%as.matrix(lin_1$coefficients,5,1)  
#'Y_0_hat <- X_mat%*%as.matrix(lin_0$coefficients,5,1) 
#'           
#'#Use oCBPS to get estimated propensity scores pi and the estimated average treatment effect mu_hat.
#'ocbps.fit <- CBPS(T_vec ~ X_mat, ATT=0, baseline.formula = ~X_mat[,c(1,3:5)], 
#'                  diff.formula = ~X_mat[,2])
#'pi_hat <- ocbps.fit$fitted.values 
#'mu_hat <- mean(((T_vec)*Y_outcome/pi_hat) - (((1-T_vec)*Y_outcome)/(1-pi_hat)))
#'           
#'#Use the AsyVar function to get the asymptotic variance of the 
#'#estimated average treatment effect and its confidence interval
#'AsyVar(Y=Y_outcome, Y_1_hat=Y_1_hat, Y_0_hat=Y_0_hat, CBPS_obj=ocbps.fit, 
#'       X=X_mat, TL=T_vec, pi=pi_hat, mu=mu_hat, CI=0.95)
#'           
AsyVar <- function(Y, Y_1_hat, Y_0_hat, CBPS_obj=NULL, X, TL, pi, mu, CI=0.95){
  
  n <- length(Y)
  L_hat <- Y_1_hat - Y_0_hat
  
  if(is.null(CBPS_obj)==TRUE){ #Need X, TL, pi, mu)
    p <- ncol(X)
    n_1 <- length(which(TL==1))
    n_0 <- length(which(TL==0))
    sigma_hat_1_squared <- sum(((Y - Y_1_hat)*TL)^2)/(n_1-p)
    sigma_hat_0_squared <- sum(((Y - Y_0_hat)*(1-TL))^2)/(n_0-p)
    result <- list()
    result[[1]] <- mean((sigma_hat_1_squared)/pi + (sigma_hat_0_squared)/(1-pi) + (L_hat - mu)^2)
    
    result[[2]] <- result[[1]]/n
    
    diff_ocbps <- qnorm(1-(1-CI)/2)*sqrt(result[[1]]/n)
    lower_ocbps <- mu - diff_ocbps
    upper_ocbps <- mu + diff_ocbps  
    
    result[[3]] <- c(lower_ocbps,upper_ocbps)
    names(result) <- c("Asymptotic Variance of sqrt(n)*mu_hat", "Asymptotic Variance of mu_hat",  
                       paste(CI, "% Confidence Interval for mu_hat", sep=""))
    print(result)
  }
  else{
    p <- ncol(CBPS_obj$x)
    n_1 <- length(which(CBPS_obj$y==1))
    n_0 <- length(which(CBPS_obj$y==0))
    sigma_hat_1_squared <- sum(((Y - Y_1_hat)*CBPS_obj$y)^2)/(n_1-p)
    sigma_hat_0_squared <- sum(((Y - Y_0_hat)*(1-CBPS_obj$y))^2)/(n_0-p)
    mu_hat <- mean(((CBPS_obj$y)*Y/CBPS_obj$fitted.values) - 
                     (((1-CBPS_obj$y)*Y)/(1-CBPS_obj$fitted.values)))
    result <- list()
    result[[1]] <- mean((sigma_hat_1_squared)/CBPS_obj$fitted.values + 
                          (sigma_hat_0_squared)/(1-CBPS_obj$fitted.values) + (L_hat - mu_hat)^2)
    
    result[[2]] <- result[[1]]/n
    
    diff_ocbps <- stats::qnorm(1-(1-CI)/2)*sqrt(result[[1]]/n)
    lower_ocbps <- mu_hat - diff_ocbps
    upper_ocbps <- mu_hat + diff_ocbps  
    
    result[[3]] <- c(lower_ocbps,upper_ocbps)
    names(result) <- c("Asymptotic Variance of sqrt(n)*mu_hat", "Asymptotic Variance of mu_hat",  
                       paste(CI, "% Confidence Interval for mu_hat", sep=""))
    print(result)
  }
}



