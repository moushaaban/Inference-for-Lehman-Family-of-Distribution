# Inference-for-Lehman-Family-of-Distribution
Inference for Lehman Family of Distribution
   remove(list=objects())
  options(warn = -1)
 library(bbmle)
  n <-7
  r <- 1

  theta  <-1.8
  alpha <-0.6
 N.sim <- 100
 Array <- array(dim = c(n,n),NA)
 Mat <- matrix(NA, nrow = n*r, ncol =n*r) 
 
 Bais_theta   <- matrix(NA, nrow = N.sim, ncol = 1)
 Bais_alpha   <- matrix(NA, nrow = N.sim, ncol = 1)
 theta_hat    <- matrix(NA, nrow = N.sim, ncol = 1)
 alpha_hat    <- matrix(NA, nrow = N.sim, ncol = 1)
  for(k in 1:N.sim){
  for (j in 1:r){
   for (i in 1:n){
     i <- 1:n
     j <- 1:r
     u <- runif(n^2)
 
     
     Data=function(u, theta, alpha){ (sqrt((-1*log(1-(u^alpha)))/theta))}
 
 Array[i,] <- sort(Data(u, theta, alpha))
 }
 Mat[j,] <-  diag(Array) 
 }
 x <- sort( unique(na.omit(c(unique(Mat)))))
 #set.seed(k)
 
 fn <- function(theta, alpha){-(sum(log(choose(n, i)))+r*n*log(theta)
                                + r*n*log(alpha)- theta*sum(x^2)+n*log(2)
                                +sum(log(x))
                                +sum((alpha-1)*log(1-(exp(-theta*(x^2)))))
                                +sum(alpha*(i-1)*log(1-(exp(-theta*(x^2)))))
                                +sum((n-i)*log(1-((1-(exp(-theta*(x^2))))^alpha)))
 )}
 
 
#fn <- function(theta, alpha){-(sum(log(choose(n, i)))+r*n*log(theta)
             #     + r*n*log(alpha)- theta*sum(x)+sum((alpha-1)*(1-(exp(-theta*x))))
              #    +sum(alpha*(i-1)*log(1-(exp(-theta*x))))
               #   +sum((n-i)*log(1-((1-(exp(-theta*x)))^alpha))))}
 
 
 fit <- mle2(minuslogl = fn , start = list(theta = theta, alpha=alpha), data = list(x) , 
             method = "CG",lower = 0.1, upper = 1.5)
  
 # method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Brent")
 #x = x[x>=1]
 
    theta_hat[k,]   <- coef(fit)[1]
    alpha_hat[k,]   <- coef(fit)[2]  
    
    # Bais_theta[k,]   <- (coef(fit)[1]-theta)/theta
    Bais_theta[k,] <- theta - coef(fit)[1]
    
     #Bais_alpha[k,]   <- (coef(fit)[2]-alpha)/alpha
     Bais_alpha[k,] <- alpha - coef(fit)[2]
 }
    
   Baisoftheta   <- mean(Bais_theta[,1])
   Baisofalpha   <- mean(Bais_alpha[,1])
   
   MSEoftheta   <- mean((theta - theta_hat[k,])^2)
   MSEofalpha   <- mean((alpha - alpha_hat[k,])^2)
   
  
   MLEoftheta    <- mean(theta_hat[,1])
   MLEofalpha    <- mean(alpha_hat[,1])
   
   
   Summary <- data.frame(
                         "N.sim"=N.sim ,
                         "MLEofalpha"=MLEofalpha,
                         "MLEoftheta"=MLEoftheta,
                         "MSEofalpha"=MSEofalpha,
                         "MSEoftheta"=MSEoftheta,
                         "Baisofalpha"=Baisofalpha,
                          "Baisoftheta"=Baisoftheta
                       )
     print(Summary)
