## generate-data function
sim_fried <- function(n_coun, # vector: the numbers of individuals at each cluster
                      P, # integer: number of covariates
                      K, # integer: number of cluster
                      Sigma, # covariance matrix of covariates
                      sigma, # variance of model noise
                      lambda, # time weights in outcome process
                      #gamma = c(0.3, 0.5, 0.7), # vector: the parameters of hidden variables u, mild
                      #gamma = c(0.7, 1.2, 2.0), # Moderate
                      #gamma = c(2.0, 3.0, 5.0), # vector: the parameters of hidden variables u,
                      gamma = c(3.0, 9.0, 16.0), # vector: the parameters of hidden variables u,
                      rho = c(3.5, 2.0, 10.0)
) {
  
  N <- sum(n_coun) # total number
  Z <- matrix(0, ncol = K, nrow = N) #### the dummy variables in linear function.
  
  #X <- mvrnorm(N, rep(0, P), Sigma)
  
  rawvars <- mvrnorm(N, rep(0, P), Sigma)
  X <- pnorm(rawvars)
  
  Z[1:n_coun[1], 1] <- 1
  
  for (i in 2:K) {
    Z[(sum(n_coun[1:(i-1)]) + 1) : sum(n_coun[1:i]), i] <-  1
  }
  
  ## cluster-level process
  Sigma_beta <- matrix(c(
    gamma[1],       rho[1],     rho[2],
    rho[1],     gamma[2],       rho[3],
    rho[2],   rho[3],     gamma[3]
  ), nrow = 3, byrow = TRUE)
  
  Beta <- MASS::mvrnorm(
    n = K,
    mu = c(0, 0, 0),
    Sigma = Sigma_beta
  )
  

  eta <- matrix(NA, nrow =nrow(Z), ncol = ncol(Z)  )
  
    
    for (j in 1:3) {
      eta[,j] <- as.vector(Z %*% Beta[, j])
    }
  
  
  
  
  r_1 <- lambda[1] * (10 * sin(pi * X[,1] * X[,2])  + 20 * (X[,3] - 0.5)^2 + 10 * X[,4] + 5 * X[,5])
  #mu_1 <- r_1 + eta ##### time point = 0
  
  r_2 <- lambda[2] * (10 * sin(pi * X[,6] * X[,7]) + 20 * (X[,8] - 0.5)^2 + 10 * X[,9] + 5 * X[,10])
  #mu_2 <- r_2 + eta ###### time point = 1
  
  r_3 <- lambda[3] * (10 * sin(pi * X[,11] * X[,12]) + 20 * (X[,13] - 0.5)^2 + 10 * X[,14] + 5 * X[,15])
  #mu_3 <- r_3 + eta ###### time point = 2
  
  Y1 <- r_1  + eta[,1] + sigma * rnorm(N)
  Y2 <- r_2 + eta[,2]  + sigma * rnorm(N) #r_1   +  r_2 + eta[,2]  + sigma * rnorm(N)
  Y3 <- r_3 + eta[,3]  + sigma * rnorm(N) # r_1  +  r_2 +  r_3 + eta[,3]  + sigma * rnorm(N)
  
  return( list(df = data.frame(X[,1:5], Y1, X[,c(1:3,9:10)], Y2, X[,c(1:3,14:15)], Y3, Z = Z),
               Beta = Beta) ) # , mu = mu ,  A = A,
  
}

