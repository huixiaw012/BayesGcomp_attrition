###### This simulation is applied to evaluate of our proposed method when the dataset is under country-level attrition setting
#### Specifically, set X13-X15 and hY3 also have 5 clusters attrition randomly.
#### Here set beta_j and beta_j+1 are non-linear relationship

num_chains <- c(1:1000)

k <- 10
n_attr <- 5

ss <- 200 ### sample size
country_idx <- rep(1:k, each = ss)

lambda <- c(0.3,0.6,1)
#lambda <- c(0.5,1)
a1 <- 0.4
a2 <- 0.2
a3 <- 0.1
sig <- 1

Sig<-matrix(c(c(1,0,0,0,0,a1,0,0,0,0,a2,0,0,0,0),
              c(0,1,0,0,0,0,a1,0,0,0,0,a2,0,0,0),
              c(0,0,1,0,0,0,0,a1,0,0,0,0,a2,0,0),
              c(0,0,0,1,0,0,0,0,a1,0,0,0,0,a2,0),
              c(0,0,0,0,1,0,0,0,0,a1,0,0,0,0,a2),
              c(a1,0,0,0,0,1,0,0,0,0,a1,0,0,0,0),
              c(0,a1,0,0,0,0,1,0,0,0,0,a1,0,0,0),
              c(0,0,a1,0,0,0,0,1,0,0,0,0,a1,0,0),
              c(0,0,0,a1,0,0,0,0,1,0,0,0,0,a1,0),
              c(0,0,0,0,a1,0,0,0,0,1,0,0,0,0,a1),
              c(a2,0,0,0,0,a1,0,0,0,0,1,0,0,0,0),
              c(0,a2,0,0,0,0,a1,0,0,0,0,1,0,0,0),
              c(0,0,a2,0,0,0,0,a1,0,0,0,0,1,0,0),
              c(0,0,0,a2,0,0,0,0,a1,0,0,0,0,1,0),
              c(0,0,0,0,a2,0,0,0,0,a1,0,0,0,0,1)), 15, 15, byrow = TRUE)


library(glBART)
library(foreach)
library(doParallel)
library(dplyr)
library(MASS)
library(locfit)
library(rBeta2009)


########## fitted model

start_time <- Sys.time()
print(paste0("Start time", start_time))

n_burn <- 100 # increase for good mixing e.g. to 3000
n_thin <- 1 # increase for good mixing e.g. to 250
n_save <- 1 #
n_tree <- 20 # may need to increase

opts <- Opts(num_burn = n_burn, num_save = n_save)

cl <- parallel::makeCluster(16)
doParallel::registerDoParallel(cl)

# 1) load package + source once per worker
invisible(parallel::clusterEvalQ(cl, { 
  library(glBART)
  library(rBeta2009)
  library(dplyr)
  library(MASS)
  library(locfit)
  library(rBeta2009)
  source("sim_data_nonlin.R")
}) )


results <- foreach(i = num_chains, .combine = rbind) %dopar% {
  
  set.seed(i)
  
  
  dat <- sim_fried(n_coun = rep(ss,k),
                   P = 15,
                   K = k,
                   Sigma = Sig[1:15, 1:15],
                   sigma = sig,
                   lambda
                   )

  ave.true_y <- colMeans(dat$df[,c("Y1", "Y2", "Y3")])

  dat_with_NA <- dat$df ### original data

  # Compute missing probability for each group 
  missing_prob <- expit(-1.0 + 3.5 * (-dat$Beta[,2]) ) # severe

  miss_countries <- sample(1:k, size = n_attr, replace = FALSE, prob = missing_prob) ### Randomly

  dat_with_NA[country_idx %in% miss_countries, 13:18] <- NA ### Let the values of X_15 and Y_3 at T3 from missing countries become NAs.
  lin.formula <- as.formula(paste("~", paste0("Z.",  1:k, collapse  = "+") ))
    
    #dat_with_exclude <- dat_with_NA[!(country_idx %in% miss_countries), ]
    #dat_with_exclude <- dat_with_exclude[,!(names(dat_with_exclude) %in% paste0("Z.", miss_countries))] ### Let the values of X_15 and Y_3 at T3 from missing countries become NAs.
    #lin.formula <- as.formula(paste("~", paste0("Z.",  (1:k)[-miss_countries], collapse  = "+") ))
   
   mod_fit<- riBMfits(data = dat_with_NA, # dat_with_exclude, # A matrix or data frame possibly containing confounder(s), exposure(s), outcome(s) and mortality indicator(s) in temporal order from left to right.
                      var.type = c(rep("X0", 3), rep("X",2), "Y", rep("X", 5), "Y", rep("X", 5), "Y"), # X0=baseline confounders (use Bayesian bootstrap), X=confounders (for example time-varying), FI=exposure/treatment, Y=outcome
                      linear_formula = lin.formula,
                      fixed.regime = NULL, # Estimate the ACE for treated (Fi=1) and controls (Fi=0).
                      J = 2500, # Size of pseudo data. Increase to e.g. 5000
                      opts = opts, # opts see SoftBart
                      Suppress = TRUE, # Indicates if the output should be suppressed. Default is TRUE
                      By = 1000, # If Suppress is set to FALSE, output is provided for the By:th iteration.
                      weighted = FALSE,
                      num_tree = n_tree,
                      remove_intercept = TRUE,
                      print_iter = FALSE)



   pred_mod <- gcomp_risoftbart(data = dat_with_NA, # dat_with_exclude, # A matrix or data frame possibly containing confounder(s), exposure(s), outcome(s) and mortality indicator(s) in temporal order from left to right.
                                var.type = c(rep("X0", 3), rep("X",2), "Y", rep("X", 5), "Y", rep("X", 5), "Y"), # X0=baseline confounders (use Bayesian bootstrap), X=confounders (for example time-varying), FI=exposure/treatment, Y=outcome
                                linear_formula = lin.formula,
                                fixed.regime = NULL, # Estimate the ACE for treated (Fi=1) and controls (Fi=0).
                                J = 2500, # Size of pseudo data. Increase to e.g. 5000
                                opts = opts, # opts see SoftBart
                                Suppress = TRUE, # Indicates if the output should be suppressed. Default is TRUE
                                By = 1000, # If Suppress is set to FALSE, output is provided for the By:th iteration.
                                weighted = FALSE,
                                num_tree = n_tree,
                                BModels = mod_fit,
                                time.type = c(rep("T1", 6), rep("T2", 6), rep("T3",6)),
                                impute_method = "MVN", # BART or MVN
                                ridge = 1e-6,
                                remove_intercept = TRUE)


   
   Y_hat <- pred_mod$y_hat[,1]
   
   out <- c(true = ave.true_y, est = Y_hat)
   
  
}


names(results) <- c(paste0("true_Y", 1:3), paste0("est_Y", 1:3))

#print(results)

# Stop the parallel backend
stopImplicitCluster()

end_time <- Sys.time()
print(paste0("End time", end_time))



