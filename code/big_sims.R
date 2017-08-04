###############################
## Filename: big_sims.r
## Created by: David Gerard
## Synopsis: Try out all of the functions that I've been working with and compare risk
##          performances against max efron-morris and the candes estimators under various scenarios
## Compare mode specific soft-thresholding, truncation, candes, min efron-morris, and stein
## 1. Mode specific soft-thresholding works best when diagonal along each mode
## 2. Truncation works best when multilinear rank is low
## 3. Candes' works best when low rank but not low multilinear rank
## 4. Min efron-morris should work best when covariance along one mode, but doesn't really
## 5. Stein works best when Theta drawn randomly.
###############################


## Load required functions ---------------------------------------------------
library(hose)
library(tensr)
library(snow)
suppressMessages(library(parallel))
source("./code/misc_functions.R")

## Simulation settings used throughout R code --------------------------------
snr <- 1 ## this controls the signal to noise
p <- c(10, 10, 10) ## dim of array
n <- length(p)
sim_max <- 500 ## number of iterations per simulation
fnorm_mean <- snr * sqrt(prod(p))

#' One repetition of the simulations scenario.
#'
#' @param current_seed The current seed to start
#' @param theta The current true mean.
#'
#' @author David Gerard
one_rep <- function(current_seed, theta) {
  source("./code/misc_functions.R")
  p <- dim(theta)
  n <- length(p)

  func_lasso <- list()
  dfunc_lasso <- list()
  func_lasso[[1]] <- hose::f_lasso_mult
  dfunc_lasso[[1]] <- hose::df_lasso_mult
  for(mode_index in 2:n) {
    func_lasso[[mode_index]] <- hose::f_lasso
    dfunc_lasso[[mode_index]] <- hose::df_lasso
  }

  sure_vec <- rep(NA, 5)
  loss_vec <- rep(NA, 11)

  X <- theta + array(stats::rnorm(prod(p)), dim = p)

  ## mode-specific soft-thresholding ---------------------------------------
  lambda_init <- c(0, 0, 0)
  c_init <- 1
  c_x <- hose::get_c(X)
  soft_out <- hose::soft_coord(c_obj = c_x, lambda_init = lambda_init,
                               c_init = c_init, print_iter = FALSE, tau2 = 1)
  lambda_lasso <- lambda_list_wrapper(soft_out$lambda, soft_out$c)
  sure_final_lasso <- hose::sure(X, func = func_lasso, dfunc = dfunc_lasso, lambda_lasso)
  sure_vec[1] <- sure_final_lasso$sure_val
  loss_vec[1] <- sum((sure_final_lasso$mean_est - theta)^2)

  ## truncated HOSVD -------------------------------------------------------
  sure_rank_fit <- hose::sure_rank(X)
  rank_trunc <- c(as.matrix(sure_rank_fit$min_rank))
  final_est_trunc <- hose::sure(X, func = list(hose::f_truncate, hose::f_truncate, hose::f_truncate),
                                dfunc = list(hose::df_truncate, hose::df_truncate, hose::df_truncate),
                                lambda = sure_rank_fit$min_rank)
  sure_vec[2] <- final_est_trunc$sure_val
  loss_vec[2] <- sum((final_est_trunc$mean_est - theta)^2)

  ## matrix soft-thresholding ----------------------------------------------
  svd_mode1 <- svd(tensr::mat(X,1))
  candes_optim <- optim(1, fn = hose::sure_matrix, d = svd_mode1$d, p_dim = c(p[1], prod(p[-1])),
                        tau2 = 1, method = "Brent", lower = 0, upper = svd_mode1$d[1])
  candes_est <- array(svd_mode1$u %*% diag(pos_part(svd_mode1$d - candes_optim$par)) %*% t(svd_mode1$v), dim = p)
  sure_vec[3] <- candes_optim$value
  loss_vec[3] <- sum((candes_est - theta)^2)

  ## Efron-Morris ----------------------------------------------------------
  ef_est <- hose::efron_morris(X)
  sure_vec[4] <- ef_est$sure
  loss_vec[4] <- sum((ef_est$ef_est - theta)^2)

  ## Stein -----------------------------------------------------------------
  stein_est <- hose::stein(X)
  sure_vec[5] <- stein_est$sure_est
  loss_vec[5] <- sum((stein_est$est - theta)^2)

  ## Just X ----------------------------------------------------------------
  loss_vec[6] <- sum((X - theta)^2)

  ## HOOI using SURE estimate of rank --------------------------------------
  hx <- tensr::hooi(X = X, r = rank_trunc)
  hooi_est <- tensr::atrans(hx$G, hx$U)
  loss_vec[7] <- sum((hooi_est - theta) ^ 2)

  ## Cichocki's SCORE method -----------------------------------------------
  cichocki_est <- hose::score_ylc(Y = X, return_est = TRUE)
  loss_vec[8] <- sum((cichocki_est$est - theta) ^ 2)

  ## mode-specific mdl -----------------------------------------------------
  mdl_est <- hose::trunc_hosvd(Y = X, method = "mdl", return_est = TRUE)
  loss_vec[9] <- sum((mdl_est$est - theta) ^ 2)

  ## mode-specific SVA followed by truncation ------------------------------
  par_analysis_est <- hose::trunc_hosvd(Y = X, method = "pa", return_est = TRUE)
  loss_vec[10] <- sum((par_analysis_est$est - theta) ^ 2)

  ## bicrossvalidation followed by truncation ------------------------------
  bcv_est <- hose::trunc_hosvd(Y = X, method = "bcv", return_est = TRUE)
  loss_vec[11] <- sum((bcv_est$est - theta) ^ 2)

  ## Return values ---------------------------------------------------------
  return_vec <- c(loss_vec, sure_vec, rank_trunc,
                  cichocki_est$rank, mdl_est$rank,
                  par_analysis_est$rank, bcv_est$rank)
  names(return_vec) <- c("loss_soft", "loss_trunc", "loss_candes", "loss_em", "loss_stein",
                         "loss_x", "loss_hooi", "loss_cichocki", "loss_mdl", "loss_par",
                         "loss_bcv",
                         "sure_soft", "sure_trunc", "sure_candes", "sure_em",
                         "sure_stein", "sure_rank1", "sure_rank2", "sure_rank3",
                         "cichocki_rank1", "cichocki_rank2", "cichocki_rank3",
                         "mdl_rank1", "mdl_rank2", "mdl_rank3",
                         "par_rank1", "par_rank2", "par_rank3",
                         "bcv_rank1", "bcv_rank2", "bcv_rank3")
  return(return_vec)
}

#################
## Case where Truncation should work best ------------------------------------
## Half of sv's are 0, half are same
#################
set.seed(283128)
r <- c(5, 5, 5)
trash <- utils::capture.output(theta <- generate_core_same(p, r, fnorm_mean))

write.table(file = "./output/generated_thetas/theta_trunc.txt",
            tensr::mat(theta, 1), col.names = FALSE, row.names = FALSE)

seed_vec <- 1:sim_max
cl <- makeCluster(detectCores() - 1)
sout <- t(snow::parSapply(cl = cl, seed_vec, FUN = one_rep, theta = theta))
stopCluster(cl)

write.csv(sout, file = "./output/sims_out/losses_trunc.csv", row.names = FALSE)

###########################
## Now, simulations where evenly dispersed -----------------------------------
## Mode specific soft-thresholding should do best here
###########################
rm(theta)
set.seed(342)

D_list <- list()
for(mode_index in 1:length(p)) {
    D_list[[mode_index]] <- diag(p[mode_index]:1)
}

theta <- tensr::atrans(array(stats::rnorm(prod(p)), dim = p), D_list)
theta <- theta / tensr::fnorm(theta) * fnorm_mean

write.table(file = "./output/generated_thetas/theta_soft.txt",
            tensr::mat(theta, 1), col.names = FALSE, row.names = FALSE)

##S.temp <- hosvd.full(theta)$S ## sv's should be approximately evenly spaced
##sqrt(diag(mat(S.temp,1) %*% t(mat(S.temp,1))))
##sqrt(diag(mat(S.temp,2) %*% t(mat(S.temp,2))))
##sqrt(diag(mat(S.temp,3) %*% t(mat(S.temp,3))))

seed_vec <- 1:sim_max
cl <- makeCluster(detectCores() - 1)
sout <- t(snow::parSapply(cl = cl, seed_vec, FUN = one_rep, theta = theta))
stopCluster(cl)

write.csv(sout, file = "./output/sims_out/losses_soft.csv", row.names = FALSE)



#############################
## case where diag mult but only along one mode ------------------------------
## I think candes should work best here
#############################
rm(theta)
set.seed(1223)

D_mult <- diag(p[1]:1)
theta <- tensr::amprod(array(stats::rnorm(prod(p)), dim = p), D_mult,1)
theta <- theta / tensr::fnorm(theta) * fnorm_mean

write.table(file = "./output/generated_thetas/theta_candes.txt",
            tensr::mat(theta, 1), col.names = FALSE, row.names = FALSE)

seed_vec <- 1:sim_max
cl <- makeCluster(detectCores() - 1)
sout <- t(snow::parSapply(cl = cl, seed_vec, FUN = one_rep, theta = theta))
stopCluster(cl)

write.csv(sout, file = "./output/sims_out/losses_candes.csv", row.names = FALSE)



##################################
## Now case where theta is just drawn at random with no tensor structure -----
## Stein should do best here
##################################
rm(theta)
set.seed(563)
theta <- array(stats::rnorm(prod(p)), dim = p)
theta <- theta / tensr::fnorm(theta) * fnorm_mean

write.table(file = "./output/generated_thetas/theta_stein.txt",
            tensr::mat(theta, 1), col.names = FALSE, row.names = FALSE)

seed_vec <- 1:sim_max
cl <- makeCluster(detectCores() - 1)
sout <- t(snow::parSapply(cl = cl, seed_vec, FUN = one_rep, theta = theta))
stopCluster(cl)

write.csv(sout, file = "./output/sims_out/losses_stein.csv", row.names = FALSE)


##############################
## Now case where one mode has low rank, but not necessarily low multilinear rank
## I think Candes should do well here
##############################
rm(theta)
set.seed(8121)
theta_temp <- array(stats::rnorm(prod(p)), dim = c(p[1], prod(p[-1])))
svd_theta_temp <- svd(theta_temp)
svd_theta_temp$d[6:p[1]] <- 0
theta <- array(svd_theta_temp$u %*% diag(svd_theta_temp$d) %*% t(svd_theta_temp$v), dim = p)
theta <- theta / tensr::fnorm(theta) * fnorm_mean

write.table(file = "./output/generated_thetas/theta_candes2.txt",
            tensr::mat(theta, 1), col.names = FALSE, row.names = FALSE)

seed_vec <- 1:sim_max
cl <- makeCluster(detectCores() - 1)
sout <- t(snow::parSapply(cl = cl, seed_vec, FUN = one_rep, theta = theta))
stopCluster(cl)

write.csv(sout, file = "./output/sims_out/losses_candes2.csv", row.names = FALSE)




################################
## Now where there is covariance along one mode ------------------------------
## Min Efron-morris should do best here
################################
rm(theta)
set.seed(35879)
theta <- array(stats::rnorm(prod(p)), dim = p)
rho <- 0.7
cov_mode1 <- rho ^ abs(outer(1:p[1], 1:p[1], "-")) ## an AR1 correlation matrix
##cov.model <- rWishart(1,df = p[1],Sigma = diag(p[1]))[,,1]
theta <- tensr::amprod(theta, cov_mode1, 1)
theta <- theta / tensr::fnorm(theta) * fnorm_mean

write.table(file = "./output/generated_thetas/theta_em.txt",
            tensr::mat(theta, 1), col.names = FALSE, row.names = FALSE)

seed_vec <- 1:sim_max
cl <- makeCluster(detectCores() - 1)
sout <- t(snow::parSapply(cl = cl, seed_vec, FUN = one_rep, theta = theta))
stopCluster(cl)

write.csv(sout, file = "./output/sims_out/losses_em.csv", row.names = FALSE)
