library(hose) ## loads west, the dataset we will analyze.
library(tensr)

west_array     <- west$west_array ## arc-sin transformed proportions of observed games
west_array_n   <- west$west_array_n ## attempts
west_array_sub <- west$west_array_sub ## made
west_true_as   <- west$west_true_as ## arc-sin transformed proportions of unseen games
west_true_p    <- west$west_true_p ## proportions of unseen games

team_names <- dimnames(west_array)[[1]]
p <- dim(west_array)
n <- length(team_names)

arc_sin <- function(p, n) {
    sqrt(n) * asin(2 * p - 1)
}

inv_arc_sin <- function(x, n) {
    (sin(x / sqrt(n)) + 1)/2
}

func_trunc <- list()
dfunc_trunc <- list()
for(mode_index in 1:length(dim(west_array))) {
    func_trunc[[mode_index]] <- f_truncate
    dfunc_trunc[[mode_index]] <- df_truncate
}

func_lasso <- list()
dfunc_lasso <- list()
func_lasso[[1]] <- f_lasso_mult
dfunc_lasso[[1]] <- df_lasso_mult
for(mode_index in 2:length(dim(west_array))) {
    func_lasso[[mode_index]] <- f_lasso
    dfunc_lasso[[mode_index]] <- df_lasso
}

## Interpolate diagonals -----------------------------------------------------
for(mode3_index in 1:p[3]) {
    for(mode4_index in 1:p[4]) {
        diag(west_array[ , , mode3_index, mode4_index]) <-
          (rowMeans(west_array[ , , mode3_index, mode4_index], na.rm = TRUE) +
             colMeans(west_array[ , , mode3_index, mode4_index], na.rm = TRUE))/2
    }
}

## Fit Different Methods -----------------------------------------------------

loss_dat <- data.frame(Method = c("X", "ANOVA", "Logistic", "SoftHOSE", "TruncatedHOSVD"),
                       Loss = NA)
fnorm_dat <- data.frame(Method = c("X", "SoftHOSE", "TruncatedHOSVD"),
                        Fnorm = NA)

## X
loss_dat$Loss[loss_dat$Method == "X"] <- sum((west_array - west_true_as) ^ 2, na.rm = TRUE)

## logistic
predictors_temp <- expand.grid(lapply(p,seq,from = 1))
data_logistic <- data.frame(success = c(west_array_sub),
                            failure = c(west_array_n - west_array_sub),
                            bask = c(west_array),
                            mode1 = factor(predictors_temp[,1]),
                            mode2 = factor(predictors_temp[,2]),
                            mode3 = factor(predictors_temp[,3]),
                            mode4 = factor(predictors_temp[,4]))

west_glm <- glm(cbind(success,failure) ~ mode1 + mode2 + mode3 + mode4,
                family = binomial(), data = data_logistic)

which_not_diag <- rep(c(diag(p[1])), prod(p[3:4])) == 0
logistic_predict <- array(NA, dim = p)
logistic_predict[which_not_diag] <- predict(west_glm, type = "response")
logistic_predict_as <- array(arc_sin(logistic_predict, west_array_sub), dim = p)
loss_dat$Loss[loss_dat$Method == "Logistic"] <-
  sum((west_true_as - logistic_predict_as) ^ 2, na.rm = TRUE) ## loss is 1481

## anova
west_lm <- lm(bask ~ mode1 + mode2 + mode3 + mode4, data = data_logistic)
west_predict <- array(predict(west_lm), dim = p)
west_resid <- array(residuals(west_lm), dim = p)
loss_dat$Loss[loss_dat$Method == "ANOVA"] <-
  sum((west_true_as - west_predict) ^ 2, na.rm = TRUE) ## loss is 1344

## soft on residuals then add back
lambda_init <- c(0, 0, 0, 0)
c_init <- 1
resid_c <- get_c(west_resid)
trash <- utils::capture.output(soft_out <- soft_coord(resid_c, lambda_init,
                                                      c_init, print_iter = FALSE))
lambda_soft_final <- list(c(soft_out$lambda[1], soft_out$c),
                          soft_out$lambda[2],
                          soft_out$lambda[3],
                          soft_out$lambda[4])
soft_west_final_resid <- sure_given_c(obj = resid_c,
                                      func_lasso,
                                      dfunc_lasso,
                                      lambda_soft_final)
## soft_west_final_resid$sure_val

soft_final_est <- soft_west_final_resid$mean_est + west_predict

loss_dat$Loss[loss_dat$Method == "SoftHOSE"] <-
  sum((west_true_as - soft_final_est) ^ 2, na.rm = TRUE) ## loss is 1327

fnorm_dat$Fnorm[fnorm_dat$Method == "X"] <- sqrt(sum(west_resid ^ 2)) ## fnorm of resids is 38.95
fnorm_dat$Fnorm[fnorm_dat$Method == "SoftHOSE"] <- sqrt(sum(soft_west_final_resid$mean_est ^ 2)) ## fnorm of shrunken resids is 6.99

## low rank resids then add back
trunc_out <- sure_rank(west_resid)
est_rank <- lapply(trunc_out$min_rank, c) ## rank of 4, 2, 1, 1
trunc_west_final_resid <- sure_given_c(obj = resid_c,
                                       func_trunc,
                                       dfunc_trunc,
                                       est_rank)
## trunc_west_final_resid$sure_val

trunc_final_est <- trunc_west_final_resid$mean_est + west_predict
loss_dat$Loss[loss_dat$Method == "TruncatedHOSVD"] <- sum((west_true_as - trunc_final_est) ^ 2, na.rm = TRUE) ## loss is 1391
fnorm_dat$Fnorm[fnorm_dat$Method == "TruncatedHOSVD"] <- sqrt(sum(trunc_west_final_resid$mean_est ^ 2))


write.csv(loss_dat, file = "./output/nba_results/west_loss.csv", row.names = FALSE)
write.csv(fnorm_dat, file = "./output/nba_results/west_fnorm.csv", row.names = FALSE)
