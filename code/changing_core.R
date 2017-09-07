#####################
## Filename: changing_core.r
## Created by: David Gerard
## Synopsis: Generate a core that has low multilinear rank but full rank along one mode.
##           See what the higher-order spectral functions are actually doing to the HOSVD.
##           Then I explore how the mode-specific SV's change when you apply a HOSE to them.
#####################

library(hose)
library(tensr)
library(ggplot2)
library(tidyr)
suppressMessages(library(dplyr))
source("./code/misc_functions.R")

#####################################
## generate an array that has very high rank in one dimension but very low in rest
#####################################
p <- c(10, 10, 10)
r <- c(10, 4, 3)
fnorm_mean <- 10000
x_temp <- array(stats::rnorm(prod(r)), dim = r)
s_temp <- array(0, dim = p)
indices_s_temp <- as.matrix(expand.grid(lapply(r, seq, from = 1)))
s_temp[indices_s_temp] <- x_temp
U_mult <- tensr::hosvd(array(stats::rnorm(prod(p)), dim = p))$U
theta <- tensr::atrans(s_temp / tensr::fnorm(s_temp) * fnorm_mean, U_mult)

trash <- utils::capture.output(theta <- generate_core_same(p = p, r = r, fnorm_mean = 10))
hosvd_theta <- hosvd_full(theta)

names(hosvd_theta$D) <- c(1, 2, 3)
dat <- as_data_frame(hosvd_theta$D) %>% dplyr::bind_cols(data_frame(Index = 1:p[1])) %>%
  tidyr::gather(key = "Mode", value = "Singular Values", 1:3)
dat$Mode <- as.numeric(dat$Mode)
dat$X = dat$Index + dat$Mode / 9 - 2 / 9
dat$Mode <- as.factor(dat$Mode)


pdf(file = "./output/figures/low_multilinear_rank.pdf", family = "Times", height = 1.4, width = 4.2)
ggplot(data = dat, mapping = aes(x = X, y = `Singular Values`, linetype = Mode)) +
  geom_segment(aes(xend = X, yend = 0)) +
  theme_bw() +
  ylab("Value") +
  xlab("Singular Value") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
dev.off()

#################
## See how the HOSVD changes with mode specific soft-thresholding
#################
func_lasso <- list()
dfunc_lasso <- list()
for(mode_index in 1:length(p)) {
    func_lasso[[mode_index]] <- f_lasso
    dfunc_lasso[[mode_index]] <- df_lasso
}

set.seed(129)
p = c(10, 10, 10)
theta <- array(stats::rnorm(prod(p)), dim = p)
theta <- tensr::atrans(theta, list(diag(seq(10, 8, length = 10)),
                                   diag(seq(10, 8, length = 10)),
                                   diag(seq(10, 8, length = 10))))
theta <- theta / 1000
hosvd_theta <- hose::hosvd_full(theta)

c_x <- hose::get_c(theta)
final_theta <- hose::sure_given_c(obj = c_x, func = func_lasso,
                                  dfunc = dfunc_lasso,
                                  lambda = list(5, 0, -100))
hosvd_final <- hose::hosvd_full(final_theta$mean_est)

dat <- data_frame(Mode = 1, Shrinkage = "Before",
                  SV = hosvd_theta$D[[1]] / sum(hosvd_theta$D[[1]]),
                  Index = 1:p[1]) %>%
  bind_rows(data_frame(Mode = 1, Shrinkage = "After",
                       SV = hosvd_final$D[[1]] / sum(hosvd_final$D[[1]]),
                       Index = 1:p[1])) %>%
  bind_rows(data_frame(Mode = 2, Shrinkage = "Before",
                       SV = hosvd_theta$D[[2]] / sum(hosvd_theta$D[[2]]),
                       Index = 1:p[2])) %>%
  bind_rows(data_frame(Mode = 2, Shrinkage = "After",
                       SV = hosvd_final$D[[2]] / sum(hosvd_final$D[[2]]),
                       Index = 1:p[2])) %>%
  bind_rows(data_frame(Mode = 3, Shrinkage = "Before",
                       SV = hosvd_theta$D[[3]] / sum(hosvd_theta$D[[3]]),
                       Index = 1:p[3])) %>%
  bind_rows(data_frame(Mode = 3, Shrinkage = "After",
                       SV = hosvd_final$D[[3]] / sum(hosvd_final$D[[3]]),
                       Index = 1:p[3]))
dat$X <- dat$Index + (dat$Shrinkage == "Before") / 3 - 1/6
dat$Mode <- as.factor(paste0("Mode ", dat$Mode))
dat$Shrinkage <- factor(dat$Shrinkage, levels = c("Before", "After"))

pdf(file = "./output/figures/effects_on_svs.pdf", width = 4.5, height = 1.4, family = "Times")
ggplot(data = dat, mapping = aes(x = X, y = SV, linetype = Shrinkage)) +
  geom_segment(aes(xend = X, yend = 0), lwd = 0.3) +
  facet_wrap(~Mode) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) +
  ylab("Normalized\nSingular Values") +
  theme(strip.background = element_rect(fill="white"))
dev.off()
