## In this file, we plot the results from the simulation study.

library(ggplot2)
library(tidyr)
suppressMessages(library(dplyr))
library(xtable)

## Read in simulations results -----------------------------------------------
lstein   <- read.csv("./output/sims_out/losses_stein.csv")
lcandes2 <- read.csv("./output/sims_out/losses_candes2.csv")
lem      <- read.csv("./output/sims_out/losses_em.csv")
lcandes  <- read.csv("./output/sims_out/losses_candes.csv")
ltrunc   <- read.csv("./output/sims_out/losses_trunc.csv")
lsoft    <- read.csv("./output/sims_out/losses_soft.csv")

lstein$theta   <- "A" ## "uncorrelated"
lcandes$theta  <- "B" ## "dispersed_one_mode"
lem$theta      <- "C" ## "ar1_one_mode"
lcandes2$theta <- "D" ## "low_rank_one_mode"
lsoft$theta    <- "E" ## "evenly_dispersed_svs"
ltrunc$theta   <- "F" ## "low_multilinear_rank"

ldat <- select(lstein, 1:7, ncol(lstein)) %>%
  bind_rows(select(lcandes2, 1:7, ncol(lcandes2))) %>%
  bind_rows(select(lem, 1:7, ncol(lem))) %>%
  bind_rows(select(lcandes, 1:7, ncol(lcandes))) %>%
  bind_rows(select(ltrunc, 1:7, ncol(ltrunc))) %>%
  bind_rows(select(lsoft, 1:7, ncol(lsoft))) %>%
  gather(key = "Method", value = "Loss", contains("loss_"))
ldat$Method <- stringr::str_replace(string = ldat$Method, pattern = "loss_", replacement = "")
ldat$Method[ldat$Method == "x"] <- "X"
ldat$Method[ldat$Method == "soft"] <- "ST"
ldat$Method[ldat$Method == "trunc"] <- "Tr"
ldat$Method[ldat$Method == "em"] <- "EM"
ldat$Method[ldat$Method == "candes"] <- "MS"
ldat$Method[ldat$Method == "stein"] <- "JS"
ldat$Method[ldat$Method == "hooi"] <- "HO"
ldat$Method <- factor(ldat$Method, levels = c("ST", "Tr", "HO", "MS", "EM", "JS", "X"))

smalldat <- ldat %>% group_by(Method, theta) %>% summarise(Mean = mean(Loss)) %>% ungroup() %>%
  group_by(theta) %>% summarise(`Minimum Risk` = min(Mean))

pdf(file = "./output/figures/different_scenario_sims_losses_combined2.pdf",
    width = 6.5, height = 7, family = "Times")
ggplot(data = ldat, mapping = aes(x = Method, y = Loss)) +
  facet_wrap(~theta) +
  geom_boxplot() +
  theme_bw() +
  theme(strip.background = element_rect(fill="white")) +
  geom_hline(yintercept = 1000, lty = 2) + ## 1000 for the product of the dimension, theoretical risk of X
  geom_hline(data = smalldat, mapping = aes(yintercept = `Minimum Risk`), lty = 3)
dev.off()


#############################
## Just look at losses of the truncated methods -----------------------------
#############################
loss_dat_trunc <- select(lstein, loss_trunc, loss_cichocki:loss_bcv, ncol(lstein)) %>%
  bind_rows(select(lcandes2, loss_trunc, loss_cichocki:loss_bcv, ncol(lcandes2))) %>%
  bind_rows(select(lem, loss_trunc, loss_cichocki:loss_bcv, ncol(lem))) %>%
  bind_rows(select(lcandes, loss_trunc, loss_cichocki:loss_bcv, ncol(lcandes))) %>%
  bind_rows(select(ltrunc, loss_trunc, loss_cichocki:loss_bcv, ncol(ltrunc))) %>%
  bind_rows(select(lsoft, loss_trunc, loss_cichocki:loss_bcv, ncol(lsoft))) %>%
  gather(key = "Method", value = "Loss", contains("loss_"))
loss_dat_trunc$Method <- stringr::str_replace(loss_dat_trunc$Method, "loss_", "")
ggplot(loss_dat_trunc, mapping = aes(x = Method, y = Loss)) +
  geom_boxplot() +
  facet_wrap(~theta)


#############################
## Looking at how well the rank is estimated. --------------------------------
#############################
ltrunc$index <- 1:nrow(ltrunc)
ldat <- select(as_data_frame(ltrunc), contains("rank"), ncol(ltrunc)) %>%
  gather(key = "method_mode", value = "rank", sure_rank1:bcv_rank3) %>%
  separate(col = "method_mode", into = c("Method", "Mode"), sep = "_")
ldat$Mode <- stringr::str_replace(ldat$Mode, "rank", "")
sumdat <- group_by(ldat, index, Method) %>%
  mutate(is_five = (rank == 5)) %>%
  summarize(allfive = all(is_five)) %>%
  ungroup() %>%
  group_by(Method) %>%
  summarize(propfive = mean(allfive))
sumdat
names(sumdat)[2] <- "Proportion Correct"
sumdat <- filter(sumdat, Method != "bcv" & Method != "cichocki")
sumdat$Method <- c("MDL", "PA", "SURE")
cat(file = "./output/figures/table_prop.txt",
print(xtable(sumdat, label = "tab:prop.five",
             caption = "Proportion of times the multilinear rank is estimated correctly using either the minimum description length criterion (MDL), parallel analysis (PA), or Stein's unbiased risk estimate (SURE)."),
      include.rownames = FALSE)
)

## Now look at scenario D ----------------------------------------------------
ddat <- select(as_data_frame(lcandes2), contains("rank"))
prop_rank <- as_data_frame(apply(ddat, 2, function(x) { table(factor(x, levels = 0:10)) } ) / 500)
prop_rank$rank <- 0:(nrow(prop_rank) - 1)
sdat <- select(prop_rank, contains("sure"), contains("mdl"), contains("par"), rank) %>%
  gather(key = "method_mode", value = "prop", sure_rank1:par_rank3) %>%
  separate(method_mode, into = c("Method", "Mode"), sep = "_")
sdat$Mode <- stringr::str_replace(sdat$Mode, "rank", "")
sdat$Method <- stringr::str_replace(sdat$Method, "sure", "SURE")
sdat$Method <- stringr::str_replace(sdat$Method, "mdl", "MDL")
sdat$Method <- stringr::str_replace(sdat$Method, "par", "PA")

dummydat <- data_frame(Mode = c(1, 2, 3), true_rank = c(5, 10, 10))

pl <- ggplot(data = sdat, mapping = aes(x = rank, y = prop, color = Method, lty = Method)) +
  geom_line() +
  facet_wrap(~Mode) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  ylab("Proportion") +
  xlab("Estimated Rank") +
  scale_x_continuous(breaks = 0:10) +
  geom_vline(data = dummydat, mapping = aes(xintercept = true_rank), lty = 2) +
  ggthemes::scale_color_colorblind()

pdf(file = "output/figures/est_rank.pdf", colormodel = "cmyk", family = "Times", width = 6.5, height = 2.2)
print(pl)
dev.off()
