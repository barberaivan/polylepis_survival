# Packages ----------------------------------------------------------------

library(tidyverse)
library(viridis)
library(readxl)
library(posterior)
library(rstan)
library(DHARMa)
library(tidybayes)
library(terra) # to check spatial correlation
theme_set(theme_bw())

# Functions ---------------------------------------------------------------

# function to estimate age from tree perimeter (cm) and rock cover (%).
tree_age <- function(perim, rock) {
  log_age <- - 0.16 + (0.85 * log10(perim)) + (0.0013 * rock)
  return(10 ^ log_age)
}

mean_ci <- function(x) {
  qq <- quantile(x, probs = c(0.025, 0.975), method = 8) %>% unname
  return(c("mean" = mean(x), "lower" = qq[1], "upper" = qq[2]))
}

# Import data -------------------------------------------------------------

data_name <- "Plots_remedicion 2024.03.10.xls"

# plot-level data
dplot <- read_excel(file.path("..", data_name), sheet = "plotsall_use")

# tree-level data
dtree <- read_excel(file.path("..", data_name), sheet = "taball_use")

# subset plots to use, with reliable data
dplot <- dplot[is.na(dplot$`delete?`), ] # with no comment in "delete?"
dtree <- dtree[dtree$plot %in% dplot$plot, ]


# Tidy data ---------------------------------------------------------------

# rename a few columns
dplot <- rename(dplot, manag = "manag2", fire = "fire01")
dtree <- rename(dtree, surv = "sur2018", fire = "fireplot01",
                height = "hei2003", perim = "peri2003", tree_id = "arbol")

dtree <- dtree[, c("plot","surv", "height", "perim", "tree_id", "prock")]
dtree$surv <- as.numeric(dtree$surv)

# make factors
dtree$plot <- factor(dtree$plot)
dplot$plot <- factor(dplot$plot)
dplot$basin <- factor(dplot$basin)
dplot$fire <- factor(as.character(dplot$fire), levels = c("0", "1"),
                     labels = c("Unburned plots", "Burned plots"))
dplot$manag <- factor(dplot$manag, levels = c("a", "b"),
                      labels = c("rangeland", "livestock exclusion"))

# add age
dtree$perim <- as.numeric(dtree$perim)
dtree$age <- tree_age(dtree$perim, dtree$prock)

# add plot variables to tree data
dplot_sub <- dplot[, c("plot", "manag", "fire", "elev", "pforest", "basin")]
dtree <- left_join(dtree, dplot_sub, by = "plot")

# resolve duplicated tree_ids
dtree$tree_id <- paste(dtree$plot, dtree$tree_id, sep = "_")
dtree$tree_id[duplicated(dtree$tree_id)] <-
  paste(dtree$tree_id[duplicated(dtree$tree_id)], "b", sep = "")
dtree$tree_id <- factor(dtree$tree_id)

# standardize continuous predictors, but save means and sd
plot_means <- aggregate(cbind(age, height, elev, pforest) ~ plot, dtree, mean)

pred_means <- apply(as.matrix(plot_means[, -1]), 2, mean)
pred_sds <- apply(as.matrix(plot_means[, -1]), 2, sd)

dtree$height_z <- (dtree$height - pred_means["height"]) / pred_sds["height"]
dtree$age_z <- (dtree$age - pred_means["age"]) / pred_sds["age"]
dtree$elev_z <- (dtree$elev - pred_means["elev"]) / pred_sds["elev"]
dtree$pforest_z <- (dtree$pforest - pred_means["pforest"]) / pred_sds["pforest"]

# remove rows with NA in any useful predictor
out <- is.na(dtree$height) | #is.na(dtree$age) |
       is.na(dtree$plot) | is.na(dtree$elev) | is.na(dtree$pforest)
dtree <- dtree[!out, ]
nrow(dtree) # 1229 available cases

# Identify complete and incomplete plots -----------------------------------

dn <- data.frame(plot = dtree$plot, surv = dtree$surv)
dtagg <- aggregate(surv ~ plot, dn, mean)
names(dtagg)[2] <- "meansurv"

dtagg$na <- numeric(nrow(dtagg))
dtagg$length <- numeric(nrow(dtagg))
for(i in 1:nrow(dtagg)) {
  dtagg$na[i] <- sum(is.na(dtree$surv[dtree$plot == dtagg$plot[i]]))
  dtagg$length[i] <- length(dtree$surv[dtree$plot == dtagg$plot[i]])
}

# all(dtagg$length > dtagg$na) # OK
dtagg <- left_join(dtagg, dplot[, c("plot", "ntrees", "nalive", "ndead")], by = "plot")
dtagg <- dtagg[complete.cases(dtagg), ]

dtagg[dtagg$ntrees < dtagg$length, ]
dtagg[dtagg$ntrees < dtagg$length, ] %>% nrow - 1 # 46 plots

dtagg$miss <- dtagg$length - dtagg$ntrees
sum(dtagg$miss) # 349 missing trees.
sum(dtagg$miss) / sum(dtagg$length) # 28.5 % missing

plot(miss ~ length, dtagg); abline(0, 1) # ok, always miss < total
plot(miss ~ na, dtagg); abline(0, 1) # ok, always na >= miss

# evaluate how many counted trees have no matching id.
# counted are ntrees, and identified are length - na
dtagg$identified <- dtagg$length - dtagg$na
dtagg$unidentified <- dtagg$ntrees - dtagg$identified

# see those plots
filter <- dtagg$unidentified > 0
(dsub <- dtagg[filter, ])

# compute number of unidentified dead and alive
dsub$ndead_unid <- numeric(nrow(dsub))
dsub$nalive_unid <- numeric(nrow(dsub))
for(i in 1:nrow(dsub)) {
  # i = 1
  xx <- dtree$surv[dtree$plot == dsub$plot[i]]
  xx <- na.omit(xx)

  ndead_id <- length(xx) - sum(xx)
  nalive_id <- sum(xx)

  dsub$ndead_unid[i] <- dsub$ndead[i] - ndead_id
  dsub$nalive_unid[i] <- dsub$nalive[i] - nalive_id
}

plot(ndead_unid ~ nalive_unid, dsub); abline(0, 1) # most unid are dead
# in cases where only dead or only alive trees are unidentified,
# the problem is much smaller.
sum(dsub$unidentified) # 35 trees
sum(dsub$unidentified) / sum(dtagg$ntrees) * 100
# 4 % of non-missing trees are unidentified

# list with ones and zeroes for unidentified trees
bin_list <- lapply(1:nrow(dsub), function(r) {
  c(rep(0, dsub$ndead_unid[r]), rep(1, dsub$nalive_unid[r]))
})

rows_list <- lapply(1:nrow(dsub), function(r) {
  dtree_local <- dtree[dtree$plot == dsub$plot[r], ]
  x <- dtree_local$tree_id[is.na(dtree_local$surv)]
  return(x)
  # return(unique(x))
})

# make non-redundant combinations
combs <- lapply(1:nrow(dsub), function(i) {
  # i = 7
  y <- bin_list[[i]]
  arbol <- rows_list[[i]]

  if(length(y) == 1) {
    combs <- matrix(arbol, nrow = 1)
    return(combs)
  }

  dead <- sum(y == 0)
  alive <- sum(y)

  combs_dead <- combn(arbol, dead)
  combs_alive <- combn(arbol, alive)

  cols_merge <- expand.grid(dd = 1:ncol(combs_dead),
                            aa = 1:ncol(combs_alive))

  combs_full <- matrix(NA, nrow = length(y), nrow(cols_merge))
  for(j in 1:nrow(cols_merge)) {
    # j = 1
    combs_full[, j] <- c(combs_dead[, cols_merge$dd[j]],
                         combs_alive[, cols_merge$aa[j]]) %>% as.character
  }
  unique_len <- apply(combs_full, 2, function(x) length(unique(x)))

  combs <- combs_full[, unique_len == length(y)]
  return(combs)
})

dsub$combinations_unique <- sapply(combs, ncol)
# In these plots we will compute the marginal likelihood with respect to all
# possible id combinations.

# https://mbjoseph.github.io/posts/2020-04-28-a-step-by-step-guide-to-marginalizing-over-discrete-parameters-for-ecologists-using-stan/

# Whithin each plot, for each possible id combination, compute the sum of
# binomial_logit_lpmf. Then, apply log_sum_exp to that vector of likelihoods,
# to get the likelihood of the plot marginal to which were the unidentified
# trees.

# dtagg has data at the parcel level indicating how many trees there are,
# how many na, and how many identified.
# for plots where unidentified == 0, we jus remove surv == NA.
# For plots with unidentified > 0, we separate the NA and the obverved.
dtagg %>% dim()
sum(dtagg$unidentified > 0) # 12 plots == nrow(dsub)
sum(dtagg$unidentified == 0) # 126 normal plots, where NA will be removed.

# get plots with no unidentified trees
data_id <- dtree[dtree$plot %in% dtagg$plot[dtagg$unidentified == 0], ]
# remove NAs
data_id <- data_id[complete.cases(data_id[, c("height", "elev", "pforest",
                                              "manag", "fire", "surv")]), ]

# get plots with unidentified trees
data_unid <- dtree[dtree$plot %in% dsub$plot, ]
# non-na trees
data_unid_obs <- data_unid[!is.na(data_unid$surv), ]
# na trees
data_unid_na <- data_unid[is.na(data_unid$surv), ]
nrow(data_unid_na) # 35 were recorded but unidentified, and there are 123 options.

# merge non-na trees
data_obs <- rbind(data_id, data_unid_obs)
nrow(data_obs) # 841

anyNA(data_obs[, -which(names(data_obs) %in% c("surv", "age", "perim", "age_z"))])
# View(data_obs) # OK

anyNA(data_unid_na[, -which(names(data_unid_na) %in% c("surv", "age", "age_z", "perim"))])
# View(data_unid_na)
# OK, but there is one NA in perim at unid_na


# Note about age data -----------------------------------------------------

# There is one NA at perim2003. If age is used, this code might break.

# Data for Stan -----------------------------------------------------------

# make three design matrices for each dataset:
# fixed effects (X), plot (Zp), basin (Zb).
# Suffix _unid will indicate the matrix for unidentified trees
fixed_formula <- formula(~ manag + fire +
                           height_z + elev_z + I(elev_z ^ 2) + pforest_z)
plot_formula <- formula(~ - 1 + plot)
basin_formula <- formula(~ - 1 + basin)

X <- model.matrix(fixed_formula, data = data_obs)
X_unid <- model.matrix(fixed_formula, data = data_unid_na)

Zp <- model.matrix(plot_formula, data = data_obs)
Zp_unid <- model.matrix(plot_formula, data = data_unid_na)

Zb <- model.matrix(basin_formula, data = data_obs)
Zb_unid <- model.matrix(basin_formula, data = data_unid_na)

# data for marginalization of unidentified trees

# this is the response variable, which has to be assigned to different combinations
# of unidetified trees
y_unid <- do.call("c", bin_list) # ahora sí son 35
y_unid_length <- sapply(bin_list, length)
y_unid_end <- cumsum(y_unid_length)
y_unid_start <- y_unid_end - y_unid_length + 1
p_unid <- length(y_unid_length)

# these are the unidentified trees
data_unid_agg <- aggregate(elev ~ plot, data_unid_na, length)
d_unid_length <- data_unid_agg$elev
d_unid_end <- cumsum(d_unid_length)
d_unid_start <- d_unid_end - d_unid_length + 1

# the probability will be computed in all data points, and then,
# different subsets of probability will be used in the bernoulli_logit_lpmf.
max_combs <- max(sapply(combs, ncol))

combs_rows_matrix <- do.call("rbind", lapply(1:p_unid, function(p) {
  # p = 1
  # print(p)
  cc <- combs[[p]]
  # tree_names <- rows_list[[p]]
  mm <- matrix(999999, y_unid_length[p], max_combs) # fill with 999999 so Stan works
  for(j in 1:ncol(cc)) {
    for(i in 1:nrow(cc)) {
      # i = 1; j = 1
      mm[i, j] <- which(data_unid_na$tree_id == cc[i, j])
    }
  }
  # mm[, 1:10]
  return(mm)
}))
# In Stan, this should be an integer array, used to subset the rows of prob_surv
# to compute the loglik.

# # test code:
# probs <- runif(nrow(data_unid_na)) # fitted probs
# marginal_logliks <- numeric(p_unid)
# for(p in 1:p_unid) {
#   # p = 1
#   nc <- ncol(combs[[p]])
#   # print(nc)
#   rows_response <- y_unid_start[p] : y_unid_end[p]
#   response <- y_unid[rows_response]
#
#   # rows_pred <- d_unid_start[p] : d_unid_end[p]
#
#   loglik_combs <- numeric(nc)
#
#   for(c in 1:nc) {
#     # print(c)
#     rows_pred <- combs_rows_matrix[rows_response, c]
#     pp <- probs[rows_pred]
#     loglik_combs[c] <- dbinom(response, prob = pp, size = 1) %>% log %>% sum
#   }
#
#   marginal_logliks[p] <- log(sum(exp(loglik_combs)))
# }
# # target += marginal_loglikds

stan_data <- list(
  N = nrow(data_obs),
  N_unid_rows = nrow(data_unid_na),
  N_unid_trees = length(y_unid),

  K = ncol(X),
  P = ncol(Zp),
  P_unid = length(combs),
  B = ncol(Zb),

  y = data_obs$surv,
  y_unid = y_unid,

  X = X, X_unid = X_unid,
  Zp = Zp, Zp_unid = Zp_unid,
  Zb = Zb, Zb_unid = Zb_unid,

  y_unid_start = y_unid_start,
  y_unid_end = y_unid_end,
  y_unid_length = y_unid_length,

  max_combs = max(sapply(combs, ncol)),
  combs_rows_matrix = combs_rows_matrix,
  combs_n = sapply(combs, ncol),

  unid = 1, # use unidentified trees in likelihood.

  prior_intercept_sd = 5,
  prior_b_sd = 3,
  prior_sigma_sd = 5 #3 # sigma_basin copies the prior (sd = 3). try larger.
)


# Model fit ---------------------------------------------------------------

stan_code <- stan_model("survival_with_basin.stan", verbose = T)

# # consider unidentified trees
#
# running with prior_sd = 5 in sigma
m1 <- sampling(stan_code, data = stan_data, seed = 1234, refresh = 10,
               # cores = 1, chains = 1, iter = 10, ## test
               cores = 6, chains = 6, iter = 2000,
               control = list(adapt_delta = 0.9))
saveRDS(m1, file.path("files", "survival_model_samples_with_basin_sd5.rds"))
# # 734.973 / 60 = 12.24 min
m1 <- readRDS(file.path("files", "survival_model_samples_with_basin.rds"))
m1 <- readRDS(file.path("files", "survival_model_samples_with_basin_sd5.rds"))
sm1 <- summary(m1)[[1]]
max(sm1[, "Rhat"]); min(sm1[, "n_eff"]) # nice
# [1] 1.005346
# [1] 1413.106

# ingnore unidentified trees
stan_data2 <- stan_data
stan_data2$unid <- 0
m2 <- sampling(stan_code, data = stan_data2, seed = 1234, refresh = 100,
               # cores = 1, chains = 1, iter = 10,
               cores = 6, chains = 6, iter = 2000,
               control = list(adapt_delta = 0.9, max_treedepth = 20))
sm2 <- summary(m2)[[1]] # anduvo!
max(sm2[, "Rhat"]); min(sm2[, "n_eff"]) # nice


## model without basin (it's hard to estimate sigma)
## consider unidentified trees
# stan_code_nb <- stan_model("survival.stan", verbose = T)
# m3 <- sampling(stan_code_nb, data = stan_data, seed = 1234, refresh = 10,
#                cores = 6, chains = 6, iter = 2000,
#                control = list(adapt_delta = 0.9))
# saveRDS(m3, file.path("files", "survival_model_samples_nobasin.rds"))
## 530.419 / 60 = 8.84 min
m3 <- readRDS(file.path("files", "survival_model_samples.rds"))
sm3 <- summary(m3)[[1]]
max(sm3[, "Rhat"]); min(sm3[, "n_eff"]) # nice
# 1.005588
# 1096.096

# parcelas con corr espacial?

# compare posteriors
bnames <- colnames(X)
par_names <- c(bnames, "sigma_plot", "sigma_basin")
bm1 <- as.matrix(m1, par = c("b", "sigma_plot", "sigma_basin"))
bm2 <- as.matrix(m3, par = c("b", "sigma_plot", "sigma_basin"))
colnames(bm1) <- colnames(bm2) <- par_names

# plots to compare
par(mfrow = c(3, 3))
for(j in 1:length(par_names)) {
  # j = 9

  if(par_names[j] %in% c("sigma_plot", "sigma_basin")) {
    d1 <- density(bm1[, j], from = 0)
    d2 <- density(bm2[, j], from = 0)
  } else {
    d1 <- density(bm1[, j])
    d2 <- density(bm2[, j])
  }

  xl <- range(c(d1$x, d2$x))
  yl <- c(0, max(c(d1$y, d2$y)))

  plot(d1, col = "black", ylim = yl * 1.01, xlim = xl * 1.01,
       main = NA, xlab = par_names[j], ylab = "density")
  lines(d2, col = "blue")
}
par(mfrow = c(1, 1))
# parameters considerably.

# the basin sd copies the prior.
# explore the raneffs
e_basin <- as.matrix(m1, "e_basin")
sigma_basin <- as.matrix(m1, "sigma_basin") %>% as.numeric

prior_raneff <- rnorm(length(sigma_basin) * 3, 0, sd = sigma_basin)
plot(density(prior_raneff), lwd = 2, ylim = c(0, 1), xlim = c(-10, 10))
for(i in 1:5) {
  lines(density(e_basin[, i]), col = i+1)
}
# although the sd copies the prior, raneffs are quite different.

# increasing the sd sd to 5 changes the prior, because the sigma copies the prior,
# but the raneff are estimated the same.

# Predictions -------------------------------------------------------------

manag_lev <- levels(dplot$manag)#c("Rangeland", "Livestock exclusion")
fire_lev <- levels(dplot$fire)#c("Unburned", "Burned")

pdata <- rbind(
  # elevation
  expand.grid(
    manag = factor(manag_lev, levels = manag_lev),
    fire = factor(fire_lev, fire_lev),
    elev_z = seq(min(dtree$elev_z), max(dtree$elev_z), length.out = 150),
    pforest_z = 0,
    height_z = 0,
    varying_var = "Elevation (m a.s.l.)"
  ),
  # pforest
  expand.grid(
    manag = factor(manag_lev, levels = manag_lev),
    fire = factor(fire_lev, fire_lev),
    elev_z = 0,
    pforest_z = seq(min(dtree$pforest_z), max(dtree$pforest_z), length.out = 150),
    height_z = 0,
    varying_var = "Forest cover (%)"
  ),
  # height
  expand.grid(
    manag = factor(manag_lev, levels = manag_lev),
    fire = factor(fire_lev, fire_lev),
    elev_z = 0,
    pforest_z = 0,
    height_z = seq(min(dtree$height_z), max(dtree$height_z), length.out = 150),
    varying_var = "Tree height (cm)"
  )
)

# tidy factors

# add variables at original scale
pdata$varying_val <- 0
pdata$manag_plot <- factor(pdata$manag, manag_lev,
                           labels = c("Rangeland", "Livestock exclusion"))

rows_replace <- pdata$varying_var == "Elevation (m a.s.l.)"
pdata$varying_val[rows_replace] <-
  pdata$elev_z[rows_replace] * pred_sds["elev"] + pred_means["elev"]

rows_replace <- pdata$varying_var == "Forest cover (%)"
pdata$varying_val[rows_replace] <-
  pdata$pforest_z[rows_replace] * pred_sds["pforest"] + pred_means["pforest"]

rows_replace <- pdata$varying_var == "Tree height (cm)"
pdata$varying_val[rows_replace] <-
  pdata$height_z[rows_replace] * pred_sds["height"] + pred_means["height"]

Xpred <- model.matrix(fixed_formula, pdata)

# choose which model to use!
bhat <- as.matrix(m1, "b") %>% t
# bhat <- as.matrix(m2, "b") %>% t
bhat <- as.matrix(m3, "b") %>% t

phat <- plogis(Xpred %*% bhat)
psumm <- apply(phat, 1, mean_ci) %>% t %>% as.data.frame %>% "*"(100)

predictions <- cbind(pdata, psumm)

vv <- unique(pdata$varying_var)
# v2 <- c("elev", "pforest", "")
plist <- vector("list", 3)

for(v in 1:3) {
  # v = 1
  d <- predictions[predictions$varying_var == vv[v], ]

  # add poitns
  # dp <- dplot
  # dp$x_class <- cut(dp[, ])

  pp <-
  ggplot(d, aes(x = varying_val, y = mean, ymin = lower, ymax = upper,
                          color = fire, fill = fire)) +
    geom_ribbon(color = NA, alpha = 0.25) +
    geom_line() +
    scale_color_viridis(discrete = T, option = "B", end = 0.5) +
    scale_fill_viridis(discrete = T, option = "B", end = 0.5) +
    facet_wrap(~ manag_plot, ncol = 2) +
    scale_y_continuous(limits = c(0, 100), expand = c(0.01, 0.01)) +
    xlab(vv[v]) +
    ylab("Survival (%)") +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(color = "white", fill = "white"),
          panel.spacing = unit(3, "mm"),
          axis.title.y = element_blank(),
          legend.title = element_blank(),
          legend.position = "none",
          axis.text = element_text(size = 9),
          axis.title = element_text(size = 10),
          strip.text = element_text(size = 10))

  pp
  if(v > 1) {
    pp <- pp + theme(strip.background = element_blank(),
                     strip.text = element_blank())

  }

  if(v == 2) {
    pp <- pp + theme(legend.position = "right",
                     axis.title.y = element_text())
  }

  plist[[v]] <- pp
}

surv_plot <- egg::ggarrange(plots = plist, nrow = 3)
ggsave("figures/survival_prediction_all_data_nobasin.png", plot = surv_plot,
       width = 14, height = 15, units = "cm")
# ggsave("figures/survival_prediction_found_tags.png", plot = surv_plot,
#        width = 14, height = 15, units = "cm")



# Residuals analyses ------------------------------------------------------

# define model to change easily
model <- m3
basin <- FALSE

# extract parameters
bhat <- as.matrix(model, "b") %>% t
sigma_plot <- as.matrix(model, "sigma_plot") %>% as.numeric
if(basin) sigma_basin <- as.matrix(model, "sigma_basin") %>% as.numeric
nsim <- length(sigma_plot)

# simulate data for new basins and plots, to compute DHARMa residuals.

# simulation for for identified trees
e_plot_sim <- matrix(rnorm(nrow(X) * nsim, sd = sigma_plot),
                     nrow(X), nsim)
if(basin) {
  e_basin_sim <- matrix(rnorm(nrow(X) * nsim, sd = sigma_basin),
                        nrow(X), nsim)
  eta <- X %*% bhat + e_plot_sim + e_basin_sim
} else {
  eta <- X %*% bhat + e_plot_sim
}

psim <- plogis(eta)
ysim <- apply(psim, 2, function(x) rbinom(nrow(X), prob = x, size = 1))

# In the case of unidentified trees, for each iteration, we compute
# the likelihood of all combinations, sample one according to the likelihood,
# and record their height value. Heights will then be averaged across replicates.

# extract stan variables, as the likelihood of all possible combinations will be
# computed borrowing the stan code.
N_unid_trees <- stan_data$N_unid_trees
N_unid_rows <- stan_data$N_unid_rows
P_unid <- stan_data$P_unid # number of plots with unid trees
combs_n <- stan_data$combs_n

# simulated data for unidentified trees
ysim_unid <- matrix(NA, N_unid_trees, nsim)

# height used to simulate the corresponding data
height_unid <- matrix(NA, N_unid_trees, nsim)

e_plot_sim <- matrix(rnorm(nrow(X_unid) * nsim, sd = sigma_plot),
                     nrow(X_unid), nsim)

if(basin) {
  e_basin_sim <- matrix(rnorm(nrow(X_unid) * nsim, sd = sigma_basin),
                        nrow(X_unid), nsim)
  eta <- X_unid %*% bhat + e_plot_sim + e_basin_sim
} else {
  eta <- X_unid %*% bhat + e_plot_sim
}

psim_options <- plogis(eta)
height_options <- data_unid$height_z

# loop over plots with unidentified trees
for(i in 1:nsim) {
  print(i)
  for(p in 1:P_unid) {
    # p = 1
    nc <- combs_n[p]
    ylen <- y_unid_length[p]
    rows_response <- y_unid_start[p] : y_unid_end[p]
    like_combs <- numeric(nc)

    # loglik for every possible combination at plot p
    for(c in 1:nc) {
      rows_pred = combs_rows_matrix[rows_response, c]
      like_combs[c] = dbinom(y_unid[rows_response], size = 1,
                             prob = psim_unid[rows_pred, i], log = TRUE) %>% sum %>% exp
    }

    # sample a combination of trees identities based on the likelihood
    chosen_comb <- sample(1:nc, size = 1, prob = like_combs)

    # identify to which rows in the data of unidentified trees (N_unid_rows = 123)
    # correspond the chosen options
    rows_pred <- combs_rows_matrix[rows_response, chosen_comb]

    # get probabilities and simulate
    psim_chosen <- psim_unid[rows_pred, i]
    ysim_unid[rows_response, i] <- rbinom(ylen, prob = psim_chosen, size = 1)

    # get associated predictors (height is the most important)
    height_unid[rows_response, i] <- height_options[rows_pred]
  }
}

# merge the unid data with the remaining predictors
data_unid_y <- data_unid[!duplicated(data_unid$plot), ]
rows_rep <- rep(1:P_unid, y_unid_length)
data_unid_y <- data_unid_y[rows_rep, ]
data_unid_y$height_z <- rowMeans(height_unid)
data_unid_y$height <- data_unid_y$height_z * pred_sds["height"] + pred_means["height"]
data_unid_y$surv <- y_unid

# merge data_obs with data_unid_y, and the corresponding ysim
data_unid_y <- data_unid_y[, colnames(data_obs)]
data_both <- rbind(data_obs, data_unid_y)
ysim_both <- rbind(ysim, ysim_unid)

# residuals
res <- createDHARMa(observedResponse = data_both$surv,
                    simulatedResponse = ysim_both, integerResponse = T)
data_both$res <- res$scaledResiduals

plot(res, rank = F)
plotResiduals(res, form = data_both$height)
plotResiduals(res, form = data_both$elev)
plotResiduals(res, form = data_both$pforest)
plotResiduals(res, form = data_both$manag)
plotResiduals(res, form = data_both$fire)
plotResiduals(res, form = data_both$basin)
plotResiduals(res, form = data_both$plot)
ggplot(data_both, aes(plot, res)) + geom_boxplot()

res_data <- list(res = res, data_both = data_both)
saveRDS(res_data, file.path("files", "survival_model_residuals_and_data.rds"))

# Spatial correlation? ----------------------------------------------------

# plots usable for residuals ("417" was lost)
plots_unique <- as.character(unique(data_both$plot)) %>% sort

# load waypoints
vplot <- vect(file.path("..", "plots_waypoints.shp"))
# plot(vplot)
vplot <- vplot[vplot$Name %in% plots_unique, ]
vplot <- vplot[order(vplot$Name), ]
# all(plots_unique == vplot$Name) # OK

spat_dist <- distance(vplot, symmetrical = TRUE) %>% as.matrix

hist(spat_dist[lower.tri(spat_dist)] / 1000, breaks = 30,
     main = "Pairwise distance distribution (km)")

# assess spat corr in plots random effects
e_plot_hat <- as.matrix(m1, par = c("e_plot")) %>% t
e_plot_hat <- e_plot_hat[levels(dtree$plot) %in% plots_unique, ]

e_plot_means <- apply(e_plot_hat, 1, mean)
e_plot_dist <- as.matrix(dist(e_plot_means))

# residuals diff
res_agg <- aggregate(res ~ plot, data_both, mean)
res_dist <- as.matrix(dist(res_agg$res))

dist_df <- data.frame(raneff_diff = e_plot_dist[lower.tri(e_plot_dist)],
                      res_diff = res_dist[lower.tri(res_dist)],
                      distance = spat_dist[lower.tri(spat_dist)])

ggplot(dist_df, aes(distance, raneff_diff)) +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5),
              method.args = list(family = Gamma(link = "log")),
              color = "red")

ggplot(dist_df, aes(distance, res_diff)) +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5),
              method.args = list(family = Gamma(link = "log")),
              color = "red")


# plot random effects variation among basins?
dplot_basins <- dplot[dplot$plot %in% levels(dplot$plot), ]
dplot_basins <- dplot_basins[order(dplot_basins$plot), ]
raneff_basin <- data.frame(raneff_plot = e_plot_means,
                           basin = dplot_basins$basin)

mran <- lm(raneff_plot ~ basin, raneff_basin)
car::Anova(mran)
summary(mran)

ggplot(raneff_basin, aes(x = basin, y = raneff_plot)) +
  geom_boxplot()
# parece que sí hay un efecto de cuenca no despreciable, pero no parece haber
# corr espacial.
# De todos modos, habría que chequear la corr espacial calculando los residuos
# dharma marginales a los raneff.

# habrá problemas para identificar efectos de cuencas?
table(data_unid$basin) # sólo hay unid en santa clara (c) y en el hueco (d)
table(dplot$basin, dplot$fire)
# las c y e (s clara y hueco) son las que más fuego tienen.

# habrá corr entre los eff de fuego y cuenca?
basin_names <- aggregate(as.numeric(as.character(plot)) ~ basin, dplot, mean)
names(basin_names)[2] <- "plot_ids"
basin_names$name <- c("Los Gigantes", "Río Condorito", "Santa Clara",
                      "Río Mina Clavero", "El Hueco")


# Fire and management -----------------------------------------------------

table(dplot$fire, dplot$manag)
dplot$fire_bin <- as.numeric(dplot$fire) - 1
mmmm <- glm(fire_bin ~ manag, data = dplot, family = "binomial")
summary(mmmm)
# pvalue 0.186
mmmm <- glm(fire_bin ~ manag - 1, data = dplot, family = "binomial")
coef(mmmm) %>% plogis
# managrangeland managlivestock exclusion
# 0.1617647                0.2535211
