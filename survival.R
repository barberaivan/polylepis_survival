# Packages ----------------------------------------------------------------

library(tidyverse)
library(viridis)
library(readxl)
library(posterior)
library(rstan)
library(tidybayes)
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

dpar <- read_excel(file.path("..", "Plots_remedicion 2024.03.06.xls"),
                   sheet = "parcelas_usar")
dpar <- dpar[-nrow(dpar), ] # remaining row
dtree <- read_excel(file.path("..", "Plots_remedicion 2024.03.06.xls"),
                    sheet = "tab")

dtree_neat <- read_excel(file.path("..", "Plots_remedicion 2024.03.06.xls"),
                         sheet = "tab_usar")


# Tidy data ---------------------------------------------------------------

dtree <- dtree[, c("Nparcela...19", "arbol", "Vivo1-M0_2018",
                   "alt_2003", "perim. (cm)_2003", "% roca copa")]
names(dtree) <- c("plot", "tree_id", "surv", "height", "perim", "rock")

dtree$plot <- factor(dtree$plot)
dpar$plot <- factor(dpar$plot)

# add variables
dtree$age <- tree_age(dtree$perim, dtree$rock)

dpar_sub <- dpar[, c("plot", "manag2", "fire01", "elev", "pforest", "basin")]
dtree <- left_join(dtree, dpar_sub, by = "plot")

# tidy factors
dtree$manag <- factor(dtree$manag2, levels = c("a", "b"),
                      labels = c("rangeland", "livestock exclusion"))
dtree$fire <- factor(dtree$fire01, levels = c("0", "1"),
                     labels = c("unburned", "burned"))
dtree$plot <- factor(dtree$plot)
dtree$basin <- factor(dtree$basin)
dtree$tree_id <- paste(dtree$plot, dtree$tree_id, sep = "_")
# resolve duplicated ids
dtree$tree_id[duplicated(dtree$tree_id)] <-
  paste(dtree$tree_id[duplicated(dtree$tree_id)], "b", sep = "")
dtree$tree_id <- factor(dtree$tree_id)

# standardize continuous predictors
dtree$height_z <- as.numeric(scale(dtree$height))
dtree$age_z <- as.numeric(scale(dtree$age))
dtree$elev_z <- as.numeric(scale(dtree$elev))
dtree$pforest_z <- as.numeric(scale(dtree$pforest))

# remove rows with NA in any useful predictor
out <- is.na(dtree$height) | #is.na(dtree$age) |
       is.na(dtree$plot) | is.na(dtree$elev) | is.na(dtree$pforest)
dtree <- dtree[!out, ] # 1230 available cases

# Identify complete and incomplete parcels --------------------------------

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
dtagg <- left_join(dtagg, dpar[, c("plot", "ntrees", "nalive", "ndead")], by = "plot")
dtagg <- dtagg[complete.cases(dtagg), ]

dtagg[dtagg$ntrees < dtagg$length, ]
dtagg[dtagg$ntrees < dtagg$length, ] %>% nrow - 1 # 49 parcels

dtagg$miss <- dtagg$length - dtagg$ntrees
sum(dtagg$miss) # 350 missing trees.
sum(dtagg$miss) / sum(dtagg$length) # 28.5 % missing

plot(miss ~ length, dtagg); abline(0, 1) # ok, always miss < total
plot(miss ~ na, dtagg); abline(0, 1) # ok, always na >= miss

# evaluate how many counted trees have no matching id.
# counted are ntrees, and identified are length - na
dtagg$identified <- dtagg$length - dtagg$na
dtagg$unidentified <- dtagg$ntrees - dtagg$identified

# see those parcels
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
# En cada una de estas parcelas, habría que marginalizar la likelihood con
# respecto a esas combinaciones posibles de predictoras.

# https://mbjoseph.github.io/posts/2020-04-28-a-step-by-step-guide-to-marginalizing-over-discrete-parameters-for-ecologists-using-stan/

# Para cada combinación posible, calcular la suma de la binomial_logit_lpmf.
# Luego, aplicar log_sum_exp() a ese vector que contiene tantas
# sum(binomial_lpmf) como combinaciones. Esa es la loglik de ese set de
# observaciones, marginal a cuáles hayan sido los ids.


# dtagg has data at the parcel level indicating how many trees there are,
# how many na, and how many identified.
# for parcels where unidentified == 0, we jus remove surv == NA.
# For parcels with unidentified > 0, we separate the NA and the obverved.

dtagg %>% dim()
sum(dtagg$unidentified > 0) # 12 parcels == nrow(dsub)
sum(dtagg$unidentified == 0) # 126 normal parcels, where na will be removed.

# get parcels with no unidentified trees
data_id <- dtree[dtree$plot %in% dtagg$plot[dtagg$unidentified == 0], ]
# remove NAs
data_id <- data_id[complete.cases(data_id), ]

# get parcels with unidentified trees
data_unid <- dtree[dtree$plot %in% dsub$plot, ]
# non-na trees
data_unid_obs <- data_unid[!is.na(data_unid$surv), ]
# na trees
data_unid_na <- data_unid[is.na(data_unid$surv), ]
nrow(data_unid_na) # 35 were recorded but unidentified, but there are 123 options.

# merge non-na trees
data_obs <- rbind(data_id, data_unid_obs)
nrow(data_obs) # 838

anyNA(data_obs)
anyNA(data_unid_na[, -which(names(data_unid_na) == "surv")])
# OK, but there is one NA in perim at unid_na

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

# va al palo!

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
  prior_sigma_sd = 3
)


# Model fit ---------------------------------------------------------------

stan_code <- stan_model("survival_model.stan", verbose = T)

# consider unidentified trees
m1 <- sampling(stan_code, data = stan_data, seed = 1234, refresh = 10,
               # cores = 1, chains = 1, iter = 10,
               cores = 6, chains = 6, iter = 2000,
               control = list(adapt_delta = 0.9))
# 604.862 / 60 = 10.08 min
sm1 <- summary(m1)[[1]]
max(sm1[, "Rhat"]); min(sm1[, "n_eff"]) # nice

# ingnore unidentified trees
stan_data2 <- stan_data
stan_data2$unid <- 0
m2 <- sampling(stan_code, data = stan_data2, seed = 1234, refresh = 100,
               # cores = 1, chains = 1, iter = 10,
               cores = 6, chains = 6, iter = 2000,
               control = list(adapt_delta = 0.9, max_treedepth = 20))
sm2 <- summary(m2)[[1]] # anduvo!
max(sm2[, "Rhat"]); min(sm2[, "n_eff"]) # nice

# compare posteriors
bnames <- colnames(X)
par_names <- c(bnames, "sigma_plot", "sigma_basin")
bm1 <- as.matrix(m1, par = c("b", "sigma_plot", "sigma_basin"))
bm2 <- as.matrix(m2, par = c("b", "sigma_plot", "sigma_basin"))
colnames(bm1) <- colnames(bm2) <- par_names

# plots to compare

par(mfrow = c(3, 3))
for(j in 1:length(par_names)) {
  # j = 1

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
# parameters change quite much.



# Predictions -------------------------------------------------------------

manag_lev <- c("Rangeland", "Livestock exclusion")
fire_lev <- c("Unburned", "Burned")

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

rows_replace <- pdata$varying_var == "Elevation (m a.s.l.)"
pdata$varying_val[rows_replace] <-
  pdata$elev_z[rows_replace] * sd(dtree$elev) + mean(dtree$elev)

rows_replace <- pdata$varying_var == "Forest cover (%)"
pdata$varying_val[rows_replace] <-
  pdata$pforest_z[rows_replace] * sd(dtree$pforest) + mean(dtree$pforest)

rows_replace <- pdata$varying_var == "Tree height (cm)"
pdata$varying_val[rows_replace] <-
  pdata$height_z[rows_replace] * sd(dtree$height) + mean(dtree$height)

Xpred <- model.matrix(fixed_formula, pdata)

# choose which model to use!
bhat <- as.matrix(m1, "b") %>% t
bhat <- as.matrix(m2, "b") %>% t

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
  # dp <- dpar
  # dp$x_class <- cut(dp[, ])

  pp <-
  ggplot(d, aes(x = varying_val, y = mean, ymin = lower, ymax = upper,
                          color = fire, fill = fire)) +
    geom_ribbon(color = NA, alpha = 0.25) +
    geom_line() +
    scale_color_viridis(discrete = T, option = "B", end = 0.5) +
    scale_fill_viridis(discrete = T, option = "B", end = 0.5) +
    facet_wrap(~ manag, ncol = 2) +
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
ggsave("figures/survival_prediction_all_data.png", plot = surv_plot,
       width = 14, height = 15, units = "cm")
ggsave("figures/survival_prediction_found_tags.png", plot = surv_plot,
       width = 14, height = 15, units = "cm")

summary(dtree)

# Fire and management -----------------------------------------------------

table(dpar$fire01, dpar$manag2)
mmmm <- glm(fire01~manag2, data = dpar, family = "binomial")
summary(mmmm)
mmmm <- glm(fire01~manag2-1, data = dpar, family = "binomial")
coef(mmmm) %>% plogis

