# Fit survival model with parameters at annual scale, so the likelihood is
# p(survival_annual) ^ 15.
# Compared to fitting at the 15-years scale, only the intercepts and sigma
# change, but not slopes.

# Packages ----------------------------------------------------------------

library(tidyverse)
library(viridis)
library(readxl)
library(posterior)
library(rstan)
library(DHARMa)
library(tidybayes)
library(terra)         # check spatial correlation
library(logitnorm)
library(foreach)       # parallelization of moments_logit_normal
library(doMC)
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

# vectorized logit_norm mean
logit_norm_vect <- function(mu, sigma) {
  if(length(sigma) == 1) sigma <- rep(sigma, length(mu))
  pmean <- numeric(length(mu))

  for(i in 1:length(mu)) {
    pmean[i] <- momentsLogitnorm(mu[i], sigma[i])["mean"]
  }
  return(pmean)
}

# compute the mean of the logit-normal distribution in parallel.
# Intended to be used with a matrix of logit-means and matrix of logit-sd,
# of the same size. the sds can also be a vector with length == ncol(mu).
# This is thought for columns to be posterior samples.
logit_norm_mean <- function(mu, sigma, cores = 15) {

  registerDoMC(cores = cores)

  sigma_mat <- !is.null(dim(sigma))

  # turn columns into list elements
  if(sigma_mat) {
    arg_list <- lapply(1:ncol(mu), function(j) {
      r <- cbind(mu[, j], sigma[, j])
      colnames(r) <- c("mu", "sigma")
      return(r)
    })
  } else {
    arg_list <- lapply(1:ncol(mu), function(j) {
      r <- cbind(mu[, j], sigma[j])
      colnames(r) <- c("mu", "sigma")
      return(r)
    })
  }

  # compute moments in parallel
  means_list <- foreach(cc = arg_list) %dopar% {
    logit_norm_vect(cc[, "mu"], cc[, "sigma"])
  }

  return(do.call("cbind", means_list))
}

# compute annual survival, from the probability (p) of surviving in
# period years, assuming independent annual survival.
surv_annual <- function(p, period = 15) p ^ (1 / period)

# function to summarize data by group, taking into account that groups should
# have at least O observations.
average <- function(data, response, var, O = 15) {
  # data = d; var = "height"; response = "surv"; O = 15

  ds <- split(data, data$manag_fire)

  rr <- do.call("rbind", lapply(ds, function(loc) {
    # loc = ds[[2]]
    ng <- floor(nrow(loc) / O)
    enes <- rep(O, ng)
    extra <- nrow(loc) - ng * O
    enes[ng] <- enes[ng] + extra
    fac <- rep(letters[1:ng], enes)

    dord <- loc[order(loc[, var, drop = T]), ]
    dord$fac <- fac
    dord$response <- dord[, response, drop = T]
    dord$var <- dord[, var, drop = T]
    dagg <- aggregate(cbind(response, var) ~ manag + fire + fac, dord, mean)
    rownames(dagg) <- NULL
    return(dagg)
  }))
  rownames(rr) <- NULL
  return(rr)
}


# Import data -------------------------------------------------------------

data_name <- "Plots_remedicion 2024.03.10.xls"

# plot-level data
dplot <- read_excel(file.path("..", data_name), sheet = "plotsall_use")
nrow(dplot) # initial number of plots = 139

# tree-level data
dtree <- read_excel(file.path("..", data_name), sheet = "taball_use")

# subset plots to use, with reliable data
dplot <- dplot[is.na(dplot$`delete?`), ] # with no comment in "delete?"
nrow(dplot) # used number of plots = 139

dtree <- dtree[dtree$plot %in% dplot$plot, ]

sum(!is.na(as.numeric(dtree$sur2018))) # 841 without NA (identified and not)
nrow(dtree) # 1229 original trees
dtree$plot[is.na(as.numeric(dtree$sur2018))] %>% unique %>% length
# 50 plots where there were trees without tags

length(unique(dtree$plot))
ddd <- aggregate(hei2003 ~ plot, dtree, length)
ddd$na_count <- numeric(nrow(ddd))
for(i in 1:nrow(ddd)) {
  ddd$na_count[i] <- sum(is.na(as.numeric(dtree$sur2018[dtree$plot == ddd$plot[i]])))
}

which(ddd$na_count >= ddd$hei2003)
ddd$plot[101] # tiene todos NA, pero se sabe que tenía 2 muertos y 2 vivos

# waypoints for spatial correlation check
vplot <- vect(file.path("..", "plots_waypoints.shp"))
# plot(vplot)

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
                      labels = c("Rangeland", "Livestock exclusion"))

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

# save means to plot later
dplot$surv <- (dplot$nalive / dplot$ntrees) ^ (1 / 15) * 100
manag_means <- aggregate(surv ~ manag, dplot, mean)
# write.csv(manag_means, "files/table_means_growth_manag.csv", row.names = F)

table(dplot$manag)

plot(surv ~ ntrees, dplot)
mean(dplot$surv)

# Identify complete and incomplete plots -----------------------------------

dn <- data.frame(plot = dtree$plot, surv = dtree$surv)
dtagg <- aggregate(surv ~ plot, dn, mean, drop = F) # do not remove plot 417!!
names(dtagg)[2] <- "meansurv"

dtagg$na <- numeric(nrow(dtagg))
dtagg$length <- numeric(nrow(dtagg))
for(i in 1:nrow(dtagg)) {
  dtagg$na[i] <- sum(is.na(dtree$surv[dtree$plot == dtagg$plot[i]]))
  dtagg$length[i] <- length(dtree$surv[dtree$plot == dtagg$plot[i]])
}
sum(dtagg$na > 0) # 50 plots with NA in survival

# all(dtagg$length >= dtagg$na) # OK
dtagg <- left_join(dtagg, dplot[, c("plot", "ntrees", "nalive", "ndead")], by = "plot")

dtagg[dtagg$ntrees < dtagg$length, ]
dtagg[dtagg$ntrees < dtagg$length, ] %>% nrow - 1 # 46 plots with unidentified trees

dtagg$miss <- dtagg$length - dtagg$ntrees
sum(dtagg$miss) # 349 missing trees.
sum(dtagg$miss) / sum(dtagg$length) # 28.5 % missing

# plot(miss ~ length, dtagg); abline(0, 1) # ok, always miss < total
# plot(miss ~ na, dtagg); abline(0, 1) # ok, always na >= miss

# evaluate how many counted trees have no matching id.
# counted are ntrees, and identified are length - na
dtagg$identified <- dtagg$length - dtagg$na
sum(dtagg$identified) # 841
dtagg$unidentified <- dtagg$ntrees - dtagg$identified
sum(dtagg$unidentified) # 39

nrow(dtree[!is.na(dtree$surv), ]) # 841 sin NA en survival
# pero tienen NA en growth y dieback

# subset plot with unidentified trees
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

# plot(ndead_unid ~ nalive_unid, dsub); abline(0, 1) # most unid are dead
# in cases where only dead or only alive trees are unidentified,
# the problem is much smaller.
sum(dsub$unidentified) # 39 trees
sum(dsub$unidentified) / sum(dtagg$ntrees) * 100
# 4.43 % of non-missing trees are unidentified

# list with ones and zeroes for unidentified trees
bin_list <- lapply(1:nrow(dsub), function(r) {
  c(rep(0, dsub$ndead_unid[r]), rep(1, dsub$nalive_unid[r]))
})

trees_available_list <- lapply(1:nrow(dsub), function(r) {
  dtree_local <- dtree[dtree$plot == dsub$plot[r], ]
  x <- dtree_local$tree_id[is.na(dtree_local$surv)]
  return(x)
  # return(unique(x))
})

# make non-redundant combinations of tree ids for each plot
trees_available_combs <- lapply(1:nrow(dsub), function(i) {
  # i = 7
  y <- bin_list[[i]]
  arbol <- trees_available_list[[i]]

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

dsub$combinations_unique <- sapply(trees_available_combs, ncol)
# In these plots we will compute the marginal likelihood with respect to all
# possible id combinations.

# https://mbjoseph.github.io/posts/2020-04-28-a-step-by-step-guide-to-marginalizing-over-discrete-parameters-for-ecologists-using-stan/

# Whithin each plot, for each possible id combination, compute the sum of
# binomial_logit_lpmf. Then, apply log_sum_exp to that vector of likelihoods,
# to get the likelihood of the plot marginal to which were the unidentified
# trees.

# dtagg has data at the parcel level indicating how many trees there are,
# how many na, and how many identified.
# for plots where unidentified == 0, we just remove surv == NA.
# For plots with unidentified > 0, we separate the NA and the obverved.
dtagg %>% dim()
sum(dtagg$unidentified > 0) # 13 plots == nrow(dsub)
sum(dtagg$unidentified == 0) # 126 normal plots, where NA will be removed.


# get plots with no unidentified trees
data_id <- dtree[dtree$plot %in% dtagg$plot[dtagg$unidentified == 0], ]
# remove NAs
data_id <- data_id[complete.cases(data_id[, c("height", "elev", "pforest",
                                              "manag", "fire", "surv")]), ]
nrow(data_id) # 729

# get plots with unidentified trees
data_unid <- dtree[dtree$plot %in% dsub$plot, ]

# non-na trees
data_unid_obs <- data_unid[!is.na(data_unid$surv), ]
nrow(data_unid_obs) # 112
sum(1-data_unid_obs$surv) # 18 unid dead
sum(data_unid_obs$surv)   # 94 unid alive

# na trees
data_unid_na <- data_unid[is.na(data_unid$surv), ]
nrow(data_unid_na) # 39 were recorded but unidentified, and there are 127 options.

# merge non-na trees
data_obs <- rbind(data_id, data_unid_obs)
nrow(data_obs) # 841

anyNA(data_obs[, -which(names(data_obs) %in% c("surv", "age", "perim", "age_z"))])
# View(data_obs) # OK

anyNA(data_unid_na[, -which(names(data_unid_na) %in% c("surv", "age", "age_z", "perim"))])
# View(data_unid_na)
# OK, but there is one NA in perim at unid_na


ddddd <- rbind(data_unid_obs, data_unid_na)
ddddd$plot %>% unique %>% length
# Note about age data -----------------------------------------------------

# There is one NA at perim2003. If age is used, this code might break.

# Data for Stan -----------------------------------------------------------

# make three design matrices for each dataset:
# fixed effects (X), plot (Zp), basin (Zb).
# Suffix _unid will indicate the matrix for unidentified trees
fixed_formula <- formula(~ manag + fire +
                           height_z + elev_z + pforest_z)
plot_formula <- formula(~ - 1 + plot)

X <- model.matrix(fixed_formula, data = data_obs)
X_unid <- model.matrix(fixed_formula, data = data_unid_na)

Zp <- model.matrix(plot_formula, data = data_obs)
Zp_unid <- model.matrix(plot_formula, data = data_unid_na)

# data for marginalization of unidentified trees

# this is the response variable, which has to be assigned to different combinations
# of unidentified trees
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

# for each plot with unidentified trees, we have the available trees at
# trees_available_combs. And to use it in Stan, we need to map the corresponding
# tree ids to rows in the data_unid_na. But Stan does not accept lists, so
# we collapse all the hypothetical lists into a matrix.

max_combs <- max(sapply(trees_available_combs, ncol)) # ncol to use in each matrix
combs_rows_matrix <- do.call("rbind", lapply(1:p_unid, function(p) {
  # p = 1
  # print(p)
  cc <- trees_available_combs[[p]]
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

stan_data <- list(
  N = nrow(data_obs),
  N_unid_rows = nrow(data_unid_na),
  N_unid_trees = length(y_unid),

  K = ncol(X),
  P = ncol(Zp),
  P_unid = length(trees_available_combs),

  y = data_obs$surv,
  y_unid = y_unid,

  X = X, X_unid = X_unid,
  Zp = Zp, Zp_unid = Zp_unid,

  y_unid_start = y_unid_start,
  y_unid_end = y_unid_end,
  y_unid_length = y_unid_length,

  max_combs = max(sapply(trees_available_combs, ncol)),
  combs_rows_matrix = combs_rows_matrix,
  combs_n = sapply(trees_available_combs, ncol),

  unid = 1, # use unidentified trees in likelihood.

  prior_intercept_sd = 5,
  prior_b_sd = 3,
  prior_sigma_sd = 3
)


# Model fit ---------------------------------------------------------------

stan_code <- stan_model("survival.stan", verbose = T)

m1 <- sampling(stan_code, data = stan_data, seed = 1234, refresh = 10,
               # cores = 1, chains = 1, iter = 10, ## test
               cores = 6, chains = 6, iter = 2000,
               control = list(adapt_delta = 0.9))
saveRDS(m1, file.path("files", "survival_model_samples_annual_elevLinear.rds"))
# 690 / 60 = 11.5  min
m1 <- readRDS(file.path("files", "survival_model_samples_annual_elevLinear.rds"))
# m1 <- readRDS(file.path("files", "survival_model_samples_annual.rds"))
sm1 <- summary(m1)[[1]]
max(sm1[, "Rhat"]); min(sm1[, "n_eff"]) # nice
# [1] 1.001151
# [1] 1242.651

ids <- c(grep("b\\[", rownames(sm1)),
         which("sigma_plot" == rownames(sm1)))

# compute proportion of the posterior above zero
bbbb <- as.matrix(m1)[, ids]
# colnames(bbbb)
pgt0 <- apply(bbbb, 2, function(x) sum(x > 0) / length(x))

summ_survival <- sm1[ids, ]
summ_survival <- cbind(variable = c(colnames(X), "sigma_plot"),
                       as.data.frame(sm1[ids, ]),
                       prob_gt0 = pgt0)

write.csv(summ_survival,
          file.path("files", "table_summary_survival.csv"),
          row.names = F)

# Residuals analyses ------------------------------------------------------

# extract parameters
bhat <- as.matrix(m1, "b") %>% t
sigma_plot <- as.matrix(m1, "sigma_plot") %>% as.numeric
nsim <- length(sigma_plot)

# simulate data for new basins and plots, to compute DHARMa residuals.

# simulation for for identified trees
e_plot_sim <- matrix(rnorm(nrow(X) * nsim, sd = sigma_plot),
                     nrow(X), nsim)
eta <- X %*% bhat + e_plot_sim
psim <- plogis(eta) ^ 15 # 15 years survival
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
eta <- X_unid %*% bhat + e_plot_sim
psim_options <- plogis(eta) ^ 15
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
                             prob = psim_options[rows_pred, i], log = TRUE) %>% sum %>% exp
    }

    # sample a combination of trees identities based on the likelihood
    chosen_comb <- sample(1:nc, size = 1, prob = like_combs)

    # identify to which rows in the data of unidentified trees (N_unid_rows = 123)
    # correspond the chosen options
    rows_pred <- combs_rows_matrix[rows_response, chosen_comb]

    # get probabilities and simulate
    psim_chosen <- psim_options[rows_pred, i]
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

res_data <- list(res = res, data_both = data_both)
saveRDS(res_data, file.path("files", "survival_model_residuals_and_data_elevLinear.rds"))
#
res_data <- readRDS(file.path("files", "survival_model_residuals_and_data_elevLinear.rds"))
data_both <- res_data$data_both
# plot model residuals in a fancy way
# manag, fire
# elev, pforest, height


# Spatial correlation check ------------------------------------------------

vplot <- vplot[vplot$Name %in% levels(dplot$plot), ]
vplot <- vplot[order(vplot$Name), ]
# all(levels(dplot$plot) == vplot$Name) # OK

# compute spatial distance across observations
spat_dist <- distance(vplot, symmetrical = TRUE) %>% as.matrix
hist(spat_dist[lower.tri(spat_dist)] / 1000, breaks = 30,
     main = "Pairwise distance distribution (km)")

# residuals diff
res_agg <- aggregate(res ~ plot, data_both, mean)
res_dist <- as.matrix(dist(res_agg$res))

dist_df <- data.frame(res_diff = res_dist[lower.tri(res_dist)],
                      distance = spat_dist[lower.tri(spat_dist)] / 1000)

ggplot(dist_df, aes(distance, res_diff)) +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 15),
              method.args = list(family = mgcv::betar()),
              color = "red")


# Residuals plots ---------------------------------------------------------

res_elev <-
ggplot(data_both, aes(x = elev, y = res)) +
  geom_hline(yintercept = c(0.25, 0.5, 0.75), linetype = "dashed",
             linewidth = 0.3, color = viridis(1)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 10, bs = "cr"),
              method.args = list(family = mgcv::betar()),
              color = viridis(1, begin = 0.4),
              fill = viridis(1, begin = 0.4)) +
  scale_y_continuous(expand = c(0.01, 0.01), limits = c(0, 1)) +
  xlab("Elevation (m a.s.l.)") +
  ylab("DHARMa residuals") +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        plot.title = element_text(hjust = 0),
        panel.grid.minor = element_blank()) +
  ggtitle("A")

res_pforest <-
ggplot(data_both, aes(x = pforest, y = res)) +
  geom_hline(yintercept = c(0.25, 0.5, 0.75), linetype = "dashed",
             linewidth = 0.3, color = viridis(1)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 10, bs = "cr"),
              method.args = list(family = mgcv::betar()),
              color = viridis(1, begin = 0.4),
              fill = viridis(1, begin = 0.4)) +
  scale_y_continuous(expand = c(0.01, 0.01), limits = c(0, 1)) +
  xlab("Forest cover (%)") +
  ylab("DHARMa residuals") +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        plot.title = element_text(hjust = 0),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(color = "white")) +
  ggtitle("B")

res_height <-
ggplot(data_both, aes(x = height, y = res)) +
  geom_hline(yintercept = c(0.25, 0.5, 0.75), linetype = "dashed",
             linewidth = 0.3, color = viridis(1)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 10, bs = "cr"),
              method.args = list(family = mgcv::betar()),
              color = viridis(1, begin = 0.4),
              fill = viridis(1, begin = 0.4)) +
  scale_y_continuous(expand = c(0.01, 0.01), limits = c(0, 1)) +
  xlab("Tree height (cm)") +
  ylab("DHARMa residuals") +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        plot.title = element_text(hjust = 0),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(color = "white")) +
  ggtitle("C")


data_both$managb <- factor(data_both$manag, levels = levels(data_both$manag),
                           labels = c("Rangeland", "Livestock\nexcusion"))
res_manag <-
ggplot(data_both, aes(x = managb, y = res)) +
  geom_hline(yintercept = c(0.25, 0.5, 0.75), linetype = "dashed",
             linewidth = 0.3, color = viridis(1)) +
  geom_boxplot(width = 0.5, fill = NA) +
  scale_y_continuous(expand = c(0.01, 0.01), limits = c(0, 1)) +
  xlab("Management") +
  ylab("DHARMa residuals") +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        plot.title = element_text(hjust = 0),
        panel.grid.minor = element_blank()) +
  ggtitle("D")

res_fire <-
ggplot(data_both, aes(x = fire, y = res)) +
  geom_hline(yintercept = c(0.25, 0.5, 0.75), linetype = "dashed",
             linewidth = 0.3, color = viridis(1)) +
  geom_boxplot(width = 0.5, fill = NA) +
  scale_y_continuous(expand = c(0.01, 0.01), limits = c(0, 1)) +
  xlab("Fire occurrence") +
  ylab("DHARMa residuals") +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        plot.title = element_text(hjust = 0),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(color = "white")) +
  ggtitle("E")

res_spat <-
ggplot(dist_df, aes(distance, res_diff)) +
  geom_point(alpha = 0.15) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 50),
              method.args = list(family = mgcv::betar()),
              color = viridis(1, begin = 0.4),
              fill = viridis(1, begin = 0.4)) +
  scale_y_continuous(expand = c(0.01, 0.01), limits = c(0, 1)) +
  xlab("Distance (km)") +
  ylab("Residuals difference") +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        plot.title = element_text(hjust = 0),
        panel.grid.minor = element_blank()) +
  ggtitle("F")


res_surv <- egg::ggarrange(res_elev, res_pforest, res_height,
                           res_manag, res_fire, res_spat,
                           ncol = 3)
ggsave("figures/survivial_residuals.png", plot = res_surv,
       width = 17, height = 11, units = "cm")

# Fire ~ management -------------------------------------------------------

# management affects survival mostly mediated by fire, which is more frequent
# under livestock exclusion.
# Here we fit a model of fire probability at the plot level as a function of
# management, to estimate fire probabilities and include their uncertainty later,
# when the fire-weighted-average of predictions will be computed.

# S <- as.numeric(table(dplot$manag))
# y <- as.numeric(table(dplot$manag[dplot$fire == "Burned plots"]))
#
# fire_data <- list(N = length(y), S = S, y = y)
# fire_model <- stan_model("fire.stan", verbose = T)
#
# f1 <- sampling(fire_model, data = fire_data,
#                cores = 6, chains = 6, iter = 2000, warmup = 1000)
#
# pfire <- as.matrix(f1, "p")
# colnames(pfire) <- levels(dplot$manag)
# saveRDS(pfire, file.path("files", "fire_prob_samples.rds"))
pfire <- readRDS(file.path("files", "fire_prob_samples.rds"))


# probability of higher fire prob at livestock exlusion
firediff_p <- sum(pfire[, 2] > pfire[, 1]) / nrow(pfire) # 0.9175
firediff <- pfire[, 2] - pfire[, 1]

fff <- cbind(pfire, diff = firediff)

format(round(apply(pfire, 2, mean_ci) %>% t, 3), nsmall = 3)
#                      mean    lower   upper
# Rangeland           "0.161" "0.085" "0.256"
# Livestock exclusion "0.254" "0.161" "0.356"

fire_table <- apply(fff, 2, mean_ci) %>% t
fire_table <- rbind(fire_table, matrix(c(firediff_p, NA, NA), 1))
rownames(fire_table)[4] <- "prob"
write.csv(fire_table, file.path("files", "table_summary_fire.csv"))

# Predictions survival ~ predictors ----------------------------------------

manag_lev <- levels(dplot$manag)
fire_lev <- levels(dplot$fire)

height_lim <- (1200 - pred_means["height"]) / pred_sds["height"]

pdata <- rbind(
  # elevation
  expand.grid(
    manag = factor(manag_lev, manag_lev),
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
    height_z = seq(min(dtree$height_z), height_lim, length.out = 150),
    varying_var = "Tree height (cm)"
  )
)

# tidy factors

# add variables at original scale
pdata$varying_val <- 0

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

# extract parameters
bhat <- as.matrix(m1, "b") %>% t
sigma_plot <- as.matrix(m1, "sigma_plot") %>% as.numeric

eta <- Xpred %*% bhat
phat <- logit_norm_mean(eta, sigma_plot, cores = 8)

# summarize posterior
psumm <- apply(phat, 1, mean_ci) %>% t %>% as.data.frame %>% "*"(100)
predictions <- cbind(pdata, psumm)

vv <- unique(pdata$varying_var)
vv_names <- paste(c("(1)", "(2)", "(3)"), unique(pdata$varying_var))
predictions$managlet <- factor(predictions$manag,
                               levels = levels(predictions$manag),
                               labels = c("A. Rangeland",
                                          "B. Livestock exclusion"))

saveRDS(predictions, file.path("files", "survival_predictions.rds"))
predictions <- readRDS(file.path("files", "survival_predictions.rds"))

plist <- vector("list", 3)

breaks <- list(
  seq(1200, 2400, by = 300),
  seq(0, 100, by = 25),
  seq(200, 1200, by = 200)
)

# summarized data to plot alongside predictions
d <- res_data$data_both
d$manag_fire <- factor(
  paste(d$manag, d$fire, sep = "_"),
  levels = unique(paste(d$manag, d$fire, sep = "_"))
)

vars <- c("elev", "pforest", "height")
dlist <- vector("list", 3)
for(i in 1:3) {
  # i = 1
  dd <- average(d, "surv", vars[i], O = 15)
  dd$response <- dd$response ^ (1/15)
  dd$managlet <- factor(dd$manag,
                        levels = levels(dd$manag),
                        labels = c("A. Rangeland",
                                   "B. Livestock exclusion"))
  dlist[[i]] <- dd
}
lapply(dlist, function(x) table(x$manag, x$fire))

for(v in 1:3) {
  # v = 1
  ddd <- predictions[predictions$varying_var == vv[v], ]

  # add poitns
  # dp <- dplot
  # dp$x_class <- cut(dp[, ])

  pp <-
    ggplot(ddd, aes(x = varying_val, y = mean, ymin = lower, ymax = upper,
                  color = fire, fill = fire)) +
    geom_ribbon(color = NA, alpha = 0.25) +
    geom_line() +
    geom_point(data = dlist[[v]],
               mapping = aes(x = var, y = response * 100,
                             color = fire, fill = fire),
               inherit.aes = F, alpha = 0.7, size = 0.8) +
    # geom_rug(mapping = aes(x = elev, color = fire, fill = fire),
    #          data = d, inherit.aes = F) +
    scale_color_viridis(discrete = T, option = "B", end = 0.35) +
    scale_fill_viridis(discrete = T, option = "B", end = 0.35) +
    facet_wrap(~ managlet, ncol = 2) +
    scale_y_continuous(limits = c(65, 100), expand = c(0.01, 0.01)) +
    scale_x_continuous(breaks = breaks[[v]],
                       limits = range(c(breaks[[v]], d$varying_val))) +
    xlab(vv_names[v]) +
    ylab("Annual survival probability (%)") +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(color = "white", fill = "white"),
          panel.spacing = unit(3, "mm"),
          axis.title.y = element_blank(),
          legend.title = element_blank(),
          legend.position = "none",
          axis.text = element_text(size = 8),
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
ggsave("figures/survival_predictions.png", plot = surv_plot,
       width = 16, height = 15, units = "cm")

# Predictions survival ~ management --------------------------------------

# average the management predictions with respect to fire, but considering that
# fire varies between management levels.
pdata_manag <- expand.grid(
  manag = factor(manag_lev, manag_lev),
  fire = factor(fire_lev, fire_lev),
  elev_z = 0,
  pforest_z = 0,
  height_z = 0
)

Xpred_manag <- model.matrix(fixed_formula, pdata_manag)
eta_manag <- Xpred_manag %*% bhat
phat_manag <- logit_norm_mean(eta_manag, sigma_plot, cores = 14)

# order pfire in the same way as pdata_manag
pfire_rep <- rbind(1 - t(pfire), t(pfire))
prods <- phat_manag * pfire_rep

means_manag <- aggregate(prods ~ manag, pdata_manag, sum)
dim(means_manag)
rownames(means_manag) <- means_manag[, 1]
means_manag_hat <- as.matrix(means_manag[, -1])
means_manag_summ <- apply(means_manag_hat, 1, mean_ci) %>% t

format(round(means_manag_summ * 100, 2), nsmall = 2)
#                       mean    lower   upper
# Rangeland           "99.18" "98.45" "99.58"
# Livestock exclusion "98.88" "98.01" "99.39"

diffmeans <- means_manag_hat[2, ] - means_manag_hat[1, ]
plot(density(diffmeans))
format(round(mean_ci(diffmeans) * 100, 2), nsmall = 2)

# diferencia en la media condicional entre manejo (exclusion - rangeland)
#   mean   lower   upper
# "-0.30" "-1.16" " 0.48"
probb <- sum(means_manag_hat[2, ] < means_manag_hat[1, ]) / ncol(means_manag_hat)
# prob(exclusion < rangeland) = 0.8011667

# save table
eeexp <- rbind(means_manag_summ, mean_ci(diffmeans),
               matrix(c(probb, NA, NA), 1))
rownames(eeexp)[3:4] <- c("diff", "prob")

write.csv(eeexp, file.path("files", "table_summary_survival_manag_marginal.csv"))

