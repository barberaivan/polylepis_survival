# Fitting a tree growth model at the annual scale.
# Inherited from survival.R, which might explain some unused code.

# Packages ----------------------------------------------------------------

library(tidyverse)
library(viridis)
library(ggh4x)
library(deeptime)
library(grid)
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

# function to summarize data by group, taking into account that groups should
# have at least O observations.
average <- function(data, response, var, O = 15) {
  # data = dtree; var = "height"; response = "growth"; O = 15

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

data_name <- "Plots_remedicion 2024.08.14.xls"

# plot-level data
dplot <- read_excel(file.path("..", data_name), sheet = "plotsall_use")

# tree-level data
dtree <- read_excel(file.path("..", data_name), sheet = "taball_use")

dtree$plot <- dtree$nplot
dtree <- dtree[dtree$plot %in% dplot$plot, ]

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
dtree <- rename(dtree, surv = "sur2018", growth = "groyr", dieback = "dieb2018",
                height = "hei2003", perim = "peri2003", tree_id = "arbol")

dtree <- dtree[, c("plot","growth", "dieback", "height", "perim",
                   "tree_id", "prock")]
dtree$growth <- as.numeric(dtree$growth)
dtree$dieback <- as.numeric(dtree$dieback)
dtree$height <- as.numeric(dtree$height)

# make factors
dtree$plot <- factor(dtree$plot)
dplot$plot <- factor(dplot$plot)
dplot$basin <- factor(dplot$basin)
dplot$fire <- factor(as.character(dplot$fire), levels = c("0", "1"),
                     labels = c("Unburned plots", "Burned plots"))
dplot$manag <- factor(dplot$manag, levels = c("a", "b"),
                      labels = c("Ranching", "Conservation"))
dplot$manag_fire <- paste(dplot$manag, dplot$fire, sep = "_")
dplot$manag_fire <- factor(dplot$manag_fire,
                           levels = c(
                             "Ranching_Unburned plots", "Ranching_Burned plots",
                             "Conservation_Unburned plots", "Conservation_Burned plots"
                           ))
levels_managfire <- levels(dplot$manag_fire)

# add age
dtree$perim <- as.numeric(dtree$perim)
dtree$age <- tree_age(dtree$perim, dtree$prock)

# add plot variables to tree data
dplot_sub <- dplot[, c("plot", "manag", "fire", "manag_fire",
                       "elev", "pforest", "basin")]
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

# note that standardization is performed with the same data as in survival model,
# including more trees.

# remove rows with NA in response or in any useful predictor
out <- is.na(dtree$growth) | is.na(dtree$height) |
       is.na(dtree$plot) | is.na(dtree$elev) | is.na(dtree$pforest)
dtree <- dtree[!out, ]
nrow(dtree) # 662 available cases

nrow(dtree[!is.na(dtree$growth), ]) # 662
nrow(dtree[!is.na(dtree$growth), ]) # 662

# save means to plot later
plot_means2 <- aggregate(growth ~ plot + manag, dtree, mean)
manag_means <- aggregate(growth ~ manag, plot_means2, mean)
write.csv(manag_means, "files/table_means_growth_manag.csv", row.names = F)

# Export processed data table to plot data in predictions
saveRDS(dtree, "files/growth_data_dtree.rds")

# Data for Stan -----------------------------------------------------------

# make three design matrices for each dataset:
# fixed effects (X) and plot (Zp)
fixed_formula <- formula(~ manag * fire +
                           height_z +
                           elev_z +
                           pforest_z + I(pforest_z ^ 2))
fixed_formula_sigma <- formula(~ manag_fire - 1 + height_z)
plot_formula <- formula(~ - 1 + plot)

X <- model.matrix(fixed_formula, data = dtree)
Xs <- model.matrix(fixed_formula_sigma, data = dtree)
Zp <- model.matrix(plot_formula, data = dtree)

# set inits for sigma_intercepts
dsig <- aggregate(growth ~ manag_fire, dtree, sd)
sigma_avg <- mean(dsig$growth) # 6.922806

stan_data <- list(
  N = nrow(dtree),
  K = ncol(X),
  Ks = ncol(Xs),
  P = ncol(Zp),

  y = dtree$growth,
  X = X,
  Xs = Xs,
  Zp = Zp,

  prior_intercept_sd = 100,
  prior_b_sd = 100,
  prior_sigma_sd = 100,

  prior_alphasig_mean = log(sigma_avg),
  prior_alphasig_sd = 1.5,
  prior_betasig_sd = 2,

  prior_nu_a = 2,
  prior_nu_b = 0.1
  # prior for nu suggested in
  # https://statmodeling.stat.columbia.edu/2015/05/17/do-we-have-any-recommendations-for-priors-for-student_ts-degrees-of-freedom-parameter/
)

# Model fit ---------------------------------------------------------------

stan_code <- stan_model("growth.stan", verbose = F)
m1 <- sampling(stan_code, data = stan_data, seed = 1234,
               cores = 6, chains = 6, iter = 2000, refresh = 100,
               # cores = 1, chains = 1, iter = 5,
               control = list(adapt_delta = 0.95))
saveRDS(m1, file.path("files", "growth_model_samples.rds"))
# 62.647 s
m1 <- readRDS(file.path("files", "growth_model_samples.rds"))
sm1 <- summary(m1)[[1]]
max(sm1[, "Rhat"]); min(sm1[, "n_eff"]) # nice
# [1] 1.007946
# [1] 1223.66
# OK

ids <- c(grep("b\\[", rownames(sm1)),
         grep("bs\\[", rownames(sm1)),
         grep("nu", rownames(sm1)),
         which("sigma_plot" == rownames(sm1)))

# compute proportion of the posterior above zero
bbbb <- as.matrix(m1)[, ids]
# colnames(bbbb)
pgt0 <- apply(bbbb, 2, function(x) sum(x > 0) / length(x))

summ <- sm1[ids, ]
summ <- cbind(variable = c(colnames(X),
                           paste(colnames(Xs), "sigma", sep = "_"),
                           "nu",
                           "sigma_plot"),
              as.data.frame(sm1[ids, ]),
              prob_gt0 = pgt0)

write.csv(summ,
          file.path("files", "table_summary_growth.csv"),
          row.names = F)


# Check random effects fit by fire ----------------------------------------

sigma_plot <- as.matrix(m1, "sigma_plot") %>% as.numeric
e_plot <- as.matrix(m1, "e_plot") %>% t
e_point <- rowMeans(e_plot)

plots_burned <- dplot$fire == "Burned plots"
plots_unburned <- dplot$fire == "Unburned plots"

e_plot_sim <- sapply(1:length(e_point), function(i) {
  rnorm(ncol(e_plot), 0, sigma_plot)
}) %>% t


plot(density(e_plot_sim[plots_burned]), ylim = c(0, 0.3), main = "burned")
lines(density(e_plot[plots_burned, ], adjust = 1), col = "red")

plot(density(e_plot_sim[plots_unburned]), ylim = c(0, 0.3), main = "unburned")
lines(density(e_plot[plots_unburned, ], adjust = 1), col = "red")
# Good fit!

# Residuals analyses ------------------------------------------------------

# extract parameters
bhat <- as.matrix(m1, "b") %>% t
bs <- as.matrix(m1, "bs") %>% t
nuhat <- as.matrix(m1, "nu") %>% as.numeric
sigma_plot <- as.matrix(m1, "sigma_plot") %>% as.numeric
e_plot <- as.matrix(m1, "e_plot") %>% t
nsim <- length(sigma_plot)

# # simulate data for new plots, to compute DHARMa residuals.
# muhat <- X %*% bhat + Zp %*% e_plot
# sighat <- exp(Xs %*% bs)
# ysim <- sapply(1:ncol(muhat), function(i) {
#   sigma_total <- sqrt(sighat[, i] ^ 2 + sigma_plot[i] ^ 2)
#   return(rt(nrow(X), df = nuhat[i]) * sigma_total + muhat[, i])
# })

# conditional on plots:
muhat <- X %*% bhat + Zp %*% e_plot
sighat <- exp(Xs %*% bs)
ysim <- sapply(1:ncol(muhat), function(i) {
  rt(nrow(X), df = nuhat[i]) * sighat[, i] + muhat[, i]
})

# residuals
res <- createDHARMa(observedResponse = dtree$growth,
                    simulatedResponse = ysim,
                    fittedPredictedResponse = rowMeans(muhat))

dtree$res <- res$scaledResiduals
dtree$res_norm <- qnorm(dtree$res)
plot(density(dtree$res_norm))     # modelando sigma con height anda muy bien!
curve(dnorm(x), add = T, col = 2)

plot(res, rank = F)
plotResiduals(res, form = dtree$height, rank = F) # quadratic??
plotResiduals(res, form = dtree$elev, rank = F) # OK
plotResiduals(res, form = dtree$pforest, rank = F) # quad?
plotResiduals(res, form = dtree$manag_fire)
plotResiduals(res, form = dtree$basin)

ggplot(dtree, aes(x = res_norm)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_density(adjust = 1.5, fill = "red", alpha = 0.3) +
  geom_rug(outside = F) +
  facet_grid(rows = vars(manag), cols = vars(fire)) +
  ggtitle("t-model residuals (normalized)")

# student-t was needed for the model not to underestimate the mean in the
# ranching-burned condition. Residuals are much better now.

# Check fitted vs obs by group

dtree$fitted <- rowMeans(muhat)
ggplot(dtree, aes(fitted, growth)) +
  # geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  facet_grid(cols = vars(manag), rows = vars(fire)) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(color = "white"),
        strip.background = element_rect(fill = "black", color = "black")) +
  ylab("Observed growth") +
  xlab("Fitted growth")

ggsave("figures/growth_residuals_03.png",
       width = 14, height = 14, units = "cm")

# Spatial correlation check ------------------------------------------------

vplot <- vplot[vplot$Name %in% dplot$plot, ]
vplot <- vplot[order(vplot$Name), ]

# compute spatial distance across observations
spat_dist <- distance(vplot, symmetrical = TRUE) %>% as.matrix

e_plot <- as.matrix(m1, "e_plot") |> colMeans()
ranef_dist <- as.matrix(dist(e_plot))

dist_df <- data.frame(raneff_diff = ranef_dist[lower.tri(ranef_dist)],
                      distance = spat_dist[lower.tri(spat_dist)] / 1000)


# Aggreagate in 100 bins
nb <- 100
ng1 <- floor(nrow(dist_df) / nb)
ng2 <- nrow(dist_df) %% nb
group <- c(rep(1:nb, each = ng1), rep(nb, ng2))
# length(group) == nrow(dist_df) OK
dist_df <- dist_df[order(dist_df$distance), ]
dist_df$group <- group
dist_agg <- aggregate(cbind(raneff_diff, distance) ~ group, dist_df, mean)

# Residuals plots ---------------------------------------------------------

res_pforest <-
  ggplot(dtree, aes(x = pforest, y = res)) +
  geom_hline(yintercept = c(0.25, 0.5, 0.75), linetype = "dashed",
             linewidth = 0.3, color = viridis(1)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 10, bs = "cr"),
              method.args = list(family = mgcv::betar()),
              color = viridis(1, begin = 0.4),
              fill = viridis(1, begin = 0.4)) +
  facet_wrap(vars(manag), nrow = 1) +
  scale_y_continuous(expand = c(0.01, 0.01), limits = c(0, 1)) +
  xlab("Tree cover (%)") +
  ylab("DHARMa residuals") +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(color = "white", fill = "white"),
        panel.spacing = unit(3, "mm"),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        strip.text = element_text(size = 10)) +
  ggtitle("A")
res_pforest

res_elev <-
  ggplot(dtree, aes(x = elev, y = res)) +
  geom_hline(yintercept = c(0.25, 0.5, 0.75), linetype = "dashed",
             linewidth = 0.3, color = viridis(1)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 10, bs = "cr"),
              method.args = list(family = mgcv::betar()),
              color = viridis(1, begin = 0.4),
              fill = viridis(1, begin = 0.4)) +
  scale_y_continuous(expand = c(0.01, 0.01), limits = c(0, 1)) +
  facet_wrap(vars(manag), nrow = 1) +
  xlab("Elevation (m a.s.l.)") +
  ylab("DHARMa residuals") +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.spacing = unit(3, "mm"),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10)) +
  ggtitle("B")
res_elev

res_height <-
  ggplot(dtree, aes(x = height, y = res)) +
  geom_hline(yintercept = c(0.25, 0.5, 0.75), linetype = "dashed",
             linewidth = 0.3, color = viridis(1)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 10, bs = "cr"),
              method.args = list(family = mgcv::betar()),
              color = viridis(1, begin = 0.4),
              fill = viridis(1, begin = 0.4)) +
  scale_y_continuous(expand = c(0.01, 0.01), limits = c(0, 1)) +
  facet_wrap(vars(manag), nrow = 1) +
  xlab("Tree height (cm)") +
  ylab("DHARMa residuals") +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.spacing = unit(3, "mm"),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10)) +
  ggtitle("C")
res_height

res_fire <-
  ggplot(dtree, aes(x = fire, y = res)) +
  geom_hline(yintercept = c(0.25, 0.5, 0.75), linetype = "dashed",
             linewidth = 0.3, color = viridis(1)) +
  geom_boxplot(width = 0.5, fill = NA) +
  scale_y_continuous(expand = c(0.01, 0.01), limits = c(0, 1)) +
  facet_wrap(vars(manag), nrow = 1) +
  xlab("Fire") +
  ylab("DHARMa residuals") +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.spacing = unit(3, "mm"),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10)) +
  ggtitle("D")
res_fire

yname <- textGrob("Growth DHARMa residuals",
                  gp = gpar(fontsize = 10), rot = 90)

res_all <-
  deeptime::ggarrange2(res_pforest, res_elev, res_height, res_fire,
                       ncol = 1, left = yname)

ggsave("figures/growth_residuals_01.png", plot = res_all,
       width = 17, height = 24, units = "cm")

res_spat <-
  ggplot(dist_agg, aes(distance, raneff_diff)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 10),
              color = viridis(1, begin = 0.4),
              fill = viridis(1, begin = 0.4),
              linewidth = 0.7,
              n = 150) +
  scale_y_continuous(expand = c(0.01, 0.01),
                     limits = c(0, max(dist_agg$raneff_diff) * 1.1)) +
  xlab("Distance (km)") +
  ylab("Random effects difference between plots") +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        plot.title = element_text(hjust = 0),
        panel.grid.minor = element_blank())
res_spat

ggsave("figures/growth_residuals_02.png", plot = res_spat,
       width = 11, height = 10, units = "cm")

# Predictions growth ~ predictors ----------------------------------------

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
    varying_var = "Tree cover (%)"
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

# add manag_fire
pdata$manag_fire <- paste(pdata$manag, pdata$fire, sep = "_")
pdata$manag_fire <- factor(
  pdata$manag_fire, levels = levels(dplot$manag_fire)
)

# add variables at original scale
pdata$varying_val <- 0

rows_replace <- pdata$varying_var == "Elevation (m a.s.l.)"
pdata$varying_val[rows_replace] <-
  pdata$elev_z[rows_replace] * pred_sds["elev"] + pred_means["elev"]

rows_replace <- pdata$varying_var == "Tree cover (%)"
pdata$varying_val[rows_replace] <-
  pdata$pforest_z[rows_replace] * pred_sds["pforest"] + pred_means["pforest"]

rows_replace <- pdata$varying_var == "Tree height (cm)"
pdata$varying_val[rows_replace] <-
  pdata$height_z[rows_replace] * pred_sds["height"] + pred_means["height"]

Xpred <- model.matrix(fixed_formula, pdata)
Xpred_s <- model.matrix(fixed_formula_sigma, pdata)

# extract parameters
bhat <- as.matrix(m1, "b") %>% t
bs <- as.matrix(m1, "bs") %>% t
nsim <- length(sigma_plot)
muhat <- Xpred %*% bhat

# summarize posterior
psumm <- apply(muhat, 1, mean_ci) %>% t %>% as.data.frame
predictions <- cbind(pdata, psumm)

vv <- unique(pdata$varying_var)
vv_names <- unique(pdata$varying_var)
predictions$managlet <- factor(predictions$manag,
                               levels = levels(predictions$manag),
                               labels = c("A. Ranching",
                                          "B. Conservation"))
saveRDS(predictions, file.path("files", "growth_predictions.rds"))
predictions <- readRDS(file.path("files", "growth_predictions.rds"))



plist <- vector("list", 3)

breaks <- list(
  seq(1200, 2400, by = 300),
  seq(0, 100, by = 25),
  seq(200, 1200, by = 200)
)

# summarized data to plot alongside predictions
dtree$manag_fire <- factor(
  paste(dtree$manag, dtree$fire, sep = "_"),
  levels = unique(paste(dtree$manag, dtree$fire, sep = "_"))
)

vars <- c("elev", "pforest", "height")
dlist <- vector("list", 3)
for(i in 1:3) {
  # i = 3
  dd <- average(dtree, response = "growth", var = vars[i], O = 15)
  dd$managlet <- factor(dd$manag,
                        levels = levels(dd$manag),
                        labels = c("A. Ranching",
                                   "B. Conservation"))
  dlist[[i]] <- dd
}
lapply(dlist, function(x) table(x$manag, x$fire))
saveRDS(dlist, "files/growth_data_summary.rds")


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
    geom_point(data = dlist[[v]],
               mapping = aes(x = var, y = response,
                             color = fire, fill = fire),
               inherit.aes = F, alpha = 0.7, size = 0.8) +
    scale_color_viridis(discrete = T, option = "B", end = 0.35) +
    scale_fill_viridis(discrete = T, option = "B", end = 0.35) +
    facet_wrap(~ managlet, ncol = 2) +
    # scale_y_continuous(limits = c(50, 100), expand = c(0.01, 0.01)) +
    scale_x_continuous(breaks = breaks[[v]],
                       limits = range(c(breaks[[v]], d$varying_val))) +
    xlab(vv_names[v]) +
    ylab("Annual growth (height, cm)") +
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

growth_plot <- ggarrange2(plots = plist, nrow = 3)
ggsave("figures/growth_predictions.png", plot = growth_plot,
       width = 16, height = 15, units = "cm")

# Predictions growth ~ management --------------------------------------

pfire <- readRDS(file.path("files", "fire_prob_samples.rds"))

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
mu_manag <- Xpred_manag %*% bhat

# order pfire in the same way as pdata_manag
pfire_rep <- rbind(1 - t(pfire), t(pfire))
prods <- mu_manag * pfire_rep

means_manag <- aggregate(prods ~ manag, pdata_manag, sum)
dim(means_manag)
rownames(means_manag) <- means_manag[, 1]
means_manag_hat <- as.matrix(means_manag[, -1])
means_manag_summ <- apply(means_manag_hat, 1, mean_ci) %>% t

format(round(means_manag_summ, 2), nsmall = 2)
# mean    lower   upper
# Ranching     " 0.71" "-0.92" " 2.33"
# Conservation "-0.54" "-1.99" " 0.89"

diffmeans <- means_manag_hat[2, ] - means_manag_hat[1, ]
plot(density(diffmeans)); abline(v = 0, lty = 2)
format(round(mean_ci(diffmeans), 2), nsmall = 2)
# mean   lower   upper
# "-1.24" "-3.24" " 0.81"

# Probability of being smaller in conservation
(probb <- sum(means_manag_hat[2, ] < means_manag_hat[1, ]) / ncol(means_manag_hat))
# prob(exclusion < rangeland) = 0.8835

# save table
eeexp <- rbind(means_manag_summ, mean_ci(diffmeans),
               matrix(c(probb, NA, NA), 1))
rownames(eeexp)[3:4] <- c("diff", "prob")

write.csv(eeexp, file.path("files", "table_summary_growth_manag_marginal.csv"))