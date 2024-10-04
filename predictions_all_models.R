library(tidyverse)
library(viridis)
library(ggh4x)
library(deeptime)
library(grid)
theme_set(theme_classic())

# Functions ---------------------------------------------------------------

nice_theme <- function() {
  theme(
    panel.border = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major =  element_line(),

    axis.line = element_line(linewidth = 0.25),

    axis.ticks = element_line(linewidth = 0.25),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 10),

    strip.text = element_text(size = 10),
    strip.background = element_rect(fill = "white", color = "white")
  )
}

# function to summarize data by group, taking into account that groups should
# have at least O observations.
# response is the response, var is the predictor
average <- function(data, response, var, O = 3) {
  # data = d_surv; var = "surv"; response = "pforest"; O = 3

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

# Load predictions and data -----------------------------------------------

pdata_list <- list(
  survival = readRDS("files/survival_predictions.rds"),
  growth = readRDS("files/growth_predictions.rds"),
  dieback =  readRDS("files/dieback_predictions.rds")
)

# survival data with predicted height for unidentified trees
surv_extra <- readRDS(file.path("files", "survival_model_residuals_and_data_elevLinear.rds"))

data_list <- list(
  survival = surv_extra$data_both,
  growth = readRDS("files/growth_data_dtree.rds"),
  dieback = readRDS("files/dieback_data_dtree.rds")
)

# Constants ---------------------------------------------------------------

response_names <- names(data_list)
response_names2 <- c("surv", "growth", "dieback")
response_ylab <- c("1. Annual survival (%)", "2. Annual growth (cm)",
                   "3. Dieback (%)")

pred_names <- c("pforest", "elev", "height")
pred_xlab <- c("Tree cover (%)", "Elevation (m a.s.l.)", "Tree height (cm)")

# breaks for predictors
breaks <- list(
  "pforest" = seq(0, 100, by = 25),
  "elevation" = seq(1200, 2400, by = 300),
  "height" = seq(200, 800, by = 200)
)

vir_option <- "C"
vir_end <- 0.4

# number of plots to aggregate by point
O <- c(3, 3, 3)

# Generalized code --------------------------------------------------------

plots_list <- lapply(1:3, function(i) {
  l <- vector("list", 3)
  names(l) <- response_names
  return(l)
})
names(plots_list) <- pred_names

for(j in 1:3) {    # predictors
  for(i in 1:3) {  # responses

    # j = 1; i = 3
    print(paste("predictor(j) =", j, "; response(i) =", i))

    data00 <- data_list[[i]][, c("plot", "manag", "fire", "manag_fire",
                                 response_names2[i], pred_names[j])]
    colnames(data00) <- c("plot", "manag", "fire", "manag_fire",
                          "response", "predictor")
    # average by plot
    data0 <- aggregate(cbind(response, predictor) ~
                         plot + manag + fire + manag_fire,
                       data00, mean)

    # group plots
    datalocal <- average(data0, "response", "predictor", O = O[i])

    # scale survival proportion to annual percentage
    if(i == 1) {
      datalocal$response <- (datalocal$response ^ (1 / 15)) * 100
    }

    datalocal$managlet <- factor(datalocal$manag,
                                 levels = levels(datalocal$manag),
                                 labels = c("A. Ranching", "B. Conservation"))

    predlocal <- pdata_list[[i]]
    predlocal <- predlocal[predlocal$varying_var == pred_xlab[j], ]

    # ignore predictions for trees that are too tall
    if(j == 3) predlocal <- predlocal[predlocal$varying_val <= 1000, ]

    ylims <- range(c(datalocal$response, predlocal$upper, predlocal$lower))
    if(response_names[i] == "dieback") ylims[1] <- 0

    plotcito <-
      ggplot(predlocal,
             aes(x = varying_val, y = mean, ymin = lower, ymax = upper,
                 color = fire, fill = fire)) +
      geom_ribbon(color = NA, alpha = 0.25) +
      geom_line() +
      geom_point(data = datalocal, # d_surv_pforest,
                 mapping = aes(x = var, y = response,
                               color = fire, fill = fire),
                 inherit.aes = F, alpha = 0.7, size = 0.8) +

      scale_color_viridis(discrete = T, option = vir_option, end = vir_end) +
      scale_fill_viridis(discrete = T, option = vir_option, end = vir_end) +

      facet_wrap(~ managlet, ncol = 2, axes = "all", axis.labels = "margins") +

      scale_y_continuous(limits = ylims) +
      scale_x_continuous(breaks = breaks[[j]],
                         limits = range(c(breaks[[j]], predlocal$varying_val)),
                         expand = c(0.01, 0.01)) +

      xlab(pred_xlab[j]) +
      ylab(response_ylab[i]) +

      nice_theme() +
      theme(panel.grid.minor = element_blank(),
            strip.background = element_rect(color = "white", fill = "white"),
            panel.spacing = unit(4, "mm"),
            legend.title = element_blank(),
            legend.position = "none",
            axis.text = element_text(size = 8),
            axis.title = element_text(size = 10),
            strip.text = element_text(size = 10),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank())
    # plotcito

    if(i %in% 2:3) {
      plotcito <- plotcito + theme(strip.background = element_blank(),
                                   strip.text = element_blank())
    }

    if(i == 3) {
      plotcito <-
        plotcito +
          theme(
            axis.text.x = element_text(),
            axis.ticks.x = element_line(),
            axis.title.x = element_text(),
            legend.text = element_text(size = 9,
                                       margin = margin(1, 3, 1, 1, unit = "mm")),
            legend.position = "bottom",
            legend.box = "horizontal",
            legend.background = element_rect(fill = rgb(0, 0, 0, 0))
          )
    }

    # plotcito
    plots_list[[j]][[i]] <- plotcito
  }
}

ww <- 15; hh <- 17

fig_pforest <- ggarrange2(plots = plots_list[[1]], ncol = 1)
ggsave("figures/predictions_pforest.png",
       plot = fig_pforest,
       width = ww, height = hh, units = "cm")

fig_elev <- ggarrange2(plots = plots_list[[2]], ncol = 1)
ggsave("figures/predictions_elev.png",
       plot = fig_elev,
       width = ww, height = hh, units = "cm")

fig_height <- ggarrange2(plots = plots_list[[3]], ncol = 1)
ggsave("figures/predictions_height.png",
       plot = fig_height,
       width = ww, height = hh, units = "cm")