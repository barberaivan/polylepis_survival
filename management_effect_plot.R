library(tidyverse)
library(viridis)
theme_set(theme_classic())

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

vir_option <- "D"
vir_end <- 0.4


surv <- read.csv(file.path("files", "table_summary_survival_manag_marginal.csv"))
growth <- read.csv(file.path("files", "table_summary_growth_manag_marginal.csv"))
die <- read.csv(file.path("files", "table_summary_dieback_manag_marginal.csv"))

bb <- rbind(surv[1:2, -1] * 100, growth[1:2, -1], die[1:2, -1] * 100) %>% as.data.frame
bb$response <- factor(rep(c("1. Annual survival (%)",
                            "2. Annual growth (cm)",
                            "3. Dieback (%)"), each = 2),
                      levels = c(
                        "1. Annual survival (%)",
                        "2. Annual growth (cm)",
                        "3. Dieback (%)"
                      ))
bb$manag <- factor(rep(c("Ranching", "Conservation"), 3),
                   levels = c("Ranching", "Conservation"))


# probs
pp <- data.frame(
  manag = 1.5,
  probs = round(c(surv[4, 2], growth[4, 2], die[4, 2]), 3),
  response = levels(bb$response),
  yval = c(99.75, 2.4, 22.5)
)


# plot
ggplot(bb, aes(manag, mean, ymin = lower, ymax = upper, color = manag, fill = manag)) +
  geom_linerange(linewidth = 0.8) +
  geom_point(size = 3) +
  geom_text(data = pp, mapping = aes(x = manag, y = yval, label = probs),
            inherit.aes = F, size = 2.8, alpha = 0.75) +
  scale_color_viridis(discrete = TRUE, end = vir_end, option = vir_option) +
  facet_wrap(vars(response), scales = "free_y",
             strip.position = "left",
             axes = "all", axis.labels = "margins") +
  xlab("Land management") +
  nice_theme() +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = "white", color = "white"),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 8),
        strip.text = element_text(size = 9, color = "black"),
        axis.title.x = element_text(size = 9,
                                    margin = margin(2, unit = "mm")))

ggsave("figures/management_effects.png",
       width = 16, height = 6, units = "cm")