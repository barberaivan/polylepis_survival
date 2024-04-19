library(tidyverse)
library(viridis)

surv <- read.csv(file.path("files", "table_summary_survival_manag_marginal.csv"))
growth <- read.csv(file.path("files", "table_summary_growth_manag_marginal.csv"))
die <- read.csv(file.path("files", "table_summary_dieback_manag_marginal.csv"))

bb <- rbind(surv[1:2, -1] * 100, growth[1:2, -1], die[1:2, -1] * 100) %>% as.data.frame
bb$response <- factor(rep(c("A. Annual survival probability (%)",
                            "B. Annual growth (height, cm)",
                            "C. Dieback (%)"), each = 2),
                      levels = c(
                        "A. Annual survival probability (%)",
                        "B. Annual growth (height, cm)",
                        "C. Dieback (%)"
                      ))
bb$manag <- factor(rep(c("Rangeland", "Livestock\nexclusion"), 3),
                   levels = c("Livestock\nexclusion", "Rangeland"))


# probs
pp <- data.frame(
  manag = 1.5,
  probs = round(c(surv[4, 2], growth[4, 2], die[4, 2]), 3),
  response = levels(bb$response),
  yval = c(99.75, 2.2, 22.5)
)


# plot
ggplot(bb, aes(manag, mean, ymin = lower, ymax = upper, color = manag, fill = manag)) +
  geom_linerange(linewidth = 0.8) +
  geom_point(size = 3) +
  geom_text(data = pp, mapping = aes(x = manag, y = yval, label = probs),
            inherit.aes = F, size = 2.8, alpha = 0.75) +
  scale_color_viridis(discrete = TRUE, end = 0.5) +
  facet_wrap(vars(response), scales = "free_y",
             strip.position = "left") +
  xlab("Land management") +
  theme_bw() +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = "white", color = "white"),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 8),
        strip.text = element_text(size = 9, color = "black"),
        axis.title.x = element_text(size = 9))

ggsave("figures/management_effects.png",
       width = 16, height = 7, units = "cm")

