# figures_H1_v3.R

library(dplyr)
library(ggplot2)

df <- read.csv("H1_FINAL.csv", stringsAsFactors = FALSE)

# ========= ORTAK AYARLAR =========
dist_levels <- c("normal", "exponential", "lognormal", "poisson", "negbin")
dist_labels <- c("Normal", "Exponential", "Lognormal", "Poisson", "Neg. Binomial")

df <- df %>%
  mutate(
    dist = factor(dist, levels = dist_levels, labels = dist_labels),
    var_type = factor(var_type, levels = c("homogeneous","heterogeneous"),
                      labels = c("Homogeneous Variance","Heterogeneous Variance")),
    balance = factor(balance, levels = c("balanced","unbalanced"),
                     labels = c("Balanced","Unbalanced")),
    k_label = factor(paste0("k = ", k), levels = paste0("k = ", c(3,5,8))),
    n_num = as.numeric(n)
  )

# ========= RENK + SEKIL + CIZGI TIPI =========
dist_colors <- c(
  "Normal"        = "#0072B2",
  "Exponential"   = "#E69F00",
  "Lognormal"     = "#009E73",
  "Poisson"       = "#CC79A7",
  "Neg. Binomial" = "#56B4E9"
)

dist_shapes <- c(
  "Normal"        = 16,
  "Exponential"   = 17,
  "Lognormal"     = 15,
  "Poisson"       = 18,
  "Neg. Binomial" = 8
)

dist_linetypes <- c(
  "Normal"        = "solid",
  "Exponential"   = "dashed",
  "Lognormal"     = "dotted",
  "Poisson"       = "dotdash",
  "Neg. Binomial" = "longdash"
)

# ========= ORTAK TEMA =========
power_theme <- theme_bw(base_size = 13) +
  theme(
    strip.background = element_rect(fill = "white", color = "grey40"),
    strip.text = element_text(face = "bold", size = 10, color = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 11),
    legend.key = element_rect(fill = "white", color = NA),
    legend.key.width = unit(1.5, "cm"),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
  )

# ========= ORTAK CIZIM FONKSIYONU =========
make_power_figure <- function(data, power_col, title_text) {
  data %>%
    rename(power = !!power_col) %>%
    ggplot(aes(x = n_num, y = power, color = dist, shape = dist, linetype = dist, group = dist)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2.8) +
    facet_grid(
      rows = vars(var_type, balance),
      cols = vars(k_label)
    ) +
    scale_color_manual(values = dist_colors, name = "Distribution") +
    scale_shape_manual(values = dist_shapes, name = "Distribution") +
    scale_linetype_manual(values = dist_linetypes, name = "Distribution") +
    scale_x_continuous(breaks = c(5,6,8,10), labels = c("5","6","8","10")) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    guides(
      color = guide_legend(override.aes = list(size = 3, linewidth = 1.2), nrow = 1),
      shape = guide_legend(nrow = 1),
      linetype = guide_legend(nrow = 1)
    ) +
    labs(
      title = title_text,
      x = "Sample size per group (n)",
      y = "Empirical Power"
    ) +
    power_theme
}

# ========= FIGURE 2: ANOVA =========
fig2 <- make_power_figure(df, "power_anova",
                          "Empirical Power of Classical One-Way ANOVA")
ggsave("Figure2_ANOVA_Power.png", fig2, width = 13, height = 9.5, dpi = 300, bg = "white")

# ========= FIGURE 3: WELCH ANOVA =========
fig3 <- make_power_figure(df, "power_welch",
                          "Empirical Power of Welch's ANOVA")
ggsave("Figure3_Welch_Power.png", fig3, width = 13, height = 9.5, dpi = 300, bg = "white")

# ========= FIGURE 4: KRUSKAL-WALLIS =========
fig4 <- make_power_figure(df, "power_kw",
                          "Empirical Power of the Kruskal-Wallis Test")
ggsave("Figure4_KW_Power.png", fig4, width = 13, height = 9.5, dpi = 300, bg = "white")

message("Figurler kaydedildi.")