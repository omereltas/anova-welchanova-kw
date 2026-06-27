# figure1_type1error.R

library(dplyr)
library(ggplot2)
library(tidyr)
library(openxlsx)

# ========= VERI OKU =========
df_raw <- read.xlsx("H0_results/H0_FINAL_RESULTS_ALL_SCENARIOS.xlsx")

# "groups" sutununu "k" ye cevir
if ("groups" %in% names(df_raw)) {
  df_raw <- df_raw %>% rename(k = groups)
}

# ========= VERI HAZIRLA =========
dist_levels <- c("normal", "exponential", "lognormal", "poisson", "negbin")
dist_labels <- c("Normal", "Exponential", "Lognormal", "Poisson", "Neg. Binomial")

df <- df_raw %>%
  mutate(
    dist = factor(dist, levels = dist_levels, labels = dist_labels),
    var_type = factor(var_type,
                      levels = c("homogeneous", "heterogeneous"),
                      labels = c("Homogeneous Variance", "Heterogeneous Variance")),
    balance = factor(balance,
                     levels = c("balanced", "unbalanced"),
                     labels = c("Balanced", "Unbalanced")),
    k_label = factor(paste0("k = ", k), levels = paste0("k = ", c(3, 5, 8))),
    n_num   = as.numeric(n)
  )

# Uzun formata cevir
df_long <- df %>%
  select(n_num, k_label, balance, var_type, dist,
         type1_anova, type1_welch, type1_kw) %>%
  pivot_longer(
    cols      = c(type1_anova, type1_welch, type1_kw),
    names_to  = "test",
    values_to = "type1_error"
  ) %>%
  mutate(
    test = factor(test,
                  levels = c("type1_anova", "type1_welch", "type1_kw"),
                  labels = c("ANOVA", "Welch's ANOVA", "Kruskal-Wallis"))
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

# ========= TEMA =========
type1_theme <- theme_bw(base_size = 13) +
  theme(
    strip.background  = element_rect(fill = "white", color = "grey40"),
    strip.text        = element_text(face = "bold", size = 10, color = "black"),
    panel.grid.minor  = element_blank(),
    panel.grid.major  = element_line(color = "grey85", linewidth = 0.3),
    panel.background  = element_rect(fill = "white", color = NA),
    plot.background   = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.position   = "bottom",
    legend.title      = element_text(face = "bold", size = 12),
    legend.text       = element_text(size = 11),
    legend.key        = element_rect(fill = "white", color = NA),
    legend.key.width  = unit(1.5, "cm"),
    axis.text         = element_text(size = 10, color = "black"),
    axis.title        = element_text(size = 12),
    plot.title        = element_text(face = "bold", size = 14, hjust = 0.5)
  )

# ========= BRADLEY SINIRLARI =========
# Nominal alpha = 0.05
# Liberal  : 0.025 - 0.075 (alpha +/- %50)
# Strict   : 0.045 - 0.055 (alpha +/- %10)

bradley_liberal_lo <- 0.025
bradley_liberal_hi <- 0.075
bradley_strict_lo  <- 0.045
bradley_strict_hi  <- 0.055

# ========= FIGURE 1 =========
fig1 <- df_long %>%
  ggplot(aes(x = n_num, y = type1_error,
             color = dist, shape = dist, linetype = dist, group = dist)) +
  
  # Bradley liberal sinir (acik gri alan)
  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = bradley_liberal_lo, ymax = bradley_liberal_hi,
           fill = "grey85", alpha = 0.6) +
  
  # Bradley strict sinir (daha koyu gri alan)
  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = bradley_strict_lo, ymax = bradley_strict_hi,
           fill = "grey70", alpha = 0.6) +
  
  # Nominal alpha cizgisi
  geom_hline(yintercept = 0.05,
             linetype = "solid", color = "black", linewidth = 0.7) +
  
  # Veriler
  geom_line(linewidth = 1) +
  geom_point(size = 2.8) +
  
  facet_grid(
    rows = vars(var_type, balance),
    cols = vars(test)
  ) +
  
  scale_color_manual(values = dist_colors,      name = "Distribution") +
  scale_shape_manual(values = dist_shapes,       name = "Distribution") +
  scale_linetype_manual(values = dist_linetypes, name = "Distribution") +
  scale_x_continuous(breaks = c(5, 6, 8, 10),
                     labels = c("5", "6", "8", "10")) +
  scale_y_continuous(limits = c(0.00, 0.22),
                     breaks = c(0.00, 0.05, 0.10, 0.15, 0.20)) +
  
  guides(
    color    = guide_legend(
      override.aes = list(size = 3, linewidth = 1.2), nrow = 1),
    shape    = guide_legend(nrow = 1),
    linetype = guide_legend(nrow = 1)
  ) +
  
  labs(
    title = "Empirical Type I Error Rates Across Distributions and Design Conditions",
    x     = "Sample size per group (n)",
    y     = "Empirical Type I Error Rate"
  ) +
  type1_theme

ggsave("Figure1_Type1Error.png", fig1,
       width = 14, height = 9.5, dpi = 300, bg = "white")

message("Kaydedildi: Figure1_Type1Error.png")