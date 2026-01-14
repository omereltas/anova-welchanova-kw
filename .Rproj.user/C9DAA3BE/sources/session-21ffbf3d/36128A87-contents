############################################################
## FULL MONTE CARLO SIMULATION
## ANOVA – Welch ANOVA – Kruskal Wallis
############################################################

## =========================
## 0. Paketler
## =========================
required_packages <- c("dplyr", "car")

for (p in required_packages) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, dependencies = TRUE)
  }
}
library(dplyr)
library(car)

## =========================
## 1. Global Ayarlar
## =========================
set.seed(123456)

ALPHA <- 0.05
TARGET_ACCEPTS <- 10000
MAX_ATTEMPTS <- 1e6
CHECKPOINT <- 1000

f_target <- 0.6929735
tol <- 0.10
f_min <- f_target * (1 - tol)
f_max <- f_target * (1 + tol)

RESULTS_DIR <- "anova_kw_results"
dir.create(RESULTS_DIR, showWarnings = FALSE)

## =========================
## 2. Grup Yapısı
## =========================
make_group_sizes <- function(n, k, balance) {
  if (balance == "balanced") {
    rep(n, k)
  } else {
    base <- rep(n, k)
    base[1] <- n - 1
    base[k] <- n + 1
    base
  }
}

## =========================
## 3. Veri Üretimi
## =========================
generate_data <- function(group_sizes, dist, var_type, f_val) {
  
  k <- length(group_sizes)
  N <- sum(group_sizes)
  
  eta2_target <- (f_val^2) / (1 + f_val^2)
  
  group_means <- seq(-2, 2, length.out = k)
  group_means <- group_means / sd(group_means)
  
  sigma_base <- 1
  sigmas <- if (var_type == "homogeneous") {
    rep(sigma_base, k)
  } else {
    sigma_base * seq(1, 2, length.out = k)
  }
  
  ss_between <- eta2_target * (N - 1)
  scaling <- sqrt(ss_between / sum(group_sizes * group_means^2))
  group_means <- group_means * scaling
  
  y <- c()
  g <- c()
  
  for (i in 1:k) {
    n_i <- group_sizes[i]
    mu <- group_means[i]
    sd_i <- sigmas[i]
    
    vals <- switch(
      dist,
      normal = rnorm(n_i, mu, sd_i),
      exponential = rexp(n_i, rate = 1 / sd_i) + mu,
      lognormal = rlnorm(n_i, meanlog = mu, sdlog = sd_i),
      poisson = rpois(n_i, lambda = abs(mu) + 3),
      negbin = rnbinom(n_i, mu = abs(mu) + 3, size = 1)
    )
    
    y <- c(y, vals)
    g <- c(g, rep(i, n_i))
  }
  
  data.frame(y = y, group = factor(g))
}

## =========================
## 4. Tek Senaryo Simülasyonu
## =========================
run_scenario <- function(n, k, balance, dist, var_type) {
  
  id <- paste(n, k, balance, dist, var_type, sep = "_")
  file_path <- file.path(RESULTS_DIR, paste0(id, ".rds"))
  
  if (file.exists(file_path)) {
    message("Loaded cached: ", id)
    return(readRDS(file_path))
  }
  
  acc <- 0
  attempts <- 0
  
  rej_anova <- logical(TARGET_ACCEPTS)
  rej_welch <- logical(TARGET_ACCEPTS)
  rej_kw <- logical(TARGET_ACCEPTS)
  
  eta_vals <- numeric(TARGET_ACCEPTS)
  
  group_sizes <- make_group_sizes(n, k, balance)
  
  while (acc < TARGET_ACCEPTS && attempts < MAX_ATTEMPTS) {
    
    attempts <- attempts + 1
    f_try <- runif(1, f_min, f_max)
    
    dat <- generate_data(group_sizes, dist, var_type, f_try)
    
    fit <- aov(y ~ group, dat)
    ss <- summary(fit)[[1]]$`Sum Sq`
    eta2 <- ss[1] / sum(ss)
    
    if (eta2 >= (f_min^2 / (1 + f_min^2)) &&
        eta2 <= (f_max^2 / (1 + f_max^2))) {
      
      acc <- acc + 1
      eta_vals[acc] <- eta2
      
      rej_anova[acc] <- summary(fit)[[1]]$`Pr(>F)`[1] < ALPHA
      rej_welch[acc] <- oneway.test(y ~ group, dat)$p.value < ALPHA
      rej_kw[acc] <- kruskal.test(y ~ group, dat)$p.value < ALPHA
      
      if (acc %% CHECKPOINT == 0) {
        saveRDS(list(
          acc = acc,
          attempts = attempts,
          eta_vals = eta_vals[1:acc],
          anova = rej_anova[1:acc],
          welch = rej_welch[1:acc],
          kw = rej_kw[1:acc]
        ), file_path)
      }
    }
  }
  
  result <- list(
    attempts_total = attempts,
    acceptance_rate = acc / attempts,
    eta_min = min(eta_vals),
    eta_mean = mean(eta_vals),
    eta_max = max(eta_vals),
    power_anova = mean(rej_anova),
    power_welch = mean(rej_welch),
    power_kw = mean(rej_kw)
  )
  
  saveRDS(result, file_path)
  result
}

## =========================
## 5. Tam Faktöriyel Tasarım
## =========================
design <- expand.grid(
  n = c(5, 6, 8, 10),
  k = c(3, 5, 8),
  balance = c("balanced", "unbalanced"),
  dist = c("normal", "exponential", "lognormal", "poisson", "negbin"),
  var_type = c("homogeneous", "heterogeneous"),
  stringsAsFactors = FALSE
)

## =========================
## 6. Tüm Senaryolar
## =========================
final_results <- design %>%
  rowwise() %>%
  mutate(
    res = list(run_scenario(n, k, balance, dist, var_type)),
    power_anova = res$power_anova,
    power_welch = res$power_welch,
    power_kw = res$power_kw,
    eta_min = res$eta_min,
    eta_max = res$eta_max,
    acceptance_rate = res$acceptance_rate
  ) %>%
  ungroup() %>%
  select(-res)

write.csv(final_results,
          file.path(RESULTS_DIR, "anova_kw_full_results.csv"),
          row.names = FALSE)

message("ALL SIMULATIONS COMPLETED SUCCESSFULLY.")
