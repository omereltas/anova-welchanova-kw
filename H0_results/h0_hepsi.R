############################################################
## H0 Monte Carlo: ANOVA – Welch – Kruskal–Wallis
## Full Design | Type I Error | eta² range
############################################################

library(car)
library(dplyr)
library(openxlsx)

set.seed(123456789)

ALPHA <- 0.05
TARGET_ACCEPTS <- 10000
MAX_ATTEMPTS <- 1e6

RESULT_DIR <- "H0_results"
dir.create(RESULT_DIR, showWarnings = FALSE)

## =========================
## Design (H1 ile AYNI)
## =========================
design <- expand.grid(
  n = c(5, 6, 8, 10),
  groups = c(3, 5, 8),
  balance = c("balanced", "unbalanced"),
  dist = c("normal", "exponential", "lognormal", "poisson", "negbin"),
  var_type = c("homogeneous", "heterogeneous"),
  stringsAsFactors = FALSE
)

## =========================
## Veri üretimi (H0)
## =========================
generate_data_H0 <- function(n, k, balance, dist, var_type) {
  
  if (balance == "balanced") {
    ns <- rep(n, k)
  } else {
    ns <- seq(n - 1, n + 1, length.out = k)
    ns <- round(ns)
  }
  
  sd_vec <- if (var_type == "homogeneous") {
    rep(1, k)
  } else {
    seq(1, 2, length.out = k)
  }
  
  y <- c()
  g <- c()
  
  for (i in 1:k) {
    ni <- ns[i]
    sdi <- sd_vec[i]
    
    yi <- switch(dist,
                 normal = rnorm(ni, 0, sdi),
                 exponential = rexp(ni, rate = 1 / sdi),
                 lognormal = rlnorm(ni, 0, sdi),
                 poisson = rpois(ni, lambda = 3),
                 negbin = rnbinom(ni, mu = 3, size = 1)
    )
    
    y <- c(y, yi)
    g <- c(g, rep(i, ni))
  }
  
  data.frame(y = y, group = factor(g))
}

## =========================
## Simülasyon (tek senaryo)
## =========================
run_H0_scenario <- function(row) {
  
  scenario_id <- paste(
    row$n, row$groups, row$balance,
    row$dist, row$var_type,
    sep = "_"
  )
  
  checkpoint_file <- file.path(
    RESULT_DIR,
    paste0("H0_", scenario_id, ".rds")
  )
  
  if (file.exists(checkpoint_file)) {
    return(readRDS(checkpoint_file))
  }
  
  accepts <- 0
  attempts <- 0
  
  rej_a <- logical(TARGET_ACCEPTS)
  rej_w <- logical(TARGET_ACCEPTS)
  rej_k <- logical(TARGET_ACCEPTS)
  eta_vals <- numeric(TARGET_ACCEPTS)
  
  while (accepts < TARGET_ACCEPTS && attempts < MAX_ATTEMPTS) {
    
    attempts <- attempts + 1
    dat <- generate_data_H0(
      n = row$n,
      k = row$groups,
      balance = row$balance,
      dist = row$dist,
      var_type = row$var_type
    )
    
    fit <- aov(y ~ group, dat)
    ss <- summary(fit)[[1]]$`Sum Sq`
    eta2 <- ss[1] / sum(ss)
    
    accepts <- accepts + 1
    eta_vals[accepts] <- eta2
    
    rej_a[accepts] <- summary(fit)[[1]]$`Pr(>F)`[1] < ALPHA
    rej_w[accepts] <- oneway.test(y ~ group, dat)$p.value < ALPHA
    rej_k[accepts] <- kruskal.test(y ~ group, dat)$p.value < ALPHA
  }
  
  out <- list(
    scenario_id = scenario_id,
    attempts = attempts,
    eta_min = min(eta_vals),
    eta_mean = mean(eta_vals),
    eta_max = max(eta_vals),
    type1_anova = mean(rej_a),
    type1_welch = mean(rej_w),
    type1_kw = mean(rej_k)
  )
  
  saveRDS(out, checkpoint_file)
  out
}

## =========================
## TÜM senaryolar
## =========================
results <- design %>%
  rowwise() %>%
  mutate(res = list(run_H0_scenario(cur_data()))) %>%
  ungroup() %>%
  tidyr::unnest_wider(res)

saveRDS(results, file.path(RESULT_DIR, "H0_FINAL_RESULTS_ALL_SCENARIOS.rds"))
write.xlsx(results, file.path(RESULT_DIR, "H0_FINAL_RESULTS_ALL_SCENARIOS.xlsx"))
View(final_results)
cat("H0 simulations completed correctly.\n")
