# rerun_incomplete_cells.R
# Sadece tamamlanmayan 2 negbin heterojen hücreyi çalıştırır

library(dplyr)

# ========= AYARLAR =========
OUT_CSV  <- "negbin_het_incomplete_fixed.csv"
CKPT_DIR <- "rerun_ckpt_fixed"  # Aynı checkpoint klasörü — kaldığı yerden devam eder

ALPHA          <- 0.05
TARGET_ACCEPTS <- 10000
MAX_ATTEMPTS   <- 15e6   # 15 milyon
SAVE_EVERY     <- 1000
SEED_BASE      <- 20260227

f_target <- 0.6929735
tol      <- 0.10
f_min    <- f_target * (1 - tol)
f_max    <- f_target * (1 + tol)

# ========= YARDIMCI FONKSİYONLAR =========
make_group_sizes <- function(n, k, balance) {
  if (balance == "balanced") {
    rep(n, k)
  } else {
    base    <- rep(n, k)
    base[1] <- n - 1
    base[k] <- n + 1
    base
  }
}

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
  
  ss_between  <- eta2_target * (N - 1)
  scaling     <- sqrt(ss_between / sum(group_sizes * group_means^2))
  group_means <- group_means * scaling
  
  y <- c()
  g <- c()
  
  for (i in 1:k) {
    n_i  <- group_sizes[i]
    mu   <- group_means[i]
    sd_i <- sigmas[i]
    
    vals <- switch(
      dist,
      normal      = rnorm(n_i, mu, sd_i),
      exponential = rexp(n_i, rate = 1 / sd_i) + mu,
      lognormal   = rlnorm(n_i, meanlog = mu, sdlog = sd_i),
      poisson     = rpois(n_i, lambda = abs(mu) + 3),
      negbin      = rnbinom(n_i, mu = abs(mu) + 3, size = 1)
    )
    
    y <- c(y, vals)
    g <- c(g, rep(i, n_i))
  }
  
  data.frame(y = y, group = factor(g))
}

make_seed <- function(n, k, balance, dist, var_type) {
  SEED_BASE        +
    as.integer(n)  * 100000 +
    as.integer(k)  *  10000 +
    nchar(balance) *   1000 +
    nchar(dist)    *    100 +
    nchar(var_type)*     10
}

# ========= ANA FONKSİYON =========
run_one_cell <- function(n, k, balance, dist, var_type) {
  
  cell_id   <- sprintf("%d_%d_%s_%s_%s", n, k, balance, dist, var_type)
  ckpt_path <- file.path(CKPT_DIR, paste0(cell_id, ".rds"))
  
  set.seed(make_seed(n, k, balance, dist, var_type))
  
  group_sizes <- make_group_sizes(n, k, balance)
  eta2_min    <- (f_min^2) / (1 + f_min^2)
  eta2_max    <- (f_max^2) / (1 + f_max^2)
  
  # Checkpoint varsa kaldığı yerden devam et
  if (file.exists(ckpt_path)) {
    st         <- readRDS(ckpt_path)
    acc        <- st$acc
    attempts   <- st$attempts
    rej_anova  <- st$rej_anova
    rej_welch  <- st$rej_welch
    rej_kw     <- st$rej_kw
    eta_vals   <- st$eta_vals
    welch_fail <- st$welch_fail
    message("Checkpoint yuklendi: ", cell_id,
            " | acc=", acc, " attempts=", attempts)
  } else {
    acc        <- 0L
    attempts   <- 0L
    rej_anova  <- logical(TARGET_ACCEPTS)
    rej_welch  <- rep(NA, TARGET_ACCEPTS)
    rej_kw     <- logical(TARGET_ACCEPTS)
    eta_vals   <- rep(NA_real_, TARGET_ACCEPTS)
    welch_fail <- 0L
    message("Yeni basliyor: ", cell_id)
  }
  
  while (acc < TARGET_ACCEPTS && attempts < MAX_ATTEMPTS) {
    
    attempts <- attempts + 1L
    f_try    <- runif(1, f_min, f_max)
    dat      <- generate_data(group_sizes, dist, var_type, f_try)
    
    fit <- tryCatch(aov(y ~ group, dat), error = function(e) NULL)
    if (is.null(fit)) next
    
    ss   <- summary(fit)[[1]]$`Sum Sq`
    eta2 <- ss[1] / sum(ss)
    if (!is.finite(eta2)) next
    
    if (eta2 >= eta2_min && eta2 <= eta2_max) {
      
      acc           <- acc + 1L
      eta_vals[acc] <- eta2
      
      p_a            <- summary(fit)[[1]]$`Pr(>F)`[1]
      rej_anova[acc] <- is.finite(p_a) && (p_a < ALPHA)
      
      p_w <- tryCatch(
        oneway.test(y ~ group, dat)$p.value,
        error = function(e) NA_real_
      )
      if (is.finite(p_w)) {
        rej_welch[acc] <- (p_w < ALPHA)
      } else {
        rej_welch[acc] <- NA
        welch_fail     <- welch_fail + 1L
      }
      
      p_k           <- tryCatch(
        kruskal.test(y ~ group, dat)$p.value,
        error = function(e) NA_real_
      )
      rej_kw[acc] <- is.finite(p_k) && (p_k < ALPHA)
      
      if (acc %% SAVE_EVERY == 0L) {
        saveRDS(
          list(acc=acc, attempts=attempts,
               rej_anova=rej_anova, rej_welch=rej_welch,
               rej_kw=rej_kw, eta_vals=eta_vals,
               welch_fail=welch_fail),
          ckpt_path
        )
        message("Checkpoint: ", cell_id,
                " | acc=", acc, " / attempts=", attempts)
      }
    }
  }
  
  # Son checkpoint
  saveRDS(
    list(acc=acc, attempts=attempts,
         rej_anova=rej_anova, rej_welch=rej_welch,
         rej_kw=rej_kw, eta_vals=eta_vals,
         welch_fail=welch_fail),
    ckpt_path
  )
  
  complete <- (acc == TARGET_ACCEPTS)
  if (!complete) warning("HALA TAMAMLANAMADI: ", cell_id,
                         " | Kabul: ", acc, "/", TARGET_ACCEPTS)
  
  welch_ok <- !is.na(rej_welch[1:acc])
  power_w  <- if (any(welch_ok)) mean(rej_welch[1:acc][welch_ok]) else NA_real_
  
  tibble(
    n=n, k=k, balance=balance, dist=dist, var_type=var_type,
    attempts_total  = attempts,
    accepted        = acc,
    complete        = complete,
    acceptance_rate = acc / attempts,
    eta_min         = min(eta_vals[1:acc],  na.rm=TRUE),
    eta_mean        = mean(eta_vals[1:acc], na.rm=TRUE),
    eta_max         = max(eta_vals[1:acc],  na.rm=TRUE),
    power_anova     = mean(rej_anova[1:acc]),
    power_welch     = power_w,
    power_kw        = mean(rej_kw[1:acc]),
    welch_fail_rate = welch_fail / acc
  )
}

# ========= SADECE 2 TAMAMLANMAYAN HÜCRE =========
target_cells <- data.frame(
  n        = c(10, 10),
  k        = c(8,   8),
  balance  = c("balanced", "unbalanced"),
  dist     = c("negbin",   "negbin"),
  var_type = c("heterogeneous", "heterogeneous"),
  stringsAsFactors = FALSE
)

message("Calistirilacak hucre sayisi: ", nrow(target_cells))

results <- bind_rows(lapply(seq_len(nrow(target_cells)), function(i) {
  r <- target_cells[i, ]
  message(sprintf("\n[%d/2] n=%d k=%d %s %s %s",
                  i, r$n, r$k, r$balance, r$dist, r$var_type))
  run_one_cell(r$n, r$k, r$balance, r$dist, r$var_type)
}))

write.csv(results, OUT_CSV, row.names = FALSE)
message("\nTamamlandi. Sonuclar: ", OUT_CSV)

# Kontrol
message("\n--- SONUC ---")
print(results[, c("n","k","balance","complete","accepted",
                  "acceptance_rate","power_anova","power_welch","power_kw")])