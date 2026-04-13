# rerun_problem_cells_checkpoint.R
# Rerun ONLY discrete cells with eta_min==0 (or complete==FALSE) WITH CHECKPOINTS

library(dplyr)

# ========= USER SETTINGS =========
IN_CSV <- "anova_kw_full_results.csv"   # <- mevcut H1 özet csv dosyan
OUT_CSV <- "h1_rerun_problem_cells_fixed.csv"

CKPT_DIR <- "rerun_ckpt"
dir.create(CKPT_DIR, showWarnings = FALSE, recursive = TRUE)

ALPHA <- 0.05
TARGET_ACCEPTS <- 10000
MAX_ATTEMPTS <- 5e6        # Poisson/NB için gerekebilir
SAVE_EVERY <- 1000          
SEED_BASE <- 20260227

f_target <- 0.6929735
tol <- 0.10
f_min <- f_target * (1 - tol)
f_max <- f_target * (1 + tol)

# ========= HELPERS =========
make_group_sizes <- function(n, k, balance) {
  if (balance == "balanced") rep(n, k) else round(seq(n - 1, n + 1, length.out = k))
}

# ======== YOUR ORIGINAL generate_data (paste/baseline) ========
generate_data <- function(group_sizes, dist, var_type, f_val) {
  
  k <- length(group_sizes)
  N <- sum(group_sizes)
  
  eta2_target <- (f_val^2) / (1 + f_val^2)
  
  group_means <- seq(-2, 2, length.out = k)
  group_means <- group_means / sd(group_means)
  
  sigma_base <- 1
  sigmas <- if (var_type == "homogeneous") rep(sigma_base, k) else sigma_base * seq(1, 2, length.out = k)
  
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

# ========= CORE: run one cell with checkpoint =========
run_one_cell_ckpt <- function(n, k, balance, dist, var_type) {
  
  cell_id <- sprintf("%d_%d_%s_%s_%s", n, k, balance, dist, var_type)
  ckpt_path <- file.path(CKPT_DIR, paste0(cell_id, ".rds"))
  
  set.seed(SEED_BASE + as.integer(n) * 10000 + as.integer(k) * 1000 + nchar(balance) * 100 + nchar(dist) * 10)
  
  group_sizes <- make_group_sizes(n, k, balance)
  
  eta2_min <- (f_min^2) / (1 + f_min^2)
  eta2_max <- (f_max^2) / (1 + f_max^2)
  
  # ---- load checkpoint if exists
  if (file.exists(ckpt_path)) {
    st <- readRDS(ckpt_path)
    acc <- st$acc
    attempts <- st$attempts
    rej_anova <- st$rej_anova
    rej_kw <- st$rej_kw
    rej_welch <- st$rej_welch
    eta_vals <- st$eta_vals
    welch_fail <- st$welch_fail
    message("Resuming from checkpoint: ", cell_id, " | acc=", acc, " attempts=", attempts)
  } else {
    acc <- 0L
    attempts <- 0L
    rej_anova <- logical(TARGET_ACCEPTS)
    rej_kw <- logical(TARGET_ACCEPTS)
    rej_welch <- rep(NA, TARGET_ACCEPTS)         # keep NA for failures
    eta_vals <- rep(NA_real_, TARGET_ACCEPTS)    # BUGFIX: NA not 0
    welch_fail <- 0L
    message("Starting new: ", cell_id)
  }
  
  # ---- main loop
  while (acc < TARGET_ACCEPTS && attempts < MAX_ATTEMPTS) {
    
    attempts <- attempts + 1L
    f_try <- runif(1, f_min, f_max)
    dat <- generate_data(group_sizes, dist, var_type, f_try)
    
    fit <- tryCatch(aov(y ~ group, dat), error = function(e) NULL)
    if (is.null(fit)) next
    
    ss <- summary(fit)[[1]]$`Sum Sq`
    eta2 <- ss[1] / sum(ss)
    if (!is.finite(eta2)) next
    
    if (eta2 >= eta2_min && eta2 <= eta2_max) {
      acc <- acc + 1L
      eta_vals[acc] <- eta2
      
      p_a <- summary(fit)[[1]]$`Pr(>F)`[1]
      rej_anova[acc] <- is.finite(p_a) && (p_a < ALPHA)
      
      p_w <- tryCatch(oneway.test(y ~ group, dat)$p.value, error = function(e) NA_real_)
      if (is.finite(p_w)) {
        rej_welch[acc] <- (p_w < ALPHA)
      } else {
        rej_welch[acc] <- NA
        welch_fail <- welch_fail + 1L
      }
      
      p_k <- tryCatch(kruskal.test(y ~ group, dat)$p.value, error = function(e) NA_real_)
      rej_kw[acc] <- is.finite(p_k) && (p_k < ALPHA)
      
      # ---- checkpoint save
      if (acc %% SAVE_EVERY == 0L) {
        saveRDS(list(
          acc = acc,
          attempts = attempts,
          rej_anova = rej_anova,
          rej_kw = rej_kw,
          rej_welch = rej_welch,
          eta_vals = eta_vals,
          welch_fail = welch_fail
        ), ckpt_path)
        message("Checkpoint saved: ", cell_id, " | acc=", acc, " attempts=", attempts)
      }
    }
  }
  
  complete <- (acc == TARGET_ACCEPTS)
  
  # ---- final save checkpoint (always)
  saveRDS(list(
    acc = acc,
    attempts = attempts,
    rej_anova = rej_anova,
    rej_kw = rej_kw,
    rej_welch = rej_welch,
    eta_vals = eta_vals,
    welch_fail = welch_fail
  ), ckpt_path)
  
  # ---- summaries using accepted only
  if (acc == 0L) {
    return(tibble(
      n = n, k = k, balance = balance, dist = dist, var_type = var_type,
      attempts_total = attempts, accepted = acc, complete = complete,
      acceptance_rate = 0,
      eta_min = NA_real_, eta_mean = NA_real_, eta_max = NA_real_,
      power_anova = NA_real_, power_welch = NA_real_, power_kw = NA_real_,
      welch_fail_rate = NA_real_
    ))
  }
  
  welch_ok <- !is.na(rej_welch[1:acc])
  power_w <- if (any(welch_ok)) mean(rej_welch[1:acc][welch_ok]) else NA_real_
  
  tibble(
    n = n, k = k, balance = balance, dist = dist, var_type = var_type,
    attempts_total = attempts,
    accepted = acc,
    complete = complete,
    acceptance_rate = acc / attempts,
    eta_min = min(eta_vals[1:acc], na.rm = TRUE),
    eta_mean = mean(eta_vals[1:acc], na.rm = TRUE),
    eta_max = max(eta_vals[1:acc], na.rm = TRUE),
    power_anova = mean(rej_anova[1:acc]),
    power_welch = power_w,
    power_kw = mean(rej_kw[1:acc]),
    welch_fail_rate = welch_fail / acc
  )
}

# ========= load current csv and find only problematic discrete cells =========
df <- read.csv(IN_CSV, stringsAsFactors = FALSE)

has_complete <- "complete" %in% names(df)

problem <- df %>%
  filter(dist %in% c("poisson", "negbin")) %>%
  filter(if (has_complete) (complete == FALSE) else (eta_min == 0))

if (nrow(problem) == 0) {
  message("No problem cells found.")
  quit(save = "no")
}

message("Problem cells to rerun: ", nrow(problem))
print(problem %>% select(n, k, balance, dist, var_type, eta_min, acceptance_rate) %>% head(50))

# ========= rerun each problem cell with checkpoint =========
rerun_res <- bind_rows(lapply(seq_len(nrow(problem)), function(i) {
  r <- problem[i, ]
  message(sprintf("Rerun %d/%d: n=%d k=%d %s %s %s",
                  i, nrow(problem), r$n, r$k, r$balance, r$dist, r$var_type))
  run_one_cell_ckpt(r$n, r$k, r$balance, r$dist, r$var_type)
}))

write.csv(rerun_res, OUT_CSV, row.names = FALSE)
message("Saved rerun results to: ", OUT_CSV)