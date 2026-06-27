# final_merge.R

library(dplyr)

main  <- read.csv("H1_FINAL_MERGED2.csv",             stringsAsFactors = FALSE)
rerun <- read.csv("poisson_negbin_rerun_clean.csv",   stringsAsFactors = FALSE)

keep_cols <- c("n","k","balance","dist","var_type",
               "power_anova","power_welch","power_kw",
               "eta_min","eta_max","acceptance_rate")

add_missing <- function(df, cols) {
  for (col in cols) if (!col %in% names(df)) df[[col]] <- NA
  df[, cols]
}

main_clean  <- add_missing(main,  keep_cols)
rerun_clean <- add_missing(rerun, keep_cols)

main_reduced <- main_clean %>%
  anti_join(rerun_clean, by = c("n","k","balance","dist","var_type"))

final <- bind_rows(main_reduced, rerun_clean) %>%
  arrange(dist, var_type, balance, k, n)

message("Satir sayisi: ", nrow(final), " (beklenen: 240)")

message("\n--- ETA_MIN=0 KONTROLU ---")
zero <- final[final$eta_min == 0, ]
if (nrow(zero) == 0) message("Temiz.") else print(zero)

message("\n--- WELCH NA KONTROLU ---")
na_w <- final[is.na(final$power_welch), ]
if (nrow(na_w) == 0) message("Temiz.") else {
  message(nrow(na_w), " adet NA var:")
  print(na_w[, c("n","k","balance","dist","var_type")])
}

write.csv(final, "H1_FINAL.csv", row.names = FALSE)
message("\nKaydedildi: H1_FINAL.csv")