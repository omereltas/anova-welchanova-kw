# merge_all_h1_results.R

library(dplyr)

# ========= DOSYALARI OKU =========
main   <- read.csv("anova_kw_results/anova_kw_full_results.csv", stringsAsFactors = FALSE)
rerun1 <- read.csv("anova_kw_results/h1_rerun_problem_cells_fixed.csv",                 stringsAsFactors = FALSE)
rerun2 <- read.csv("negbin_het_rerun_fixed.csv",                 stringsAsFactors = FALSE)
rerun3 <- read.csv("negbin_het_incomplete_fixed.csv",            stringsAsFactors = FALSE)

# ========= SUTUN STANDARDIZASYONU =========
if ("groups" %in% names(main)) {
  main <- main %>% rename(k = groups)
}

message("Ana dosya sutunlari: ", paste(names(main), collapse=", "))

# ========= ORTAK SUTUNLARA INDIR =========
keep_cols <- c("n","k","balance","dist","var_type",
               "power_anova","power_welch","power_kw",
               "eta_min","eta_max","acceptance_rate")

add_missing <- function(df, cols) {
  for (col in cols) {
    if (!col %in% names(df)) df[[col]] <- NA
  }
  df[, cols]
}

main_clean   <- add_missing(main,   keep_cols)
rerun2_clean <- add_missing(rerun2, keep_cols)
rerun3_clean <- add_missing(rerun3, keep_cols)

# ========= BIRLESTIRME =========
overrides <- bind_rows(rerun2_clean, rerun3_clean) %>%
  distinct(n, k, balance, dist, var_type, .keep_all = TRUE)

main_reduced <- main_clean %>%
  anti_join(overrides, by = c("n","k","balance","dist","var_type"))

final <- bind_rows(main_reduced, overrides) %>%
  arrange(dist, var_type, balance, k, n)

# ========= KONTROL =========
message("\nSatir sayilari:")
message("  Ana (orijinal):   ", nrow(main_clean))
message("  Override sayisi:  ", nrow(overrides))
message("  Ana (azaltilmis): ", nrow(main_reduced))
message("  Final:            ", nrow(final))

if (nrow(final) != nrow(main_clean)) {
  warning("SATIR SAYISI UYUSMUYOR! Beklenen: ", nrow(main_clean),
          " | Bulunan: ", nrow(final))
} else {
  message("Satir sayisi dogru.")
}

# Negbin heterojen kontrolu
message("\n--- NEGBIN HOMOJEN vs HETEROJEN (k=8, balanced) ---")
final %>%
  filter(dist == "negbin", k == 8, balance == "balanced") %>%
  select(n, var_type, power_anova, power_welch, power_kw) %>%
  print()

# Eksik deger kontrolu
na_check <- final %>%
  summarise(across(everything(), ~sum(is.na(.)))) %>%
  tidyr::pivot_longer(everything(),
                      names_to  = "sutun",
                      values_to = "na_sayisi") %>%
  filter(na_sayisi > 0)

if (nrow(na_check) == 0) {
  message("\nEksik deger yok.")
} else {
  message("\nEksik deger iceren sutunlar:")
  print(na_check)
}

# ========= KAYDET =========
write.csv(final, "H1_FINAL_MERGED2.csv", row.names = FALSE)
message("\nKaydedildi: H1_FINAL_MERGED2.csv")