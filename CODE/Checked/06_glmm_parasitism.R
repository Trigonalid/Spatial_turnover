# ============================================================
# 06_glmm_parasitism.R   >>> ANALYSIS ONLY <<<
# Study: Approaching the limit of diversity: spatial turnover
#        amplifies with trophic level (Libra et al.)
#
# Analysis 6 (core GLMM). Tests whether parasitism rate on a caterpillar
# species at a locality increases with the number of parasitoid species
# attacking it (reviewer comment OL18.1).
#
# Unit: locality x caterpillar species. Binomial GLMM:
#   cbind(N_parasitized, N_total - N_parasitized) ~
#       N_par_species + log(N_total) + (1|locality) + (1|CAT_sp)
# log(N_total) controls for sampling effort. Main threshold N_total >= 10;
# robustness at >= 5 and >= 20.
#
# Outputs: tables (Word + csv) + output/rds/ for the plotting script.
# Figures (A/B/C variants) are built in 06_glmm_parasitism_plots.R.
# ============================================================

# ------------------------------------------------------------
# 0. Libraries
# ------------------------------------------------------------
library(dplyr)
library(tidyr)
library(tibble)
library(readxl)
library(here)
library(lme4)        # glmer()
library(broom.mixed) # tidy()
library(flextable)
library(officer)

# ------------------------------------------------------------
# 1. Load data (robust load; Ohu2 excluded)
# PAR_sp MUST be read correctly here: N_parasitized and N_par_species
# both derive from it, so the guess_max / as.character fix is essential.
# ------------------------------------------------------------
MASTER <- read_excel(here("DATA/MASTER.xlsx"), guess_max = 1048576) %>%
  as_tibble() %>%
  rename(
    PLANT_sp = PLANT_species_code,
    CAT_sp   = CAT_scientific_name,
    PAR_sp   = PAR_species_code
  ) %>%
  mutate(PAR_sp = na_if(as.character(PAR_sp), "")) %>%
  filter(locality != "Ohu2")

dir.create("output",     showWarnings = FALSE, recursive = TRUE)
dir.create("output/rds", showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# 2. Aggregate per (locality x CAT_sp)
# Each MASTER row = one reared caterpillar; PAR_sp non-NA = parasitized.
# ------------------------------------------------------------
agg_data <- MASTER %>%
  filter(!is.na(CAT_sp)) %>%
  group_by(locality, CAT_sp) %>%
  summarise(
    N_total       = n(),
    N_parasitized = sum(!is.na(PAR_sp)),
    N_par_species = n_distinct(PAR_sp[!is.na(PAR_sp)]),
    .groups = "drop"
  ) %>%
  mutate(parasitism_rate = N_parasitized / N_total)

cat("\n--- Aggregated data summary ---\n")
cat(sprintf("  Total (locality x CAT_sp) combinations: %d\n", nrow(agg_data)))
cat(sprintf("  With >= 5 individuals:  %d\n", sum(agg_data$N_total >= 5)))
cat(sprintf("  With >= 10 individuals: %d\n", sum(agg_data$N_total >= 10)))
cat(sprintf("  With >= 20 individuals: %d\n\n", sum(agg_data$N_total >= 20)))

# ------------------------------------------------------------
# 3. Descriptive summary by N_par_species category (N_total >= 10)
# ------------------------------------------------------------
agg_main <- agg_data %>% filter(N_total >= 10)

descriptive_summary <- agg_main %>%
  mutate(n_par_cat = factor(case_when(
    N_par_species == 0 ~ "0 parasitoids",
    N_par_species == 1 ~ "1 parasitoid",
    N_par_species == 2 ~ "2 parasitoids",
    N_par_species >= 3 ~ "3+ parasitoids"),
    levels = c("0 parasitoids", "1 parasitoid", "2 parasitoids", "3+ parasitoids"))) %>%
  group_by(n_par_cat) %>%
  summarise(
    n_combinations    = n(),
    mean_parasitism   = mean(parasitism_rate, na.rm = TRUE),
    sd_parasitism     = sd(parasitism_rate,   na.rm = TRUE),
    median_parasitism = median(parasitism_rate, na.rm = TRUE),
    mean_N_total      = mean(N_total),
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

cat("--- Parasitism rate by number of attacking parasitoid species (N_total >= 10) ---\n")
print(descriptive_summary); cat("\n")

# ------------------------------------------------------------
# 4. GLMM fit (binomial; random intercepts locality + CAT_sp)
# ------------------------------------------------------------
fit_glmm <- function(df, label) {
  df_fit <- df %>% mutate(N_unparasitized = N_total - N_parasitized, logN = log(N_total))
  form <- cbind(N_parasitized, N_unparasitized) ~ N_par_species + logN +
    (1 | locality) + (1 | CAT_sp)
  ctrl <- lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))

  m <- tryCatch(
    lme4::glmer(form, data = df_fit, family = binomial("logit"), control = ctrl),
    error   = function(e) NULL,
    warning = function(w) suppressWarnings(
      lme4::glmer(form, data = df_fit, family = binomial("logit"), control = ctrl)))

  if (is.null(m)) { cat(sprintf("[GLMM %s] Model failed to fit.\n", label)); return(NULL) }

  coef <- broom.mixed::tidy(m, effects = "fixed") %>%
    mutate(threshold = label, across(where(is.numeric), ~ round(.x, 4)))
  list(model = m, coef = coef)
}

res_main    <- fit_glmm(agg_main,                              "N_total >= 10 (main)")
res_liberal <- fit_glmm(agg_data %>% filter(N_total >= 5),     "N_total >= 5 (liberal)")
res_strict  <- fit_glmm(agg_data %>% filter(N_total >= 20),    "N_total >= 20 (strict)")

glmm_results <- bind_rows(res_main$coef, res_liberal$coef, res_strict$coef) %>%
  select(threshold, term, estimate, std.error, statistic, p.value)

cat("--- GLMM results ---\n"); print(glmm_results); cat("\n")

# ------------------------------------------------------------
# 5. Sample sizes by threshold
# ------------------------------------------------------------
sample_at <- function(thr) {
  agg_data %>% filter(N_total >= thr) %>%
    summarise(threshold = paste0("N_total >= ", thr),
              n_combinations = n(), n_localities = n_distinct(locality),
              n_caterpillars = n_distinct(CAT_sp),
              total_reared = sum(N_total), total_parasit = sum(N_parasitized))
}
sample_summary <- bind_rows(sample_at(5), sample_at(10), sample_at(20))
cat("--- Sample sizes by threshold ---\n"); print(sample_summary); cat("\n")

# ------------------------------------------------------------
# 6. Plotting support: main-model predictions (fixed effects only)
# ------------------------------------------------------------
plot_df <- agg_main %>% mutate(N_unparasitized = N_total - N_parasitized, logN = log(N_total))

pred_df <- tibble(
  N_par_species = seq(0, max(plot_df$N_par_species), by = 1),
  logN = mean(plot_df$logN), locality = NA_character_, CAT_sp = NA_character_)
if (!is.null(res_main$model)) {
  pred_df$pred <- predict(res_main$model, newdata = pred_df, re.form = NA, type = "response")
}

saveRDS(plot_df,      "output/rds/glmm_plot_df.rds")
saveRDS(pred_df,      "output/rds/glmm_pred_df.rds")
saveRDS(glmm_results, "output/rds/glmm_results.rds")

# ------------------------------------------------------------
# 7. Export tables (csv + Word) — Table S11
# ------------------------------------------------------------
write.csv(agg_data,            "output/rds/glmm_aggregated_data.csv",  row.names = FALSE)
write.csv(sample_summary,      "output/rds/Table_S11a_sample_sizes.csv",   row.names = FALSE)
write.csv(descriptive_summary, "output/rds/Table_S11b_descriptive.csv",    row.names = FALSE)
write.csv(glmm_results,        "output/rds/Table_S11c_glmm_results.csv",   row.names = FALSE)

doc_pd <- officer::read_docx() %>%
  officer::body_add_par("Table S11: Local parasitism rate increases with the number of attacking parasitoid species",
                        style = "heading 1") %>%
  officer::body_add_par("a) Sample sizes across thresholds", style = "heading 2") %>%
  flextable::body_add_flextable(flextable(sample_summary)) %>%
  officer::body_add_par("b) Parasitism rate by N_par_species (N_total >= 10)", style = "heading 2") %>%
  flextable::body_add_flextable(flextable(descriptive_summary)) %>%
  officer::body_add_par("c) GLMM: parasitism ~ N_par_species + log(N_total) + (1|locality) + (1|CAT_sp)",
                        style = "heading 2") %>%
  flextable::body_add_flextable(flextable(glmm_results))
print(doc_pd, target = "output/Table_S11_parasitism_diversity.docx")

# ------------------------------------------------------------
# 8. Console interpretation
# ------------------------------------------------------------
cat("=== PARASITISM RATE vs PARASITOID DIVERSITY ===\n")
cat("  N_par_species positive & significant -> more parasitoid species means\n")
cat("    higher parasitism (supports 'limited top-down control').\n")
cat("  Non-significant / near zero -> compensatory effects; soften the claim.\n")
