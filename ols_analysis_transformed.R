# =============================================================================
# Advanced OLS Analysis with Transformations and Stepwise Selection
# =============================================================================
# This script performs:
#   1. Log transformations to improve OLS assumptions
#   2. Forward stepwise selection (building from 0 predictors)
#   3. Comprehensive diagnostics comparing linear vs log-log models
#   4. Model comparison and recommendation
#
# Response: gdp_per_capita (log-transformed)
# Predictors: Selected via forward stepwise based on R² improvement
# =============================================================================

# Load Required Packages ------------------------------------------------------
required_packages <- c("tidyverse", "car", "lmtest", "MASS", "ggplot2", "scales", "gridExtra")
missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(missing_packages) > 0) {
  cat("Installing missing packages:", paste(missing_packages, collapse = ", "), "\n")
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(car)
  library(lmtest)
  library(MASS)
  library(ggplot2)
  library(scales)
  library(gridExtra)
})

select <- dplyr::select

# Create output directory
output_dir <- "plots/ols_analysis_transformed"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cat("\n", strrep("=", 80), "\n")
cat("ADVANCED OLS ANALYSIS: TRANSFORMATIONS + STEPWISE SELECTION\n")
cat(strrep("=", 80), "\n\n")

# =============================================================================
# 1. DATA PREPARATION AND TRANSFORMATION
# =============================================================================

cat("Step 1: Loading and transforming data...\n")

# Load raw data
raw <- read.csv("owid-co2-data.csv", stringsAsFactors = FALSE)

# Apply same filtering as eda.R
data0 <- raw %>%
  filter(year >= 1851) %>%
  filter(iso_code != "", !is.na(iso_code)) %>%
  filter(!is.na(gdp), !is.na(population))

data1 <- data0 %>%
  mutate(
    iso_code = factor(iso_code),
    gdp_per_capita = gdp / population
  ) %>%
  drop_na(gdp_per_capita)

# Select the chosen variables
wanted_columns <- c(
  "iso_code",
  "year",
  "gdp_per_capita",
  "cement_co2_per_capita",
  "coal_co2_per_capita",
  "oil_co2_per_capita",
  "gas_co2_per_capita",
  "flaring_co2_per_capita",
  "land_use_change_co2_per_capita",
  "energy_per_capita",
  "methane_per_capita",
  "nitrous_oxide_per_capita"
)

data_selected <- data1 %>% select(all_of(wanted_columns))

# Get complete cases first
numeric_vars <- wanted_columns[!wanted_columns %in% c("iso_code")]
data_complete <- data_selected %>% 
  select(all_of(numeric_vars)) %>% 
  drop_na()

cat("  - Original data: ", nrow(data_selected), " rows\n")
cat("  - Complete cases: ", nrow(data_complete), " rows\n")

# Apply log transformations (log(x + 1) to handle zeros)
cat("\n  Applying log(x + 1) transformations...\n")

data_log <- data_complete %>%
  mutate(
    log_gdp = log(gdp_per_capita + 1),
    log_cement = log(cement_co2_per_capita + 1),
    log_coal = log(coal_co2_per_capita + 1),
    log_oil = log(oil_co2_per_capita + 1),
    log_gas = log(gas_co2_per_capita + 1),
    log_flaring = log(flaring_co2_per_capita + 1),
    log_land_use = log(land_use_change_co2_per_capita + 1),
    log_energy = log(energy_per_capita + 1),
    log_methane = log(methane_per_capita + 1),
    log_nitrous = log(nitrous_oxide_per_capita + 1)
  )

# Check for infinite or NA values after transformation
if (any(!is.finite(as.matrix(data_log %>% select(starts_with("log_")))))) {
  cat("  ⚠ Warning: Some log transformations produced non-finite values\n")
}

cat("  ✓ Transformations complete\n\n")

# =============================================================================
# 2. FORWARD STEPWISE SELECTION (Building from 0 predictors)
# =============================================================================

cat("Step 2: Forward stepwise selection starting from null model...\n")

# Available log-transformed predictors (excluding response)
available_predictors <- c(
  "year",
  "log_cement",
  "log_coal",
  "log_oil",
  "log_gas",
  "log_flaring",
  "log_land_use",
  "log_energy",
  "log_methane",
  "log_nitrous"
)

# Initialize
selected_predictors <- character(0)
remaining_predictors <- available_predictors
r_squared_history <- data.frame(
  Step = integer(),
  Predictor_Added = character(),
  R_squared = numeric(),
  Adj_R_squared = numeric(),
  Delta_R_squared = numeric(),
  stringsAsFactors = FALSE
)

# Null model (intercept only)
null_model <- lm(log_gdp ~ 1, data = data_log)
current_r2 <- 0
current_adj_r2 <- 0

cat(sprintf("  Starting R² = %.4f (null model)\n\n", current_r2))

# Threshold for stopping (when R² improvement < 0.005, i.e., < 0.5%)
r2_threshold <- 0.005
step <- 0

while (length(remaining_predictors) > 0) {
  step <- step + 1
  cat(sprintf("  --- Step %d ---\n", step))
  
  # Try adding each remaining predictor
  best_r2 <- current_r2
  best_predictor <- NULL
  best_model <- NULL
  
  for (pred in remaining_predictors) {
    # Build formula with current + this predictor
    test_predictors <- c(selected_predictors, pred)
    formula_test <- as.formula(paste("log_gdp ~", paste(test_predictors, collapse = " + ")))
    model_test <- lm(formula_test, data = data_log)
    test_r2 <- summary(model_test)$r.squared
    
    if (test_r2 > best_r2) {
      best_r2 <- test_r2
      best_predictor <- pred
      best_model <- model_test
    }
  }
  
  # Check if improvement is significant enough
  delta_r2 <- best_r2 - current_r2
  
  if (is.null(best_predictor) || delta_r2 < r2_threshold) {
    cat(sprintf("  Stopping: R² improvement (%.4f) below threshold (%.4f)\n", delta_r2, r2_threshold))
    break
  }
  
  # Add the best predictor
  selected_predictors <- c(selected_predictors, best_predictor)
  remaining_predictors <- setdiff(remaining_predictors, best_predictor)
  
  adj_r2 <- summary(best_model)$adj.r.squared
  
  # Record history
  r_squared_history <- rbind(r_squared_history, data.frame(
    Step = step,
    Predictor_Added = best_predictor,
    R_squared = best_r2,
    Adj_R_squared = adj_r2,
    Delta_R_squared = delta_r2
  ))
  
  cat(sprintf("  Added: %s\n", best_predictor))
  cat(sprintf("  R² = %.4f (Δ = %.4f), Adj R² = %.4f\n\n", best_r2, delta_r2, adj_r2))
  
  current_r2 <- best_r2
  current_adj_r2 <- adj_r2
}

cat("\n  Forward selection complete!\n")
cat(sprintf("  Final model: %d predictors selected\n", length(selected_predictors)))
cat(sprintf("  Final R² = %.4f, Adj R² = %.4f\n\n", current_r2, current_adj_r2))

# Display selection history
cat("  Selection History:\n")
print(r_squared_history, row.names = FALSE)

# Save selection history
write.csv(r_squared_history, file.path(output_dir, "stepwise_selection_history.csv"), row.names = FALSE)

# Plot R² progression
p_r2 <- ggplot(r_squared_history, aes(x = Step)) +
  geom_line(aes(y = R_squared, color = "R²"), size = 1.2) +
  geom_point(aes(y = R_squared, color = "R²"), size = 3) +
  geom_line(aes(y = Adj_R_squared, color = "Adj R²"), size = 1.2, linetype = "dashed") +
  geom_point(aes(y = Adj_R_squared, color = "Adj R²"), size = 3) +
  geom_text(aes(y = R_squared, label = Predictor_Added), vjust = -0.5, size = 3, angle = 15) +
  scale_color_manual(values = c("R²" = "#3498DB", "Adj R²" = "#E74C3C")) +
  labs(
    title = "Forward Stepwise Selection: R² Progression",
    subtitle = paste("Final model:", length(selected_predictors), "predictors"),
    x = "Step",
    y = "R²",
    color = "Metric"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(output_dir, "stepwise_r2_progression.png"), p_r2, width = 10, height = 6, dpi = 150)

cat("\n")

# =============================================================================
# 3. FIT FINAL SELECTED MODEL
# =============================================================================

cat("Step 3: Fitting final stepwise-selected log-log model...\n")

if (length(selected_predictors) == 0) {
  cat("  ⚠ No predictors selected! Using null model.\n")
  final_model <- null_model
} else {
  formula_final <- as.formula(paste("log_gdp ~", paste(selected_predictors, collapse = " + ")))
  final_model <- lm(formula_final, data = data_log)
  
  cat("\n  Final Model Formula:\n")
  print(formula_final)
  cat("\n")
  print(summary(final_model))
  cat("\n")
}

# Save model summary
sink(file.path(output_dir, "final_model_summary.txt"))
cat("FINAL STEPWISE-SELECTED LOG-LOG MODEL\n")
cat(strrep("=", 80), "\n\n")
cat("Selected Predictors:\n")
cat(paste("  -", selected_predictors, collapse = "\n"))
cat("\n\n")
print(summary(final_model))
sink()

# =============================================================================
# 4. VIF ANALYSIS ON SELECTED MODEL
# =============================================================================

cat("\nStep 4: Checking multicollinearity (VIF) on selected model...\n")

if (length(selected_predictors) > 1) {
  vif_values <- vif(final_model)
  vif_df <- data.frame(
    Variable = names(vif_values),
    VIF = as.numeric(vif_values)
  ) %>% arrange(desc(VIF))
  
  cat("\n  VIF values:\n")
  print(vif_df, row.names = FALSE)
  
  max_vif <- max(vif_df$VIF)
  if (max_vif > 10) {
    cat(sprintf("\n  ⚠ WARNING: Max VIF = %.2f (high multicollinearity)\n", max_vif))
  } else if (max_vif > 5) {
    cat(sprintf("\n  ⚠ CAUTION: Max VIF = %.2f (moderate multicollinearity)\n", max_vif))
  } else {
    cat(sprintf("\n  ✓ Max VIF = %.2f (no multicollinearity concerns)\n", max_vif))
  }
  
  # VIF plot
  p_vif <- ggplot(vif_df, aes(x = reorder(Variable, VIF), y = VIF)) +
    geom_col(fill = ifelse(vif_df$VIF > 10, "#E74C3C", 
                           ifelse(vif_df$VIF > 5, "#F39C12", "#27AE60"))) +
    geom_hline(yintercept = 10, linetype = "dashed", color = "red", linewidth = 1) +
    geom_hline(yintercept = 5, linetype = "dashed", color = "orange", linewidth = 0.8) +
    coord_flip() +
    labs(title = "Variance Inflation Factors (VIF) - Transformed Model",
         subtitle = "Red: VIF = 10, Orange: VIF = 5",
         x = "Predictor", y = "VIF") +
    theme_minimal(base_size = 12)
  
  ggsave(file.path(output_dir, "vif_final_model.png"), p_vif, width = 8, height = 6, dpi = 150)
  write.csv(vif_df, file.path(output_dir, "vif_final_model.csv"), row.names = FALSE)
  
} else {
  cat("  Only one predictor - VIF not applicable\n")
}

# =============================================================================
# 5. ASSUMPTION CHECKS: HETEROSCEDASTICITY
# =============================================================================

cat("\nStep 5: Checking constant variance (homoscedasticity)...\n")

# Breusch-Pagan test
bp_test <- bptest(final_model)
cat(sprintf("\n  Breusch-Pagan Test:\n"))
cat(sprintf("    BP statistic = %.4f, p-value = %.4f\n", bp_test$statistic, bp_test$p.value))

if (bp_test$p.value < 0.05) {
  cat("    ⚠ WARNING: Evidence of heteroscedasticity (p < 0.05)\n")
} else {
  cat("    ✓ No evidence of heteroscedasticity (p >= 0.05)\n")
}

# Residuals vs Fitted plot
fitted_vals <- fitted(final_model)
residuals_vals <- residuals(final_model)

p_resid_fitted <- ggplot(data.frame(Fitted = fitted_vals, Residuals = residuals_vals),
                         aes(x = Fitted, y = Residuals)) +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
  geom_smooth(method = "loess", color = "blue", se = TRUE) +
  labs(title = "Residuals vs Fitted Values",
       subtitle = "Check for constant variance and linearity",
       x = "Fitted Values", y = "Residuals") +
  theme_minimal(base_size = 12)

ggsave(file.path(output_dir, "residuals_vs_fitted.png"), p_resid_fitted, width = 8, height = 6, dpi = 150)

# Scale-Location plot
sqrt_abs_resid <- sqrt(abs(scale(residuals_vals)))
p_scale_loc <- ggplot(data.frame(Fitted = fitted_vals, SqrtResid = sqrt_abs_resid),
                      aes(x = Fitted, y = SqrtResid)) +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_smooth(method = "loess", color = "red", se = TRUE) +
  labs(title = "Scale-Location Plot",
       subtitle = "Check for homoscedasticity",
       x = "Fitted Values", y = "√|Standardized Residuals|") +
  theme_minimal(base_size = 12)

ggsave(file.path(output_dir, "scale_location.png"), p_scale_loc, width = 8, height = 6, dpi = 150)

cat("  ✓ Variance diagnostics complete\n")

# =============================================================================
# 6. ASSUMPTION CHECKS: NORMALITY
# =============================================================================

cat("\nStep 6: Checking normality of residuals...\n")

# Shapiro-Wilk test (sample if n > 5000)
n <- length(residuals_vals)
if (n > 5000) {
  set.seed(42)
  sample_indices <- sample(1:n, 5000)
  sw_test <- shapiro.test(residuals_vals[sample_indices])
  cat(sprintf("\n  Shapiro-Wilk Test (n = 5000 sample):\n"))
} else {
  sw_test <- shapiro.test(residuals_vals)
  cat(sprintf("\n  Shapiro-Wilk Test:\n"))
}

cat(sprintf("    W statistic = %.4f, p-value = %.4e\n", sw_test$statistic, sw_test$p.value))

if (sw_test$p.value < 0.05) {
  cat("    ⚠ WARNING: Residuals may not be normally distributed (p < 0.05)\n")
  cat("    Note: With large samples, minor deviations can be significant.\n")
} else {
  cat("    ✓ Residuals appear normally distributed (p >= 0.05)\n")
}

# Q-Q plot
p_qq <- ggplot(data.frame(Sample = residuals_vals), aes(sample = Sample)) +
  stat_qq(alpha = 0.4, size = 1.5) +
  stat_qq_line(color = "red", linewidth = 1) +
  labs(title = "Normal Q-Q Plot",
       subtitle = "Check for normality of residuals",
       x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme_minimal(base_size = 12)

ggsave(file.path(output_dir, "qq_plot.png"), p_qq, width = 8, height = 6, dpi = 150)

# Histogram of residuals
p_hist <- ggplot(data.frame(Residuals = residuals_vals), aes(x = Residuals)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50, fill = "#3498DB", alpha = 0.7) +
  geom_density(color = "red", linewidth = 1) +
  labs(title = "Distribution of Residuals",
       subtitle = "Check for normality",
       x = "Residuals", y = "Density") +
  theme_minimal(base_size = 12)

ggsave(file.path(output_dir, "residuals_histogram.png"), p_hist, width = 8, height = 6, dpi = 150)

cat("  ✓ Normality diagnostics complete\n")

# =============================================================================
# 7. INFLUENCE DIAGNOSTICS
# =============================================================================

cat("\nStep 7: Checking for influential observations...\n")

# Calculate diagnostics
n_obs <- nrow(data_log)
p_params <- length(coef(final_model))

leverage <- hatvalues(final_model)
cooks_d <- cooks.distance(final_model)
dffits_vals <- dffits(final_model)
std_resid <- rstandard(final_model)

# Thresholds
leverage_threshold <- 2 * p_params / n_obs
cooks_threshold <- 4 / n_obs
dffits_threshold <- 2 * sqrt(p_params / n_obs)

# Count influential points
high_leverage <- sum(leverage > leverage_threshold)
high_cooks <- sum(cooks_d > cooks_threshold)
high_dffits <- sum(abs(dffits_vals) > dffits_threshold)

cat(sprintf("\n  Thresholds:\n"))
cat(sprintf("    Leverage: %.4f (2p/n)\n", leverage_threshold))
cat(sprintf("    Cook's D: %.4f (4/n)\n", cooks_threshold))
cat(sprintf("    DFFITS: %.4f (2√(p/n))\n", dffits_threshold))

cat(sprintf("\n  Influential points:\n"))
cat(sprintf("    High leverage: %d (%.2f%%)\n", high_leverage, 100*high_leverage/n_obs))
cat(sprintf("    High Cook's D: %d (%.2f%%)\n", high_cooks, 100*high_cooks/n_obs))
cat(sprintf("    High DFFITS: %d (%.2f%%)\n", high_dffits, 100*high_dffits/n_obs))

# Create influence dataframe
influence_df <- data.frame(
  Index = 1:n_obs,
  Leverage = leverage,
  Cooks_D = cooks_d,
  DFFITS = dffits_vals,
  Std_Residual = std_resid,
  High_Leverage = leverage > leverage_threshold,
  High_Cooks = cooks_d > cooks_threshold,
  High_DFFITS = abs(dffits_vals) > dffits_threshold
) %>% arrange(desc(Cooks_D))

# Save top influential points
write.csv(head(influence_df, 50), file.path(output_dir, "top_influential_points.csv"), row.names = FALSE)

# Cook's D plot
p_cooks <- ggplot(influence_df, aes(x = Index, y = Cooks_D)) +
  geom_segment(aes(xend = Index, yend = 0), alpha = 0.5) +
  geom_point(aes(color = High_Cooks), alpha = 0.6) +
  geom_hline(yintercept = cooks_threshold, linetype = "dashed", color = "red") +
  scale_color_manual(values = c("FALSE" = "gray40", "TRUE" = "red")) +
  labs(title = "Cook's Distance",
       subtitle = paste("Influential points:", high_cooks),
       x = "Observation Index", y = "Cook's Distance") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

ggsave(file.path(output_dir, "cooks_distance.png"), p_cooks, width = 10, height = 6, dpi = 150)

# Leverage plot
p_leverage <- ggplot(influence_df, aes(x = Index, y = Leverage)) +
  geom_segment(aes(xend = Index, yend = 0), alpha = 0.5) +
  geom_point(aes(color = High_Leverage), alpha = 0.6) +
  geom_hline(yintercept = leverage_threshold, linetype = "dashed", color = "red") +
  scale_color_manual(values = c("FALSE" = "gray40", "TRUE" = "red")) +
  labs(title = "Leverage Values",
       subtitle = paste("High leverage points:", high_leverage),
       x = "Observation Index", y = "Leverage") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

ggsave(file.path(output_dir, "leverage.png"), p_leverage, width = 10, height = 6, dpi = 150)

cat("  ✓ Influence diagnostics complete\n")

# =============================================================================
# 8. OUTLIER DETECTION
# =============================================================================

cat("\nStep 8: Identifying outliers...\n")

outliers <- influence_df %>% filter(abs(Std_Residual) > 3)
n_outliers <- nrow(outliers)

cat(sprintf("  Outliers (|std residual| > 3): %d (%.2f%%)\n", n_outliers, 100*n_outliers/n_obs))

if (n_outliers > 0) {
  cat("\n  Top outliers:\n")
  print(head(outliers %>% select(Index, Std_Residual, Leverage, Cooks_D) %>% arrange(desc(abs(Std_Residual))), 10), row.names = FALSE)
  write.csv(outliers, file.path(output_dir, "outliers.csv"), row.names = FALSE)
}

cat("  ✓ Outlier detection complete\n")

# =============================================================================
# 9. INFLUENTIAL POINT REMOVAL AND MODEL REFITTING
# =============================================================================

cat("\nStep 9: Removing influential points and refitting model...\n")

# Strategy: Remove points that are BOTH highly influential AND outliers
# This is a conservative approach to avoid removing too many observations

# Define removal criteria (more conservative than detection thresholds)
removal_criteria <- influence_df %>%
  filter(
    (High_Cooks & abs(Std_Residual) > 3) |  # High Cook's D AND outlier
    (abs(Std_Residual) > 4)                  # Extreme outliers (|resid| > 4)
  )

n_removed <- nrow(removal_criteria)
pct_removed <- 100 * n_removed / n_obs

cat(sprintf("  Points to remove: %d (%.2f%% of data)\n", n_removed, pct_removed))
cat("  Removal criteria:\n")
cat("    - High Cook's D AND |std resid| > 3, OR\n")
cat("    - Extreme outliers with |std resid| > 4\n\n")

if (n_removed > 0) {
  # Remove influential points
  data_log_clean <- data_log[-removal_criteria$Index, ]
  
  cat(sprintf("  Original data: %d observations\n", n_obs))
  cat(sprintf("  Cleaned data: %d observations\n", nrow(data_log_clean)))
  
  # Refit model with cleaned data
  cat("\n  Refitting model with cleaned data...\n")
  final_model_clean <- lm(formula_final, data = data_log_clean)
  
  cat("\n  Cleaned Model Summary:\n")
  print(summary(final_model_clean))
  cat("\n")
  
  # Compare original vs cleaned model
  cat("  MODEL COMPARISON (Before vs After Removal):\n\n")
  
  comparison_clean <- data.frame(
    Metric = c("N", "R²", "Adj R²", "Residual SE", "F-statistic", 
               "BP p-value", "SW p-value", "Max VIF", "Outliers (>3)", "High Cook's D"),
    Before_Removal = c(
      n_obs,
      summary(final_model)$r.squared,
      summary(final_model)$adj.r.squared,
      summary(final_model)$sigma,
      summary(final_model)$fstatistic[1],
      bp_test$p.value,
      sw_test$p.value,
      ifelse(length(selected_predictors) > 1, max(vif(final_model)), NA),
      n_outliers,
      high_cooks
    ),
    After_Removal = c(
      nrow(data_log_clean),
      summary(final_model_clean)$r.squared,
      summary(final_model_clean)$adj.r.squared,
      summary(final_model_clean)$sigma,
      summary(final_model_clean)$fstatistic[1],
      bptest(final_model_clean)$p.value,
      ifelse(nrow(data_log_clean) > 5000,
             shapiro.test(sample(residuals(final_model_clean), 5000))$p.value,
             shapiro.test(residuals(final_model_clean))$p.value),
      ifelse(length(selected_predictors) > 1, max(vif(final_model_clean)), NA),
      sum(abs(rstandard(final_model_clean)) > 3),
      sum(cooks.distance(final_model_clean) > (4/nrow(data_log_clean)))
    )
  )
  
  print(comparison_clean, row.names = FALSE, digits = 4)
  
  # Calculate improvements
  r2_improvement <- comparison_clean$After_Removal[2] - comparison_clean$Before_Removal[2]
  se_improvement <- comparison_clean$Before_Removal[4] - comparison_clean$After_Removal[4]
  bp_improvement <- comparison_clean$After_Removal[6] / comparison_clean$Before_Removal[6]
  sw_improvement <- comparison_clean$After_Removal[7] / comparison_clean$Before_Removal[7]
  
  cat("\n  IMPROVEMENTS:\n")
  cat(sprintf("    R² change: %+.4f (%.2f%%)\n", r2_improvement, 100*r2_improvement/comparison_clean$Before_Removal[2]))
  cat(sprintf("    Residual SE reduction: %.4f (%.2f%%)\n", se_improvement, 100*se_improvement/comparison_clean$Before_Removal[4]))
  cat(sprintf("    BP p-value multiplier: %.2fx %s\n", bp_improvement, 
              ifelse(bp_improvement > 1, "✓ (better)", "(worse)")))
  cat(sprintf("    SW p-value multiplier: %.2fx %s\n", sw_improvement,
              ifelse(sw_improvement > 1, "✓ (better)", "(worse)")))
  cat(sprintf("    Outliers reduced: %d → %d\n", 
              as.integer(comparison_clean$Before_Removal[9]), 
              as.integer(comparison_clean$After_Removal[9])))
  cat(sprintf("    High Cook's D reduced: %d → %d\n",
              as.integer(comparison_clean$Before_Removal[10]),
              as.integer(comparison_clean$After_Removal[10])))
  
  # Update final model to cleaned version
  cat("\n  ✓ Using CLEANED model for remaining diagnostics\n")
  final_model <- final_model_clean
  data_log <- data_log_clean
  n_obs <- nrow(data_log_clean)
  
  # Recalculate diagnostics for cleaned model
  leverage <- hatvalues(final_model)
  cooks_d <- cooks.distance(final_model)
  dffits_vals <- dffits(final_model)
  std_resid <- rstandard(final_model)
  fitted_vals <- fitted(final_model)
  residuals_vals <- residuals(final_model)
  
  # Update thresholds
  p_params <- length(coef(final_model))
  leverage_threshold <- 2 * p_params / n_obs
  cooks_threshold <- 4 / n_obs
  dffits_threshold <- 2 * sqrt(p_params / n_obs)
  
  # Save cleaned model summary
  sink(file.path(output_dir, "final_model_summary.txt"))
  cat("FINAL STEPWISE-SELECTED LOG-LOG MODEL (AFTER INFLUENTIAL POINT REMOVAL)\n")
  cat(strrep("=", 80), "\n\n")
  cat(sprintf("Removed %d influential points (%.2f%% of original data)\n", n_removed, pct_removed))
  cat("\nRemoval criteria:\n")
  cat("  - High Cook's D AND |std residual| > 3, OR\n")
  cat("  - Extreme outliers with |std residual| > 4\n\n")
  cat("Selected Predictors:\n")
  cat(paste("  -", selected_predictors, collapse = "\n"))
  cat("\n\n")
  print(summary(final_model))
  sink()
  
} else {
  cat("  No points meet removal criteria - keeping original model\n")
}

cat("\n")

# =============================================================================
# 10. COMPREHENSIVE DIAGNOSTIC PLOTS (CLEANED MODEL)
# =============================================================================

cat("\nStep 10: Creating diagnostic plot panels for cleaned model...\n")

# Standard 4-panel diagnostic plot
png(file.path(output_dir, "diagnostic_panel_4.png"), width = 1200, height = 1000, res = 150)
par(mfrow = c(2, 2))
plot(final_model)
dev.off()

# Custom 6-panel plot
png(file.path(output_dir, "diagnostic_panel_6.png"), width = 1800, height = 1200, res = 150)
par(mfrow = c(2, 3))
plot(final_model, which = 1:6)
dev.off()

cat("  ✓ Diagnostic panels created\n")

# =============================================================================
# 11. MODEL COMPARISON: LINEAR VS LOG-LOG
# =============================================================================

cat("\nStep 11: Comparing linear vs log-log models (both with cleaned data)...\n")

# Fit linear model with same selected predictors (use original variables)
if (length(selected_predictors) > 0) {
  # Map log predictors back to original names
  original_predictors <- gsub("^log_", "", selected_predictors)
  original_predictors <- gsub("^cement$", "cement_co2_per_capita", original_predictors)
  original_predictors <- gsub("^coal$", "coal_co2_per_capita", original_predictors)
  original_predictors <- gsub("^oil$", "oil_co2_per_capita", original_predictors)
  original_predictors <- gsub("^gas$", "gas_co2_per_capita", original_predictors)
  original_predictors <- gsub("^flaring$", "flaring_co2_per_capita", original_predictors)
  original_predictors <- gsub("^land_use$", "land_use_change_co2_per_capita", original_predictors)
  original_predictors <- gsub("^energy$", "energy_per_capita", original_predictors)
  original_predictors <- gsub("^methane$", "methane_per_capita", original_predictors)
  original_predictors <- gsub("^nitrous$", "nitrous_oxide_per_capita", original_predictors)
  
  # Need to align data_complete with cleaned indices
  # Get the row indices that were kept in data_log
  original_indices <- as.numeric(rownames(data_log))
  data_complete_clean <- data_complete[original_indices, ]
  
  formula_linear <- as.formula(paste("gdp_per_capita ~", paste(original_predictors, collapse = " + ")))
  model_linear <- lm(formula_linear, data = data_complete_clean)
  
  # Compare diagnostics
  comparison <- data.frame(
    Metric = c("N", "R²", "Adj R²", "Residual SE", "AIC", "BIC", 
               "BP Test p-value", "SW Test p-value", "Max VIF"),
    Linear_Model = c(
      nrow(data_complete_clean),
      summary(model_linear)$r.squared,
      summary(model_linear)$adj.r.squared,
      summary(model_linear)$sigma,
      AIC(model_linear),
      BIC(model_linear),
      bptest(model_linear)$p.value,
      ifelse(nrow(data_complete_clean) > 5000, 
             shapiro.test(sample(residuals(model_linear), 5000))$p.value,
             shapiro.test(residuals(model_linear))$p.value),
      ifelse(length(original_predictors) > 1, max(vif(model_linear)), NA)
    ),
    LogLog_Model = c(
      nrow(data_log),
      summary(final_model)$r.squared,
      summary(final_model)$adj.r.squared,
      summary(final_model)$sigma,
      AIC(final_model),
      BIC(final_model),
      bptest(final_model)$p.value,
      ifelse(nrow(data_log) > 5000,
             shapiro.test(sample(residuals(final_model), 5000))$p.value,
             shapiro.test(residuals(final_model))$p.value),
      ifelse(length(selected_predictors) > 1, max(vif(final_model)), NA)
    )
  )
  
  cat("\n  Model Comparison:\n")
  print(comparison, row.names = FALSE, digits = 4)
  
  write.csv(comparison, file.path(output_dir, "model_comparison.csv"), row.names = FALSE)
  
  # Recommendation
  cat("\n  RECOMMENDATION:\n")
  if (comparison$LogLog_Model[7] > comparison$Linear_Model[7]) {
    cat("    ✓ Log-log model shows BETTER homoscedasticity (higher BP p-value)\n")
  } else {
    cat("    Linear model shows better homoscedasticity\n")
  }
  
  if (comparison$LogLog_Model[8] > comparison$Linear_Model[8]) {
    cat("    ✓ Log-log model shows BETTER normality (higher SW p-value)\n")
  } else {
    cat("    Linear model shows better normality\n")
  }
  
  if (!is.na(comparison$LogLog_Model[9]) && !is.na(comparison$Linear_Model[9])) {
    if (comparison$LogLog_Model[9] < comparison$Linear_Model[9]) {
      cat("    ✓ Log-log model shows LOWER multicollinearity (lower max VIF)\n")
    } else {
      cat("    Linear model shows lower multicollinearity\n")
    }
  }
  
  if (comparison$LogLog_Model[5] < comparison$Linear_Model[5]) {
    cat("    ✓ Log-log model has LOWER AIC (better fit)\n")
  } else {
    cat("    Linear model has lower AIC\n")
  }
}

# =============================================================================
# FINAL SUMMARY
# =============================================================================

cat("\n", strrep("=", 80), "\n")
cat("ANALYSIS SUMMARY\n")
cat(strrep("=", 80), "\n\n")

cat("STEPWISE SELECTION RESULTS:\n")
cat(sprintf("  Predictors selected: %d out of %d available\n", length(selected_predictors), length(available_predictors)))
cat(sprintf("  Final R²: %.4f\n", current_r2))
cat(sprintf("  Final Adj R²: %.4f\n", current_adj_r2))
cat("\n  Selected predictors:\n")
for (i in seq_along(selected_predictors)) {
  cat(sprintf("    %d. %s\n", i, selected_predictors[i]))
}

cat("\n\nASSUMPTION CHECKS (LOG-LOG MODEL - AFTER CLEANING):\n")

# Recalculate final diagnostics
bp_test_final <- bptest(final_model)
sw_test_final <- ifelse(n_obs > 5000,
                        shapiro.test(sample(residuals(final_model), 5000))$p.value,
                        shapiro.test(residuals(final_model))$p.value)
high_cooks_final <- sum(cooks_d > cooks_threshold)
n_outliers_final <- sum(abs(std_resid) > 3)

cat(sprintf("  1. Homoscedasticity: BP p-value = %.4e %s\n", 
            bp_test_final$p.value, 
            ifelse(bp_test_final$p.value >= 0.05, "✓", "⚠")))
cat(sprintf("  2. Normality: SW p-value = %.4e %s\n", 
            sw_test_final,
            ifelse(sw_test_final >= 0.05, "✓", "⚠")))
if (length(selected_predictors) > 1) {
  max_vif_val <- max(vif(final_model))
  cat(sprintf("  3. Multicollinearity: Max VIF = %.2f %s\n", 
              max_vif_val,
              ifelse(max_vif_val < 5, "✓", ifelse(max_vif_val < 10, "⚠", "✗"))))
}
cat(sprintf("  4. Influential points: %.2f%% with high Cook's D %s\n",
            100*high_cooks_final/n_obs,
            ifelse(high_cooks_final/n_obs < 0.05, "✓", "⚠")))
cat(sprintf("  5. Outliers: %.2f%% with |std resid| > 3 %s\n",
            100*n_outliers_final/n_obs,
            ifelse(n_outliers_final/n_obs < 0.05, "✓", "⚠")))

cat("\n\nOUTPUT FILES:\n")
cat(sprintf("  All results saved to: %s\n", output_dir))
cat("  - stepwise_selection_history.csv (selection process)\n")
cat("  - stepwise_r2_progression.png (R² plot)\n")
cat("  - final_model_summary.txt (model details)\n")
cat("  - model_comparison.csv (linear vs log-log)\n")
cat("  - vif_final_model.csv/png (multicollinearity)\n")
cat("  - Diagnostic plots (residuals, QQ, leverage, Cook's D, etc.)\n")
cat("  - Influence and outlier CSVs\n")

cat("\n", strrep("=", 80), "\n")
cat("ANALYSIS COMPLETE!\n")
cat(strrep("=", 80), "\n\n")
