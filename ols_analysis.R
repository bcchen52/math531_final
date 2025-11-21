# =============================================================================
# Comprehensive OLS Analysis for GDP Per Capita Prediction
# =============================================================================
# This script performs thorough Ordinary Least Squares regression analysis
# with complete diagnostic checking as per MA 455/531 requirements.
#
# Analysis includes:
#   1. Data preparation and correlation analysis
#   2. Initial full model fit
#   3. VIF calculation and multicollinearity assessment
#   4. Iterative VIF reduction to handle collinearity
#   5. Residual diagnostics (constant variance, normality)
#   6. Leverage and influence diagnostics (Cook's D, DFFITS, DFBETAS)
#   7. Outlier detection
#   8. Comprehensive diagnostic plots
#
# Response: gdp_per_capita
# Predictors: 10 variables (year + per-capita emission metrics)
# =============================================================================

# Load Required Packages ------------------------------------------------------
# Install missing packages if needed
required_packages <- c("tidyverse", "car", "lmtest", "MASS", "ggplot2", "scales", "gridExtra")
missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(missing_packages) > 0) {
  cat("Installing missing packages:", paste(missing_packages, collapse = ", "), "\n")
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(car)        # For VIF, residual plots
  library(lmtest)     # For Breusch-Pagan test
  library(MASS)       # For robust regression comparison
  library(ggplot2)
  library(scales)
  library(gridExtra)  # For arranging multiple plots
})

# Ensure dplyr::select is used (avoid MASS::select conflict)
select <- dplyr::select

# Create output directory
output_dir <- "plots/ols_analysis"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cat("\n", strrep("=", 80), "\n")
cat("COMPREHENSIVE OLS REGRESSION ANALYSIS\n")
cat("Response: gdp_per_capita\n")
cat(strrep("=", 80), "\n\n")

# =============================================================================
# 1. DATA PREPARATION
# =============================================================================

cat("Step 1: Loading and preparing data...\n")

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

# Select the 11 chosen variables (gdp_per_capita is response, removed co2_per_capita)
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

# For regression, use complete cases only
numeric_vars <- wanted_columns[!wanted_columns %in% c("iso_code")]
data_complete <- data_selected %>% 
  select(all_of(numeric_vars)) %>% 
  drop_na()

cat("  - Original data: ", nrow(data_selected), " rows\n")
cat("  - Complete cases: ", nrow(data_complete), " rows\n")
cat("  - Variables: ", ncol(data_complete), " (1 response + ", ncol(data_complete)-1, " predictors)\n\n")

# =============================================================================
# 2. CORRELATION ANALYSIS
# =============================================================================

cat("Step 2: Analyzing correlations between predictors...\n")

# Compute correlation matrix for predictors only (exclude response gdp_per_capita)
predictors <- setdiff(names(data_complete), "gdp_per_capita")
cor_mat <- cor(data_complete[, predictors])

# Find highly correlated pairs (|r| > 0.7)
high_cor_pairs <- which(abs(cor_mat) > 0.7 & abs(cor_mat) < 1, arr.ind = TRUE)
if (nrow(high_cor_pairs) > 0) {
  cat("\n  High correlations (|r| > 0.7) between predictors:\n")
  for (i in 1:nrow(high_cor_pairs)) {
    row_idx <- high_cor_pairs[i, 1]
    col_idx <- high_cor_pairs[i, 2]
    if (row_idx < col_idx) {  # Avoid duplicates
      var1 <- rownames(cor_mat)[row_idx]
      var2 <- colnames(cor_mat)[col_idx]
      corr_val <- cor_mat[row_idx, col_idx]
      cat(sprintf("    %s <-> %s: %.3f\n", var1, var2, corr_val))
    }
  }
}

# Save correlation matrix
write.csv(cor_mat, file.path(output_dir, "predictor_correlations.csv"))

# Correlation heatmap
cor_long <- as.data.frame(cor_mat) %>%
  mutate(var1 = rownames(cor_mat)) %>%
  pivot_longer(-var1, names_to = "var2", values_to = "corr")

p_cor <- ggplot(cor_long, aes(var1, var2, fill = corr)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", corr)), size = 2.5) +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1)) +
  labs(title = "Predictor Correlation Matrix",
       subtitle = "Values show Pearson correlation coefficients",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave(file.path(output_dir, "predictor_correlation_heatmap.png"), p_cor,
       width = 10, height = 8, dpi = 150)

cat("  ✓ Correlation analysis complete\n\n")

# =============================================================================
# 3. INITIAL FULL MODEL
# =============================================================================

cat("Step 3: Fitting initial full model with all predictors...\n")

# Fit full model
formula_full <- as.formula(paste("gdp_per_capita ~", paste(predictors, collapse = " + ")))
model_full <- lm(formula_full, data = data_complete)

cat("\n")
print(summary(model_full))
cat("\n")

# Save full model summary
sink(file.path(output_dir, "model_full_summary.txt"))
print(summary(model_full))
sink()

# =============================================================================
# 4. VIF ANALYSIS AND MULTICOLLINEARITY
# =============================================================================

cat("\nStep 4: Checking for multicollinearity using VIF...\n")

# Calculate VIF for all predictors
vif_values <- vif(model_full)
vif_df <- data.frame(
  Variable = names(vif_values),
  VIF = as.numeric(vif_values)
) %>% arrange(desc(VIF))

cat("\n  Variance Inflation Factors:\n")
print(vif_df, row.names = FALSE)

# Identify problematic VIFs (typically VIF > 10 indicates serious multicollinearity)
high_vif <- vif_df %>% filter(VIF > 10)
if (nrow(high_vif) > 0) {
  cat("\n  ⚠ WARNING: Variables with VIF > 10 (high multicollinearity):\n")
  print(high_vif, row.names = FALSE)
}

# VIF plot
p_vif <- ggplot(vif_df, aes(x = reorder(Variable, VIF), y = VIF)) +
  geom_col(fill = ifelse(vif_df$VIF > 10, "#E74C3C", 
                         ifelse(vif_df$VIF > 5, "#F39C12", "#27AE60"))) +
  geom_hline(yintercept = 10, linetype = "dashed", color = "red", size = 1) +
  geom_hline(yintercept = 5, linetype = "dashed", color = "orange", size = 0.8) +
  coord_flip() +
  labs(title = "Variance Inflation Factors (VIF)",
       subtitle = "Red line: VIF = 10 (high concern), Orange line: VIF = 5 (moderate concern)",
       x = "Predictor", y = "VIF") +
  theme_minimal(base_size = 12)

ggsave(file.path(output_dir, "vif_initial.png"), p_vif, width = 8, height = 6, dpi = 150)

write.csv(vif_df, file.path(output_dir, "vif_initial.csv"), row.names = FALSE)

# =============================================================================
# 5. ITERATIVE VIF REDUCTION (if needed)
# =============================================================================

cat("\nStep 5: Iterative VIF reduction to handle multicollinearity...\n")

current_predictors <- predictors
iteration <- 0
max_iterations <- 10

while (TRUE) {
  iteration <- iteration + 1
  
  # Fit model with current predictors
  formula_current <- as.formula(paste("gdp_per_capita ~", 
                                      paste(current_predictors, collapse = " + ")))
  model_current <- lm(formula_current, data = data_complete)
  
  # Calculate VIF
  vif_current <- vif(model_current)
  max_vif <- max(vif_current)
  
  cat(sprintf("  Iteration %d: Max VIF = %.2f\n", iteration, max_vif))
  
  # Stop if all VIFs are acceptable or max iterations reached
  if (max_vif < 5 || iteration >= max_iterations) {
    break
  }
  
  # Remove predictor with highest VIF
  var_to_remove <- names(which.max(vif_current))
  cat(sprintf("    Removing: %s (VIF = %.2f)\n", var_to_remove, max_vif))
  current_predictors <- setdiff(current_predictors, var_to_remove)
}

# Final reduced model
model_reduced <- model_current
vif_reduced <- vif(model_reduced)

cat("\n  Final Model Predictors:\n")
cat("   ", paste(current_predictors, collapse = ", "), "\n\n")

cat("  Final VIF values:\n")
vif_reduced_df <- data.frame(
  Variable = names(vif_reduced),
  VIF = as.numeric(vif_reduced)
) %>% arrange(desc(VIF))
print(vif_reduced_df, row.names = FALSE)

# Save reduced model VIF
write.csv(vif_reduced_df, file.path(output_dir, "vif_reduced.csv"), row.names = FALSE)

cat("\n")
print(summary(model_reduced))

# Save reduced model summary
sink(file.path(output_dir, "model_reduced_summary.txt"))
print(summary(model_reduced))
sink()

# Use reduced model for all subsequent diagnostics
model_final <- model_reduced

cat("\n  ✓ Model selection complete. Using reduced model for diagnostics.\n\n")

# =============================================================================
# 6. RESIDUAL DIAGNOSTICS - CONSTANT VARIANCE
# =============================================================================

cat("Step 6: Checking constant variance assumption...\n")

# Breusch-Pagan test for heteroscedasticity
bp_test <- bptest(model_final)
cat("\n  Breusch-Pagan Test for Heteroscedasticity:\n")
cat(sprintf("    BP statistic = %.4f, p-value = %.4f\n", bp_test$statistic, bp_test$p.value))
if (bp_test$p.value < 0.05) {
  cat("    ⚠ WARNING: Evidence of heteroscedasticity (p < 0.05)\n")
} else {
  cat("    ✓ No significant evidence of heteroscedasticity (p >= 0.05)\n")
}

# Residuals vs Fitted plot
residuals_std <- rstandard(model_final)
fitted_vals <- fitted(model_final)

p_resid_fitted <- ggplot(data.frame(fitted = fitted_vals, residual = residuals_std),
                         aes(fitted, residual)) +
  geom_point(alpha = 0.5, color = "#3498DB") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_smooth(se = TRUE, color = "#E74C3C", method = "loess") +
  labs(title = "Residuals vs Fitted Values",
       subtitle = "Checking for constant variance and linearity",
       x = "Fitted Values", y = "Standardized Residuals") +
  theme_minimal(base_size = 12)

ggsave(file.path(output_dir, "residuals_vs_fitted.png"), p_resid_fitted,
       width = 8, height = 6, dpi = 150)

# Scale-Location plot (sqrt of standardized residuals)
p_scale_location <- ggplot(data.frame(fitted = fitted_vals, 
                                      sqrt_std_resid = sqrt(abs(residuals_std))),
                           aes(fitted, sqrt_std_resid)) +
  geom_point(alpha = 0.5, color = "#9B59B6") +
  geom_smooth(se = TRUE, color = "#E74C3C", method = "loess") +
  labs(title = "Scale-Location Plot",
       subtitle = "Checking homoscedasticity (constant variance)",
       x = "Fitted Values", y = "√|Standardized Residuals|") +
  theme_minimal(base_size = 12)

ggsave(file.path(output_dir, "scale_location.png"), p_scale_location,
       width = 8, height = 6, dpi = 150)

cat("  ✓ Variance diagnostics complete\n\n")

# =============================================================================
# 7. NORMALITY OF RESIDUALS
# =============================================================================

cat("Step 7: Checking normality assumption for residuals...\n")

# Shapiro-Wilk test (use sample if n > 5000)
n_resid <- length(residuals_std)
if (n_resid > 5000) {
  set.seed(42)
  sample_indices <- sample(1:n_resid, 5000)
  sw_test <- shapiro.test(residuals_std[sample_indices])
  cat("\n  Shapiro-Wilk Test (on sample of 5000):\n")
} else {
  sw_test <- shapiro.test(residuals_std)
  cat("\n  Shapiro-Wilk Test:\n")
}

cat(sprintf("    W statistic = %.4f, p-value = %.4e\n", sw_test$statistic, sw_test$p.value))
if (sw_test$p.value < 0.05) {
  cat("    ⚠ WARNING: Residuals may not be normally distributed (p < 0.05)\n")
  cat("    Note: With large samples, minor deviations can be significant.\n")
} else {
  cat("    ✓ Residuals appear normally distributed (p >= 0.05)\n")
}

# Q-Q plot
p_qq <- ggplot(data.frame(sample = residuals_std), aes(sample = sample)) +
  stat_qq(color = "#3498DB", alpha = 0.6) +
  stat_qq_line(color = "#E74C3C", size = 1) +
  labs(title = "Normal Q-Q Plot",
       subtitle = "Checking normality of residuals",
       x = "Theoretical Quantiles", y = "Standardized Residuals") +
  theme_minimal(base_size = 12)

ggsave(file.path(output_dir, "qq_plot.png"), p_qq, width = 8, height = 6, dpi = 150)

# Histogram of residuals with normal overlay
p_hist_resid <- ggplot(data.frame(residual = residuals_std), aes(residual)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50, fill = "#3498DB", alpha = 0.7) +
  stat_function(fun = dnorm, args = list(mean = 0, sd = 1),
                color = "#E74C3C", size = 1.2) +
  labs(title = "Histogram of Standardized Residuals",
       subtitle = "Red curve shows standard normal distribution",
       x = "Standardized Residuals", y = "Density") +
  theme_minimal(base_size = 12)

ggsave(file.path(output_dir, "residual_histogram.png"), p_hist_resid,
       width = 8, height = 6, dpi = 150)

cat("  ✓ Normality diagnostics complete\n\n")

# =============================================================================
# 8. LEVERAGE AND INFLUENCE DIAGNOSTICS
# =============================================================================

cat("Step 8: Checking for high leverage points and influential observations...\n")

# Calculate diagnostic measures
n <- nrow(data_complete)
p <- length(coef(model_final))  # number of parameters including intercept
leverage <- hatvalues(model_final)
cooks_d <- cooks.distance(model_final)
dffits_vals <- dffits(model_final)
dfbetas_vals <- dfbetas(model_final)

# Thresholds
leverage_threshold <- 2 * p / n
cooks_threshold <- 4 / n
dffits_threshold <- 2 * sqrt(p / n)

# Identify problematic points
high_leverage <- which(leverage > leverage_threshold)
high_cooks <- which(cooks_d > cooks_threshold)
high_dffits <- which(abs(dffits_vals) > dffits_threshold)

cat(sprintf("\n  Leverage threshold (2p/n): %.4f\n", leverage_threshold))
cat(sprintf("  Cook's D threshold (4/n): %.4f\n", cooks_threshold))
cat(sprintf("  DFFITS threshold (2√(p/n)): %.4f\n", dffits_threshold))

cat(sprintf("\n  High leverage points: %d (%.2f%%)\n", 
            length(high_leverage), 100 * length(high_leverage) / n))
cat(sprintf("  High Cook's D points: %d (%.2f%%)\n", 
            length(high_cooks), 100 * length(high_cooks) / n))
cat(sprintf("  High DFFITS points: %d (%.2f%%)\n", 
            length(high_dffits), 100 * length(high_dffits) / n))

# Cook's Distance plot
p_cooks <- ggplot(data.frame(index = 1:n, cooks_d = cooks_d), aes(index, cooks_d)) +
  geom_point(alpha = 0.5, color = "#3498DB") +
  geom_hline(yintercept = cooks_threshold, linetype = "dashed", color = "red") +
  labs(title = "Cook's Distance",
       subtitle = sprintf("Red line: threshold = 4/n = %.4f", cooks_threshold),
       x = "Observation Index", y = "Cook's Distance") +
  theme_minimal(base_size = 12)

ggsave(file.path(output_dir, "cooks_distance.png"), p_cooks,
       width = 10, height = 6, dpi = 150)

# Residuals vs Leverage plot
p_resid_lev <- ggplot(data.frame(leverage = leverage, residual = residuals_std,
                                 cooks = cooks_d),
                      aes(leverage, residual)) +
  geom_point(aes(size = cooks, color = cooks > cooks_threshold), alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = leverage_threshold, linetype = "dashed", color = "red") +
  scale_color_manual(values = c("TRUE" = "#E74C3C", "FALSE" = "#3498DB"),
                     labels = c("TRUE" = "High Cook's D", "FALSE" = "Normal")) +
  labs(title = "Residuals vs Leverage",
       subtitle = "Point size shows Cook's distance; Red line: high leverage threshold",
       x = "Leverage", y = "Standardized Residuals",
       color = "Influence", size = "Cook's D") +
  theme_minimal(base_size = 12)

ggsave(file.path(output_dir, "residuals_vs_leverage.png"), p_resid_lev,
       width = 10, height = 6, dpi = 150)

# DFFITS plot
p_dffits <- ggplot(data.frame(index = 1:n, dffits = dffits_vals), aes(index, dffits)) +
  geom_point(alpha = 0.5, color = "#9B59B6") +
  geom_hline(yintercept = c(-dffits_threshold, dffits_threshold), 
             linetype = "dashed", color = "red") +
  labs(title = "DFFITS",
       subtitle = sprintf("Red lines: threshold = ±%.4f", dffits_threshold),
       x = "Observation Index", y = "DFFITS") +
  theme_minimal(base_size = 12)

ggsave(file.path(output_dir, "dffits.png"), p_dffits, width = 10, height = 6, dpi = 150)

# Save influence diagnostics
influence_df <- data.frame(
  Index = 1:n,
  Leverage = leverage,
  Cooks_D = cooks_d,
  DFFITS = dffits_vals,
  High_Leverage = leverage > leverage_threshold,
  High_Cooks = cooks_d > cooks_threshold,
  High_DFFITS = abs(dffits_vals) > dffits_threshold
) %>% arrange(desc(Cooks_D))

write.csv(influence_df, file.path(output_dir, "influence_diagnostics.csv"), row.names = FALSE)

# Top 10 most influential points
cat("\n  Top 10 most influential observations (by Cook's D):\n")
print(head(influence_df, 10), row.names = FALSE)

cat("\n  ✓ Influence diagnostics complete\n\n")

# =============================================================================
# 9. OUTLIER DETECTION
# =============================================================================

cat("Step 9: Identifying outliers...\n")

# Outliers based on standardized residuals (|r| > 3)
outliers_resid <- which(abs(residuals_std) > 3)
cat(sprintf("\n  Outliers (|standardized residual| > 3): %d (%.2f%%)\n", 
            length(outliers_resid), 100 * length(outliers_resid) / n))

if (length(outliers_resid) > 0) {
  outlier_df <- data.frame(
    Index = outliers_resid,
    Std_Residual = residuals_std[outliers_resid],
    Fitted = fitted_vals[outliers_resid],
    Leverage = leverage[outliers_resid],
    Cooks_D = cooks_d[outliers_resid]
  ) %>% arrange(desc(abs(Std_Residual)))
  
  cat("\n  Top outliers by residual magnitude:\n")
  print(head(outlier_df, 10), row.names = FALSE)
  
  write.csv(outlier_df, file.path(output_dir, "outliers.csv"), row.names = FALSE)
}

cat("  ✓ Outlier detection complete\n\n")

# =============================================================================
# 10. COMPREHENSIVE DIAGNOSTIC PLOT PANEL
# =============================================================================

cat("Step 10: Creating comprehensive diagnostic plot panel...\n")

# Standard 4-panel diagnostic plot
png(file.path(output_dir, "diagnostic_panel_4plot.png"), 
    width = 12, height = 10, units = "in", res = 150)
par(mfrow = c(2, 2))
plot(model_final, which = 1:4)
dev.off()

# Create custom 6-panel plot
png(file.path(output_dir, "diagnostic_panel_6plot.png"), 
    width = 15, height = 10, units = "in", res = 150)
par(mfrow = c(2, 3))
plot(model_final, which = 1)  # Residuals vs Fitted
plot(model_final, which = 2)  # Q-Q
plot(model_final, which = 3)  # Scale-Location
plot(model_final, which = 4)  # Cook's distance
plot(model_final, which = 5)  # Residuals vs Leverage
plot(model_final, which = 6)  # Cook's vs Leverage
dev.off()

cat("  ✓ Diagnostic panels created\n\n")

# =============================================================================
# 11. FINAL SUMMARY AND RECOMMENDATIONS
# =============================================================================

cat(strrep("=", 80), "\n")
cat("ANALYSIS SUMMARY AND RECOMMENDATIONS\n")
cat(strrep("=", 80), "\n\n")

cat("MODEL SPECIFICATION:\n")
cat(sprintf("  Response: gdp_per_capita\n"))
cat(sprintf("  Predictors: %d variables\n", length(current_predictors)))
cat(sprintf("  Sample size: %d observations\n", n))
cat(sprintf("  R-squared: %.4f\n", summary(model_final)$r.squared))
cat(sprintf("  Adjusted R-squared: %.4f\n", summary(model_final)$adj.r.squared))
cat(sprintf("  Residual SE: %.4f\n", summary(model_final)$sigma))

cat("\nASSUMPTION CHECKS:\n")

# Linearity
cat("  1. Linearity: Check residuals vs fitted plot\n")
cat("     ✓ Plot saved: residuals_vs_fitted.png\n")

# Constant variance
cat("\n  2. Constant Variance (Homoscedasticity):\n")
cat(sprintf("     Breusch-Pagan test p-value: %.4f\n", bp_test$p.value))
if (bp_test$p.value < 0.05) {
  cat("     ⚠ Evidence of heteroscedasticity - consider transformation or WLS\n")
} else {
  cat("     ✓ Homoscedasticity assumption satisfied\n")
}

# Normality
cat("\n  3. Normality of Residuals:\n")
cat(sprintf("     Shapiro-Wilk test p-value: %.4e\n", sw_test$p.value))
if (sw_test$p.value < 0.05) {
  cat("     ⚠ Residuals may not be perfectly normal\n")
  cat("     Note: With large n, OLS is robust to moderate non-normality (CLT)\n")
} else {
  cat("     ✓ Normality assumption satisfied\n")
}

# Independence (note: can't test fully without time series context)
cat("\n  4. Independence:\n")
cat("     Note: Data spans multiple years and countries\n")
cat("     Consider: Panel data methods, clustering, or time effects\n")

# Multicollinearity
cat("\n  5. Multicollinearity:\n")
max_vif_final <- max(vif_reduced)
cat(sprintf("     Max VIF in final model: %.2f\n", max_vif_final))
if (max_vif_final < 5) {
  cat("     ✓ No multicollinearity concerns (all VIF < 5)\n")
} else if (max_vif_final < 10) {
  cat("     ⚠ Moderate multicollinearity (VIF between 5-10)\n")
} else {
  cat("     ⚠⚠ High multicollinearity (VIF > 10) - consider further reduction\n")
}

# Influential points
cat("\n  6. Influential Observations:\n")
cat(sprintf("     High Cook's D points: %d (%.2f%%)\n", 
            length(high_cooks), 100 * length(high_cooks) / n))
cat(sprintf("     Outliers (|r| > 3): %d (%.2f%%)\n", 
            length(outliers_resid), 100 * length(outliers_resid) / n))
if (length(high_cooks) > 0.05 * n) {
  cat("     ⚠ >5% influential points - investigate further\n")
} else {
  cat("     ✓ Influential points within acceptable range\n")
}

cat("\nOUTPUT FILES:\n")
cat("  All results saved to:", output_dir, "\n")
cat("  - Model summaries (full and reduced)\n")
cat("  - VIF analysis (initial and reduced)\n")
cat("  - Correlation matrices\n")
cat("  - Diagnostic plots (12 files)\n")
cat("  - Influence diagnostics CSV\n")
cat("  - Outliers CSV\n")

cat("\n", strrep("=", 80), "\n")
cat("ANALYSIS COMPLETE!\n")
cat(strrep("=", 80), "\n\n")

# Save session info
sink(file.path(output_dir, "session_info.txt"))
cat("R Session Information\n")
cat(strrep("=", 80), "\n")
print(sessionInfo())
sink()
