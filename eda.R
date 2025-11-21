# -----------------------------------------------------------------------------
# Exploratory Data Analysis for OWID CO2 Dataset
# -----------------------------------------------------------------------------
# This script performs comprehensive EDA on owid-co2-data.csv, focusing on
# the 11 variables selected for modeling (per-capita emission metrics + year).
#
# Analysis Pipeline:
#   1. Load data and apply filtering (year >= 1851, valid iso_code, non-NA gdp_per_capita)
#   2. Create derived variable: gdp_per_capita
#   3. Profile missingness for selected variables
#   4. Generate distributions (histograms, log-scale) for all numeric predictors
#   5. Create correlation heatmap for the 10 numeric variables
#   6. Explore temporal trends for key per-capita metrics
#   7. Examine emission composition across fuel types (per-capita basis)
#   8. Identify outliers and data quality issues
#   9. Save all plots and summary tables to ./plots/
#
# Selected Variables (11 total):
#   - iso_code (country identifier)
#   - year (temporal predictor)
#   - gdp_per_capita (RESPONSE VARIABLE - created from gdp/population)
#   - cement_co2_per_capita, coal_co2_per_capita, oil_co2_per_capita,
#     gas_co2_per_capita, flaring_co2_per_capita
#   - land_use_change_co2_per_capita
#   - energy_per_capita
#   - methane_per_capita, nitrous_oxide_per_capita
#
# Removed variables:
#   - co2_per_capita (now excluded, was response in previous version)
#   - trade_co2_per_capita (82% missing)
#   - other_co2_per_capita (92% missing)
# -----------------------------------------------------------------------------

# Load packages ---------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
})

# Ensure plots directory exists
plots_dir <- "plots"
if (!dir.exists(plots_dir)) dir.create(plots_dir)

# Data Load and Cleaning -----------------------------------------------------
# Load raw data
raw <- read.csv("owid-co2-data.csv", stringsAsFactors = FALSE)

# Apply filtering logic:
# 1. Year >= 1851 (where 252 countries have data)
# 2. Valid iso_code (exclude empty strings and NA)
# 3. Non-NA gdp and population (to create gdp_per_capita, our response variable)
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

# Define the 11 wanted columns (gdp_per_capita is response, removed co2_per_capita)
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

# Select only the wanted columns
data_selected <- data1 %>% select(all_of(wanted_columns))

message("Data filtered and selected:")
message("  Rows: ", nrow(data_selected))
message("  Columns: ", ncol(data_selected))
message("  Unique countries: ", n_distinct(data_selected$iso_code))
message("  Year range: ", min(data_selected$year), " to ", max(data_selected$year))

# Missingness Profiling ------------------------------------------------------
# Calculate missingness for each of the 10 numeric variables (excluding iso_code and year)
numeric_vars <- wanted_columns[!wanted_columns %in% c("iso_code", "year")]

miss_tbl <- data_selected %>%
  select(all_of(numeric_vars)) %>%
  summarise(across(everything(), ~ mean(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "pct_missing") %>%
  arrange(desc(pct_missing))

message("\nMissingness in selected variables:")
print(miss_tbl, n = Inf)

# Visualize missingness
miss_plot <- ggplot(miss_tbl, aes(x = reorder(variable, pct_missing), y = pct_missing)) +
  geom_col(fill = "#2E86AB") +
  scale_y_continuous(labels = percent_format()) +
  labs(
    title = "Missingness in Selected Per-Capita Variables",
    subtitle = "10 numeric variables (excluding year) - After filtering for year >= 1851 and non-NA gdp_per_capita",
    x = "Variable",
    y = "% Missing"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

ggsave(file.path(plots_dir, "selected_vars_missingness.png"), miss_plot,
       width = 10, height = 6, dpi = 150)

# Save missingness table
write.csv(miss_tbl, file.path(plots_dir, "selected_vars_missing_table.csv"),
          row.names = FALSE)

# Save missingness table
write.csv(miss_tbl, file.path(plots_dir, "selected_vars_missing_table.csv"),
          row.names = FALSE)

# Distribution Plots ---------------------------------------------------------
# Create histograms for all 13 numeric variables on original scale
long_num <- data_selected %>%
  select(all_of(numeric_vars)) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "value")

# Original scale (positive values only to avoid log(0) issues)
long_num_pos <- long_num %>% filter(!is.na(value), value > 0)

hist_plot <- ggplot(long_num_pos, aes(value)) +
  geom_histogram(bins = 50, fill = "#0072B2", color = "white") +
  facet_wrap(~ variable, scales = "free", ncol = 3) +
  scale_x_continuous(labels = comma) +
  labs(
    title = "Distributions of Selected Per-Capita Variables",
    subtitle = "Positive values only (original scale)",
    x = "Value",
    y = "Count"
  ) +
  theme_minimal(base_size = 10) +
  theme(strip.text = element_text(size = 8))

ggsave(file.path(plots_dir, "selected_vars_histograms.png"), hist_plot,
       width = 14, height = 12, dpi = 150)

# Log10 scale for better visibility of distributions
log_hist_plot <- ggplot(long_num_pos, aes(log10(value))) +
  geom_histogram(bins = 50, fill = "#D55E00", color = "white") +
  facet_wrap(~ variable, scales = "free", ncol = 3) +
  labs(
    title = "Log10 Distributions of Selected Per-Capita Variables",
    subtitle = "Positive values only - log scale reveals skewness patterns",
    x = "log10(Value)",
    y = "Count"
  ) +
  theme_minimal(base_size = 10) +
  theme(strip.text = element_text(size = 8))

ggsave(file.path(plots_dir, "selected_vars_histograms_log10.png"), log_hist_plot,
       width = 14, height = 12, dpi = 150)

# Correlation Matrix and Heatmap ---------------------------------------------
# Compute correlation for the 13 numeric variables using pairwise complete obs
data_numeric <- data_selected %>%
  select(all_of(numeric_vars)) %>%
  drop_na()

message("\nComplete cases for correlation analysis: ", nrow(data_numeric))

cor_mat <- cor(data_numeric, use = "pairwise.complete.obs")

# Convert to long format for ggplot heatmap
cor_long <- as.data.frame(cor_mat) %>%
  mutate(var1 = rownames(cor_mat)) %>%
  pivot_longer(-var1, names_to = "var2", values_to = "corr")

# Create correlation heatmap
heat_plot <- ggplot(cor_long, aes(var1, var2, fill = corr)) +
  geom_tile(color = "white", size = 0.3) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
    midpoint = 0, limits = c(-1, 1),
    name = "Correlation"
  ) +
  labs(
    title = "Correlation Heatmap: Selected Per-Capita Variables",
    subtitle = "10 numeric variables for GDP per capita modeling (year analyzed separately)",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9),
    axis.text.y = element_text(size = 9),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave(file.path(plots_dir, "selected_vars_correlation_heatmap.png"), heat_plot,
       width = 12, height = 10, dpi = 150)

# Save correlation matrix as CSV
write.csv(cor_mat, file.path(plots_dir, "selected_vars_correlation_matrix.csv"),
          row.names = TRUE)

message("Correlation heatmap created with ", nrow(cor_mat), " variables")

message("Correlation heatmap created with ", nrow(cor_mat), " variables")

# Temporal Trends ------------------------------------------------------------
# Track how key per-capita metrics evolve over time (global average)
# Year is already included in data_selected
data_with_year <- data_selected

# Calculate global averages per year for selected metrics
key_percap_metrics <- c("gdp_per_capita", "energy_per_capita", "coal_co2_per_capita",
                        "methane_per_capita", "nitrous_oxide_per_capita")

temporal_trends <- data_with_year %>%
  select(year, all_of(key_percap_metrics)) %>%
  pivot_longer(-year, names_to = "metric", values_to = "value") %>%
  group_by(year, metric) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    median_value = median(value, na.rm = TRUE),
    .groups = "drop"
  )

trend_plot <- ggplot(temporal_trends, aes(year, mean_value, color = metric)) +
  geom_line(size = 0.8, alpha = 0.8) +
  facet_wrap(~ metric, scales = "free_y", ncol = 2) +
  labs(
    title = "Temporal Trends: Global Mean Per-Capita Metrics",
    subtitle = "Average across all countries per year",
    x = "Year",
    y = "Mean Value (per capita)"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none", strip.text = element_text(face = "bold"))

ggsave(file.path(plots_dir, "selected_vars_temporal_trends.png"), trend_plot,
       width = 12, height = 8, dpi = 150)

ggsave(file.path(plots_dir, "selected_vars_temporal_trends.png"), trend_plot,
       width = 12, height = 8, dpi = 150)

# Emission Composition Analysis (Per-Capita) --------------------------------
# Examine the relative contribution of different fuel types to per-capita emissions
fuel_types_percap <- c("cement_co2_per_capita", "coal_co2_per_capita",
                       "oil_co2_per_capita", "gas_co2_per_capita",
                       "flaring_co2_per_capita")

# Get latest year data
latest_year <- max(data_with_year$year, na.rm = TRUE)
composition_latest <- data_with_year %>%
  filter(year == latest_year) %>%
  select(iso_code, all_of(fuel_types_percap)) %>%
  pivot_longer(-iso_code, names_to = "fuel_type", values_to = "value") %>%
  group_by(fuel_type) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    median_value = median(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(fuel_type = str_remove(fuel_type, "_per_capita"))

comp_plot <- ggplot(composition_latest, aes(x = reorder(fuel_type, mean_value),
                                             y = mean_value, fill = fuel_type)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  scale_y_continuous(labels = comma) +
  labs(
    title = paste("Mean Per-Capita Emission Composition (", latest_year, ")", sep = ""),
    subtitle = "Average across all countries",
    x = "Fuel Type",
    y = "Mean Per-Capita CO2"
  ) +
  theme_minimal(base_size = 12)

ggsave(file.path(plots_dir, "selected_vars_fuel_composition.png"), comp_plot,
       width = 8, height = 6, dpi = 150)

# Save composition table
write.csv(composition_latest, file.path(plots_dir, "selected_vars_fuel_composition_table.csv"),
          row.names = FALSE)

write.csv(composition_latest, file.path(plots_dir, "selected_vars_fuel_composition_table.csv"),
          row.names = FALSE)

# Top and Bottom Countries by GDP Per Capita --------------------------------
# Identify highest GDP per capita countries in the latest year
latest_data <- data_with_year %>% filter(year == latest_year)

# Top 15 by GDP per capita
top_percap <- latest_data %>%
  select(iso_code, gdp_per_capita) %>%
  filter(!is.na(gdp_per_capita)) %>%
  arrange(desc(gdp_per_capita)) %>%
  slice_head(n = 15)

bar_top_percap <- ggplot(top_percap, aes(x = reorder(iso_code, gdp_per_capita),
                                          y = gdp_per_capita)) +
  geom_col(fill = "#6A3D9A") +
  coord_flip() +
  scale_y_continuous(labels = comma) +
  labs(
    title = paste("Top 15 Countries by GDP Per Capita in", latest_year),
    x = NULL,
    y = "GDP Per Capita"
  ) +
  theme_minimal(base_size = 12)

ggsave(file.path(plots_dir, "selected_vars_top_percap_gdp.png"), bar_top_percap,
       width = 8, height = 6, dpi = 150)

# Save top GDP table
write.csv(top_percap, file.path(plots_dir, "selected_vars_top_percap_table.csv"),
          row.names = FALSE)

# Pairwise Scatterplots: Key Relationships -----------------------------------
# Examine relationships between response (GDP) and key predictors
latest_complete <- latest_data %>%
  select(gdp_per_capita, energy_per_capita,
         coal_co2_per_capita, oil_co2_per_capita) %>%
  drop_na()

# Energy vs GDP per capita
scatter_energy <- ggplot(latest_complete, aes(energy_per_capita, gdp_per_capita)) +
  geom_point(alpha = 0.5, color = "#2C3E50") +
  geom_smooth(method = "lm", color = "#E74C3C", se = TRUE, alpha = 0.2) +
  scale_x_log10(labels = comma) +
  scale_y_log10(labels = comma) +
  labs(
    title = "Energy Per Capita vs GDP Per Capita",
    subtitle = paste("Latest year:", latest_year),
    x = "Energy Per Capita (log10)",
    y = "GDP Per Capita (log10)"
  ) +
  theme_minimal(base_size = 11)

ggsave(file.path(plots_dir, "selected_vars_energy_vs_gdp.png"), scatter_energy,
       width = 8, height = 6, dpi = 150)

# Coal CO2 vs GDP per capita
scatter_coal <- ggplot(latest_complete, aes(coal_co2_per_capita, gdp_per_capita)) +
  geom_point(alpha = 0.5, color = "#2C3E50") +
  geom_smooth(method = "lm", color = "#3498DB", se = TRUE, alpha = 0.2) +
  scale_x_log10(labels = comma) +
  scale_y_log10(labels = comma) +
  labs(
    title = "Coal CO2 Per Capita vs GDP Per Capita",
    subtitle = paste("Latest year:", latest_year),
    x = "Coal CO2 Per Capita (log10)",
    y = "GDP Per Capita (log10)"
  ) +
  theme_minimal(base_size = 11)

ggsave(file.path(plots_dir, "selected_vars_coal_vs_gdp.png"), scatter_coal,
       width = 8, height = 6, dpi = 150)

# Data Quality Summary -------------------------------------------------------
# Generate summary statistics for all numeric variables
summary_stats <- data_selected %>%
  select(all_of(numeric_vars)) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
  group_by(variable) %>%
  summarise(
    n = sum(!is.na(value)),
    mean = mean(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    min = min(value, na.rm = TRUE),
    max = max(value, na.rm = TRUE),
    q25 = quantile(value, 0.25, na.rm = TRUE),
    q75 = quantile(value, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(variable)

write.csv(summary_stats, file.path(plots_dir, "selected_vars_summary_statistics.csv"),
          row.names = FALSE)

message("\nSummary Statistics:")
print(summary_stats, n = Inf)

# Final Summary and Notes ----------------------------------------------------
message("\n", strrep("=", 80))
message("EDA COMPLETE - SELECTED PER-CAPITA VARIABLES")
message(strrep("=", 80))
message("\nKEY FINDINGS:")
message("  - Dataset filtered: year >= 1851, non-NA gdp_per_capita")
message("  - Response variable: gdp_per_capita")
message("  - Variables analyzed: 11 total (iso_code + year + 9 numeric predictors)")
message("  - Removed: co2_per_capita, trade_co2_per_capita, other_co2_per_capita")
message("  - Added: year (temporal predictor)")
message("  - Unique countries: ", n_distinct(data_selected$iso_code))
message("  - Total observations: ", nrow(data_selected))
message("  - Year range: ", min(data_selected$year), " to ", max(data_selected$year))
message("\nMODELING NOTES:")
message("  - GDP per capita is response; emission metrics are predictors")
message("  - All per-capita metrics are right-skewed (log transform recommended)")
message("  - Year can capture temporal trends in economic development")
message("  - Strong correlations among fuel types (multicollinearity risk)")
message("  - Some variables have substantial missingness (handle carefully)")
message("  - Energy and emissions show relationships with GDP development")
message("\nOUTPUTS SAVED to ./", plots_dir, "/")
message("  Plots: 9 files")
message("  Tables: 5 CSV files")
message(strrep("=", 80), "\n")
