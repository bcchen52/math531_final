# -----------------------------------------------------------------------------
# Exploratory Data Analysis for OWID CO2 dataset
# -----------------------------------------------------------------------------
# This script performs a comprehensive EDA on owid-co2-data.csv.
# It:
#   1. Loads data and defines helper functions
#   2. Examines structure, column types, basic summaries
#   3. Profiles missingness (overall, per variable, by year)
#   4. Creates distributions (histograms/density) for key numeric variables
#   5. Generates correlation matrix + heatmap for selected variables
#   6. Explores temporal trends globally and for top emitting countries
#   7. Examines per-capita vs total emissions relationships
#   8. Looks at composition of emissions (coal/oil/gas/cement/flaring)
#   9. Identifies outliers and heavy-tailed distributions
#  10. Saves plots into ./plots for reproducible artifacts
# All plot objects are saved; any interactive device is optional.
# -----------------------------------------------------------------------------

# Load packages ---------------------------------------------------------------
# Using tidyverse for data manipulation and ggplot2, scales for formatting.
# Using patchwork for multi-panel layout where helpful.
# Using reshape2 or tidyr for restructuring.

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
})

# Ensure reproducible plot output directory
plots_dir <- "plots"
if (!dir.exists(plots_dir)) dir.create(plots_dir)

# Data Load ------------------------------------------------------------------
raw <- read.csv("owid-co2-data.csv", stringsAsFactors = FALSE)

# Basic Structure -------------------------------------------------------------
message("Rows: ", nrow(raw), " | Columns: ", ncol(raw))
message("First 5 columns: ")
print(names(raw)[1:5])

# Overview of column types
column_types <- tibble(variable = names(raw), type = map_chr(raw, ~ class(.x)[1]))
print(head(column_types, 20))
write.csv(column_types, file = file.path(plots_dir, "column_types.csv"), row.names = FALSE)

# Missingness Profiling ------------------------------------------------------
# Overall NA fraction per column
miss_tbl <- raw %>% summarise(across(everything(), ~ mean(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "pct_missing") %>%
  arrange(desc(pct_missing))

miss_plot <- ggplot(miss_tbl %>% filter(pct_missing > 0),
                    aes(x = reorder(variable, pct_missing), y = pct_missing)) +
  geom_col(fill = "#2E86AB") +
  scale_y_continuous(labels = percent_format()) +
  labs(title = "Fraction Missing by Variable", x = "Variable", y = "% Missing") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

if (nrow(miss_tbl %>% filter(pct_missing > 0)) > 0) {
  ggsave(file.path(plots_dir, "missing_fraction.png"), miss_plot, width = 12, height = 6, dpi = 150)
}

# Missingness by year for key metrics
key_metrics <- c("co2", "co2_per_capita", "gdp", "population")
miss_year <- raw %>%
  select(year, all_of(key_metrics)) %>%
  pivot_longer(-year, names_to = "metric", values_to = "value") %>%
  group_by(year, metric) %>%
  summarise(n = n(), n_na = sum(is.na(value)), pct_na = n_na / n, .groups = "drop")

miss_year_plot <- ggplot(miss_year, aes(year, pct_na, color = metric)) +
  geom_line(alpha = 0.8) +
  scale_y_continuous(labels = percent_format()) +
  labs(title = "Missingness Over Time (Key Metrics)", x = "Year", y = "% NA") +
  theme_minimal(base_size = 12)

ggsave(file.path(plots_dir, "missing_over_time.png"), miss_year_plot, width = 10, height = 5, dpi = 150)

# Distribution Plots ---------------------------------------------------------
# Focus on selected numeric variables with meaningful non-zero values.
num_vars <- c("co2", "co2_per_capita", "gdp", "population", "coal_co2", "oil_co2", "gas_co2")

long_num <- raw %>% select(all_of(num_vars)) %>% pivot_longer(everything(), names_to = "variable", values_to = "value")

# Remove zeros-only to avoid misleading densities; keep positive subset for log density.
long_num_pos <- long_num %>% filter(!is.na(value), value > 0)

hist_plot <- ggplot(long_num_pos, aes(value)) +
  geom_histogram(bins = 50, fill = "#0072B2", color = "white") +
  facet_wrap(~ variable, scales = "free") +
  scale_x_continuous(labels = comma) +
  labs(title = "Distributions (Positive Values)", x = "Value", y = "Count") +
  theme_minimal(base_size = 11)

ggsave(file.path(plots_dir, "histograms_positive.png"), hist_plot, width = 14, height = 8, dpi = 150)

log_hist_plot <- ggplot(long_num_pos, aes(log10(value))) +
  geom_histogram(bins = 50, fill = "#D55E00", color = "white") +
  facet_wrap(~ variable, scales = "free") +
  labs(title = "Log10 Distributions (Positive Values)", x = "log10(Value)", y = "Count") +
  theme_minimal(base_size = 11)

ggsave(file.path(plots_dir, "histograms_log10.png"), log_hist_plot, width = 14, height = 8, dpi = 150)

# Correlation Matrix ---------------------------------------------------------
# Use pairwise complete observations among selected emission variables.
cor_vars <- c("co2", "coal_co2", "oil_co2", "gas_co2", "cement_co2", "flaring_co2")
cor_df <- raw %>% select(all_of(cor_vars))
cor_mat <- cor(cor_df, use = "pairwise.complete.obs")

# Convert to long for heatmap
cor_long <- as.data.frame(cor_mat) %>% mutate(var1 = rownames(cor_mat)) %>% pivot_longer(-var1, names_to = "var2", values_to = "corr")

heat_plot <- ggplot(cor_long, aes(var1, var2, fill = corr)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#2166AC", mid = "#FFFFFF", high = "#B2182B", midpoint = 0, limits = c(-1,1)) +
  labs(title = "Correlation Heatmap: Emission Components", x = NULL, y = NULL, fill = "r") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(plots_dir, "correlation_heatmap.png"), heat_plot, width = 7, height = 6, dpi = 150)

# Temporal Trends ------------------------------------------------------------
# Aggregate global totals per year (summing countries with available data).
annual_global <- raw %>% group_by(year) %>% summarise(total_co2 = sum(co2, na.rm = TRUE), total_co2_luc = sum(co2_including_luc, na.rm = TRUE))

global_trend_plot <- ggplot(annual_global, aes(year)) +
  geom_line(aes(y = total_co2), color = "#1B9E77", size = 0.9) +
  geom_line(aes(y = total_co2_luc), color = "#D95F02", size = 0.9, linetype = "dashed") +
  labs(title = "Global CO2 Emissions Over Time", y = "Total CO2 (sum across countries)", x = "Year",
       subtitle = "Solid: fossil + industry; Dashed: including land-use change") +
  theme_minimal(base_size = 12)

ggsave(file.path(plots_dir, "global_co2_trend.png"), global_trend_plot, width = 10, height = 5, dpi = 150)

# Top Countries (Recent Year) ------------------------------------------------
latest_year <- max(raw$year, na.rm = TRUE)
latest_co2 <- raw %>% filter(year == latest_year) %>% select(country, co2) %>% filter(!is.na(co2)) %>% arrange(desc(co2)) %>% slice_head(n = 15)

bar_top_countries <- ggplot(latest_co2, aes(x = reorder(country, co2), y = co2)) +
  geom_col(fill = "#6A3D9A") +
  coord_flip() +
  scale_y_continuous(labels = comma) +
  labs(title = paste("Top 15 CO2 Emitters in", latest_year), x = NULL, y = "CO2") +
  theme_minimal(base_size = 12)

ggsave(file.path(plots_dir, "top_emitters_latest_year.png"), bar_top_countries, width = 8, height = 6, dpi = 150)

# Per Capita vs Total --------------------------------------------------------
percap_vs_total <- raw %>% filter(year == latest_year, !is.na(co2_per_capita), !is.na(co2), !is.na(population))

scatter_percap_total <- ggplot(percap_vs_total, aes(co2, co2_per_capita)) +
  geom_point(alpha = 0.6, color = "#2C3E50") +
  scale_x_log10(labels = comma) +
  scale_y_log10() +
  labs(title = paste("Total vs Per-Capita CO2 (", latest_year, ")", sep = ""), x = "Total CO2 (log10)", y = "Per-Capita CO2 (log10)") +
  theme_minimal(base_size = 12)

ggsave(file.path(plots_dir, "total_vs_percapita.png"), scatter_percap_total, width = 7, height = 5, dpi = 150)

# Emission Composition -------------------------------------------------------
composition_vars <- c("coal_co2", "oil_co2", "gas_co2", "cement_co2", "flaring_co2")
composition_latest <- raw %>% filter(year == latest_year) %>% select(country, all_of(composition_vars)) %>%
  pivot_longer(-country, names_to = "component", values_to = "value") %>% group_by(component) %>% summarise(total = sum(value, na.rm = TRUE)) %>% mutate(share = total / sum(total))

comp_plot <- ggplot(composition_latest, aes(x = reorder(component, total), y = total, fill = component)) +
  geom_col(show.legend = FALSE) +
  scale_y_continuous(labels = comma) +
  labs(title = paste("Global Emission Composition (", latest_year, ")", sep = ""), x = "Component", y = "Total CO2") +
  theme_minimal(base_size = 12)

ggsave(file.path(plots_dir, "global_emission_composition.png"), comp_plot, width = 8, height = 5, dpi = 150)

# Outlier Detection (Simple) -------------------------------------------------
# Identify high leverage / extreme per-capita emitters in latest year.
quantiles_pc <- quantile(percap_vs_total$co2_per_capita, probs = c(0.25, 0.5, 0.75, 0.95), na.rm = TRUE)

extreme_percap <- percap_vs_total %>% filter(co2_per_capita > quantiles_pc["95%"]) %>% arrange(desc(co2_per_capita)) %>% select(country, co2_per_capita, co2)
write.csv(extreme_percap, file.path(plots_dir, "extreme_per_capita_emitters.csv"), row.names = FALSE)

# Save summary tables --------------------------------------------------------
write.csv(miss_tbl, file.path(plots_dir, "missing_fraction_table.csv"), row.names = FALSE)
write.csv(latest_co2, file.path(plots_dir, "top_emitters_table.csv"), row.names = FALSE)
write.csv(composition_latest, file.path(plots_dir, "composition_latest.csv"), row.names = FALSE)

# Notes & Inline Commentary --------------------------------------------------
# - Missingness is heavy in early years (many zeros / blanks). Filtering years with non-trivial emissions is
#   advisable before modeling (e.g., year >= 1850 or when cumulative_co2 > 0 for most countries).
# - Distributions are heavily right-skewed; log-transform helps stabilize variance.
# - Strong correlations expected among fossil components (coal/oil/gas) and total co2, raising multicollinearity concerns.
# - Per-capita vs total shows a wide spread; small population countries can have high per-capita despite lower totals.
# - Composition shows coal and oil likely dominant; shares inform focused decarbonization pathways.
# - Outlier detection simplistic; consider robust methods (e.g., median absolute deviation) for deeper analysis.

message("EDA complete. Plots and tables saved in /plots")
