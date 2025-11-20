#options(repos = c(CRAN = "https://cloud.r-project.org"))
#install.packages(c("tidyverse", "caret", "car", "dplyr"))

library(tidyverse)
library(caret)
library(car)
library(dplyr)

raw_data = read.csv("owid-co2-data.csv")

str(raw_data)
#summary(raw_data)
tail(raw_data)

#many of the early years have no values, remove them
year_na_summary <- raw_data %>%
  group_by(year) %>%
  summarise(
    n_rows   = n(),                                   # number of country-year rows
    n_na     = sum(is.na(co2_per_capita)),            # how many are NA
    n_non_na = n_rows - n_na,                         # how many are not NA
    frac_na  = n_na / n_rows                          # fraction missing
  ) %>%
  arrange(year)

#print(year_na_summary, n = Inf)

# ---------------------------------------------------------------------------
# Missing value counts per variable + bar chart
# ---------------------------------------------------------------------------
# We build a tidy summary of NA counts for every column and plot a bar chart
# showing both raw NA counts and percentage missing via a fill gradient.
# Adjust `vars_to_check` if you want to restrict to a subset.

vars_to_check <- names(raw_data) # you can narrow this vector if desired

na_counts <- raw_data %>%
  summarise(across(all_of(vars_to_check), ~ sum(is.na(.)))) %>%
  tidyr::pivot_longer(everything(), names_to = "variable", values_to = "n_na") %>%
  mutate(
    n_total  = nrow(raw_data),
    pct_na   = n_na / n_total
  ) %>%
  arrange(desc(n_na))

# Optional: filter out variables with zero missing values for a cleaner plot
na_counts_nonzero <- na_counts %>% filter(n_na > 0)

if (nrow(na_counts_nonzero) == 0) {
  message("All variables have 0 missing values; skipping NA plot.")
} else {
  na_counts_plot <- ggplot(na_counts_nonzero, aes(x = reorder(variable, n_na), y = n_na, fill = pct_na)) +
    geom_col() +
    scale_fill_viridis_c(labels = scales::percent_format(accuracy = 1), name = "% NA") +
    labs(
      title = "Missing Values per Variable",
      subtitle = paste("Total rows:", nrow(raw_data)),
      x = "Variable",
      y = "Count of NA"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )
  print(na_counts_plot)
  # Save the plot to a file for non-interactive runs
  ggsave(filename = "missing_values_plot.png", plot = na_counts_plot, width = 10, height = 6, dpi = 150)
  message("Saved NA plot to missing_values_plot.png")
}

# If you also want a wide table of counts (year-wise), the earlier code shows a template.
# View the table directly:
na_counts

# # 1851 is where there are 252 countries total
# data0 = raw_data %>% filter(year >= 1851)

# data1 = data0 %>% drop_na(co2_per_capita)

# #this is the cleaned data
# data1 = data1 %>% mutate(country = factor(country), gdp_per_capita = gdp/population, trade_co2_per_capita = trade_co2/population)

# #country
# #population
# #gdp
# #gdp per capita (create)
# #cement_co2_per_capita
# #coal_co2_per_capita
# #oil_co2_per_capita
# #gas_co2_per_capita
# #flaring_co2_per_capita
# #other_co2_per_capita
# #land_use_change_co2_per_capita
# #energy_per_capita
# #trade_co2_per_capita
# #ghg_per_capita

# wanted_columns = c(
#     "country", 
#     "co2_per_capita", 
#     "gdp_per_capita",
#     "cement_co2_per_capita",
#     "coal_co2_per_capita",
#     "oil_co2_per_capita",
#     "gas_co2_per_capita",
#     "flaring_co2_per_capita",
#     "other_co2_per_capita",
#     "land_use_change_co2_per_capita",
#     "energy_per_capita",
#     "trade_co2_per_capita",
#     "ghg_per_capita"
#     )

# kaya_identity = c(
#     "co2",
#     "population",
#     "gdp_per_capita",
#     "energy_per_gdp",
#     "co2_per_unit_energy"
# )

# data2 = data1 %>% select(all_of(wanted_columns)) %>% drop_na()

# data3 = data1 %>% select(all_of(kaya_identity)) %>% drop_na()

# model0 = lm(co2_per_capita ~ ., data=data2)

# model1 = lm(co2 ~ ., data=data3)

# summary(data2)

# summary(model1)

# data2_numeric = data2 %>% select(where(is.numeric))

# cor_mat <- cor(data2_numeric)
# cor_mat