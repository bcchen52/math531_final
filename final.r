#options(repos = c(CRAN = "https://cloud.r-project.org"))
#install.packages(c("tidyverse", "caret", "car", "dplyr", "ggplot2"))

library(tidyverse)
library(caret)
library(car)
library(dplyr)
library(ggplot2)

raw_data = read.csv("owid-co2-data.csv")

#str(raw_data)
#summary(raw_data)
#tail(raw_data)

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

#unique_countries = unique(raw_data$country)
unique_countries = unique(raw_data$iso_code)

# 1851 is where there are 252 countries total
data0 = raw_data %>% filter(year >= 1851) %>% filter(iso_code != "", !is.na(iso_code))
# there is an empty string for iso_code as well as NA

unique_countries = unique(data0$country)

data1 = data0 %>% drop_na(coal_co2_per_capita)

#this is the cleaned data
data1 = data1 %>% mutate(iso_code = factor(iso_code), gdp_per_capita = gdp/population, trade_co2_per_capita = trade_co2/population)

#country
#population
#gdp
#gdp per capita (create)
#cement_co2_per_capita
#coal_co2_per_capita
#oil_co2_per_capita
#gas_co2_per_capita
#flaring_co2_per_capita
#other_co2_per_capita
#land_use_change_co2_per_capita
#energy_per_capita
#trade_co2_per_capita
#ghg_per_capita

wanted_columns = c(
    "iso_code",
    "co2_including_luc_per_capita",
    "co2_per_capita",
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

cause_cols <- c(
  "cement_co2_per_capita",
  "coal_co2_per_capita",
  "oil_co2_per_capita",
  "gas_co2_per_capita",
  "flaring_co2_per_capita",
  "land_use_change_co2_per_capita"
)

kaya_identity = c(
    "co2",
    "population",
    "gdp_per_capita",
    "energy_per_gdp",
    "co2_per_unit_energy"
)

data2 = data1 %>% select(all_of(wanted_columns))

#=====Top 20 CO2 Emissions=====
top_co2_2024 <- data2 %>% 
  filter(
    year == 2024,
    !is.na(co2_per_capita),
    !is.na(iso_code),
    iso_code != ""
  ) %>% 
  arrange(desc(co2_per_capita)) %>% 
  slice_head(n = 10)

# Bar plot: top CO2 per capita in 2024 by iso_code
ggplot(top_co2_2024, aes(x = reorder(iso_code, co2_per_capita),
                         y = co2_per_capita,
                         fill = co2_per_capita)) +
  geom_col() +
  coord_flip() +  # horizontal bars for readability
  scale_fill_viridis_c(option = "plasma") +
  labs(
    title = "Top 10 CO₂ Emissions per Capita in 2024",
    subtitle = "By ISO country code",
    x = "ISO code",
    y = "CO₂ per capita (tonnes per person)",
    fill = "t/person"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 10)
  )

#=====Total Ratios=====
# 2. Sum over ALL years and ALL countries
totals <- data2 %>%
  filter(!is.na(co2_including_luc_per_capita)) %>%
  summarise(
    total_co2_inc_luc_pc = sum(co2_including_luc_per_capita, na.rm = TRUE),
    across(all_of(cause_cols), ~ sum(.x, na.rm = TRUE))
  )

# 3. Compute ratios: total_X / total_CO2_including_LUC
ratios <- totals %>%
  pivot_longer(
    cols = all_of(cause_cols),
    names_to = "cause",
    values_to = "total_cause"
  ) %>%
  mutate(
    ratio     = total_cause / total_co2_inc_luc_pc,
    ratio_pct = 100 * ratio,
    cause = case_when(
      cause == "cement_co2_per_capita"          ~ "Cement",
      cause == "coal_co2_per_capita"            ~ "Coal",
      cause == "oil_co2_per_capita"             ~ "Oil",
      cause == "gas_co2_per_capita"             ~ "Gas",
      cause == "flaring_co2_per_capita"         ~ "Flaring",
      cause == "land_use_change_co2_per_capita" ~ "Land-use change",
      TRUE                                      ~ cause
    )
  )

# 4. Bar plot of overall ratios
ggplot(ratios, aes(x = cause, y = ratio_pct)) +
  geom_col() +
  labs(
    title = "Aggregated contributions relative to CO2 including land-use change",
    subtitle = "Totals summed over all years and countries",
    x = "Source",
    y = "Total X / total CO2 (including LUC) (%)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

data2_numeric = data2 %>% select(where(is.numeric)) %>% drop_na()

#str(data2_numeric)
#head(data2_numeric)
#tail(data2_numeric)

cor_mat <- cor(data2_numeric)
cor_mat

#data3 = data1 %>% select(all_of(kaya_identity)) %>% drop_na()

data_c = unique(data2$country)

#data_c

model0 = lm(co2_per_capita ~ ., data=data2)

#model1 = lm(co2 ~ ., data=data3)

#summary(data2)

#summary(model0)