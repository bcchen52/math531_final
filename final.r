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

#unique_countries = unique(raw_data$country)
unique_countries = unique(raw_data$iso_code)

# 1851 is where there are 252 countries total
data0 = raw_data %>% filter(year >= 1851) %>% filter(iso_code != "", !is.na(iso_code))
# there is an empty string for iso_code as well as NA

unique_countries = unique(data0$country)


data1 = data0 %>% drop_na(co2_per_capita)

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
    "co2_per_capita", 
    "gdp_per_capita",
    "gdp",
    "population",
    "cement_co2_per_capita",
    "coal_co2_per_capita",
    "oil_co2_per_capita",
    "gas_co2_per_capita",
    "flaring_co2_per_capita",
    "other_co2_per_capita",
    "land_use_change_co2_per_capita",
    "energy_per_capita",
    "trade_co2_per_capita",
    "ghg_per_capita"
    )

kaya_identity = c(
    "co2",
    "population",
    "gdp_per_capita",
    "energy_per_gdp",
    "co2_per_unit_energy"
)

data2 = data1 %>% select(all_of(wanted_columns))

data2_numeric = data2 %>% select(where(is.numeric)) %>% drop_na()

str(data2_numeric)

head(data2_numeric)

tail(data2_numeric)

cor_mat <- cor(data2_numeric)
cor_mat

#data3 = data1 %>% select(all_of(kaya_identity)) %>% drop_na()

data_c = unique(data2$country)

data_c

model0 = lm(co2_per_capita ~ ., data=data2)

#model1 = lm(co2 ~ ., data=data3)

#summary(data2)

#summary(model0)