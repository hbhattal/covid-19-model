## Objective: record case-hospitalization fraction in each age group into csv file to be read in inf_asc.R

## >>> Import data and relevant global variables; set parameters ---------------

# import 'raw_data' and 'date_pulled' into workspace 

# import libraries and clear workspace (but for variables 'raw_data' and 'date_pulled')
library(date)
library(lubridate)
rm(list = ls()[!ls() %in% c("raw_data", "date_pulled")]) 

# set parameters of interest
region <- 'BC'
offset <- 21 # how many days of recent data to ignore

# population data for each of eleven age cohorts (0-9, 10-19, ..., 90+) from 2016 Census
pop_age <- c(457525, 492840, 590560, 607340, 617410, 709300, 611615, 347010, 172765, 40360 + 1325) 

## >>> Pre-process raw data  ---------------------------------------------------

# ignore recent data due to data discrepancies
raw_data <- subset(raw_data, raw_data$reported_date <= (date_pulled-offset))

# subset to HA of interest
if (region != 'BC') {
  raw_data <- subset(raw_data, raw_data$ha_res == region)
}

## >>> Summarize infections and hospitalizations for each age cohort -----------

# function which calculates infections and hospitalizations in age cohort of interest (over date period of interest)
inf_hosp <- function (raw_data, start_age, end_age, start_date, end_date) {
  # subset infection data by date period of interest
  start_date <- as.Date(start_date)
  end_date <- as.Date(end_date)
  raw_data <- subset(raw_data, raw_data$reported_date >= start_date & raw_data$reported_date <= end_date)
  
  # subset infection data by age cohort of interest
  raw_data <- subset(raw_data, raw_data$age_combined >= start_age & raw_data$age_combined <= end_age)
  
  num_I <- nrow(raw_data)
  
  # subset hospitalized cases
  hospitalized_data <- subset(raw_data, raw_data$ever_hospitalized == 'Yes')
  
  # delete entries with incomplete data 
  hospitalized_data <- subset(hospitalized_data, hospitalized_data$admission_date != '')

  num_H <- nrow(hospitalized_data)
  
  # hospitalization rate
  H_e <- num_H/num_I
  
  return(c(num_I, num_H, H_e))
}

# create vectors describing age cohorts (by lower and upper limits)
num_age_groups <- 3
age_low <- c(0, 50, 70)
age_high <- c(49, 69, 150)

# create vectors describing start and end dates for each month from Jan 2020 to Feb 2021
num_months <- 14
start_date_loop <- vector(mode="character", length=num_months)
end_date_loop <- vector(mode="character", length=num_months)

for (i in 1:12) {
  start_date_loop[i] <- paste('2020-', i, '-01', sep='')
  end_date_loop[i] <- paste('2020-', i, '-', as.character(days_in_month(i)), sep='')
}

for (i in 13:14) {
  start_date_loop[i] <- paste('2021-', (i-12), '-01', sep='')
  end_date_loop[i] <- paste('2021-', (i-12), '-', as.character(days_in_month(i-12)), sep='')
}


# functions calculating monthly hospitalization fraction  for age cohort of interest for each month
monthly_hosp_rate <- function (start_age, end_age) {
  monthly_hospitalization_rate <- numeric(num_months)
  for (i in 1:num_months) {
    monthly_hospitalization_rate[i] <- inf_hosp(raw_data, start_age, end_age, start_date_loop[i], end_date_loop[i])[3]
  }
  return(monthly_hospitalization_rate)
}

# create data frame monthly hospitalization fraction for each age cohort each month of 2020
monthly_hosp_rate_by_age <- data.frame(matrix(ncol = num_age_groups, nrow = num_months))

for (i in 1:num_age_groups) {
  monthly_hosp_rate_by_age[,i] <- monthly_hosp_rate(age_low[i], age_high[i])
}

rownames(monthly_hosp_rate_by_age) <- start_date_loop

age_cohorts <- vector(mode="character", length=num_age_groups)
for (i in 1:(num_age_groups)) {
  age_cohorts[i] <- paste(age_low[i], '-', age_high[i], sep='')
}
age_cohorts[num_age_groups] <-  paste(age_low[num_age_groups], '+', sep='')

colnames(monthly_hosp_rate_by_age) <- age_cohorts

## Plot infection and hospitalization for each age group by month of 2020.
library(ggplot2)
library(reshape2)
library(directlabels)

# add vector for months (by index 1, 2, ..., 11, 12)
monthly_hosp_rate_by_age$Months <- as.Date(start_date_loop)

# melt data frames to plot with ggplot
hosp_rate_melt <- melt(monthly_hosp_rate_by_age, id.vars="Months")

# plot with ggplot
hosp_rate_plot <- ggplot(hosp_rate_melt, aes(Months,value, col=variable)) + geom_point() + geom_line() + scale_x_date(date_breaks = "1 month", date_labels = "%b %y", limits = as.Date(c('2020-03-01','2021-02-01'))) + xlab('Month') + ylab('Case-Hospitalization Fraction') + labs(col = "Age Group") + theme(text = element_text(size = 10)) + scale_y_continuous(labels = scales::percent, limits=c(0,0.5))
print(hosp_rate_plot)
ggsave('Case_Hospitalization_Fraction.png', hosp_rate_plot)

write.csv(monthly_hosp_rate_by_age, 'pathway.csv', row.names = FALSE) # specify pathway to save csv file

