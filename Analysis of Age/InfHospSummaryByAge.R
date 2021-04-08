## Objective: summarize per capita infection and hospitalization data over course of epidemic in 2020.
## Last updated: 2021-03-07.

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
num_age_groups <- length(pop_age)
age_low <- numeric(num_age_groups)
age_high <- numeric(num_age_groups)

for (i in 1:num_age_groups) {
  age_low[i] <- 10*(i-1)
  age_high[i] <- 10*i-1
}

# create vectors describing start and end dates for each month of 2020
num_months <- 12
start_date_loop <- vector(mode="character", length=num_months)
end_date_loop <- vector(mode="character", length=num_months)

for (i in 1:num_months) {
  start_date_loop[i] <- paste('2020-', i, '-01', sep='')
  end_date_loop[i] <- paste('2020-', i, '-', as.character(days_in_month(i)), sep='')
}

# functions calculating infection and hospitalization count for age cohort of interest each month of 2020 
monthly_inf <- function (start_age, end_age) {
  monthly_infections <- numeric(num_months)
  for (i in 1:num_months) {
    monthly_infections[i] <- inf_hosp(raw_data, start_age, end_age, start_date_loop[i], end_date_loop[i])[1]
  }
  return(monthly_infections)
}

monthly_hosp <- function (start_age, end_age) {
  monthly_hospitalizations <- numeric(num_months)
  for (i in 1:num_months) {
    monthly_hospitalizations[i] <- inf_hosp(raw_data, start_age, end_age, start_date_loop[i], end_date_loop[i])[2]
  }
  return(monthly_hospitalizations)
}

monthly_hosp_rate <- function (start_age, end_age) {
  monthly_hospitalization_rate <- numeric(num_months)
  for (i in 1:num_months) {
    monthly_hospitalization_rate[i] <- inf_hosp(raw_data, start_age, end_age, start_date_loop[i], end_date_loop[i])[3]
  }
  return(monthly_hospitalization_rate)
}

# create data frames summarizing per capita infections and hospitalizations for each age cohort each month of 2020
monthly_inf_by_age <- data.frame(matrix(ncol = num_age_groups+1, nrow = num_months))
monthly_hosp_by_age <- data.frame(matrix(ncol = num_age_groups+1, nrow = num_months))
monthly_hosp_rate_by_age <- data.frame(matrix(ncol = num_age_groups+1, nrow = num_months))

for (i in 1:num_age_groups) {
  monthly_inf_by_age[,i] <- monthly_inf(age_low[i], age_high[i])/pop_age[i]
  monthly_hosp_by_age[,i] <- monthly_hosp(age_low[i], age_high[i])/pop_age[i]
  monthly_hosp_rate_by_age[,i] <- monthly_hosp_rate(age_low[i], age_high[i])
}
monthly_inf_by_age[,(num_age_groups+1)] <- monthly_inf(0, 150)/sum(pop_age)
monthly_hosp_by_age[,(num_age_groups+1)] <- monthly_hosp(0, 150)/sum(pop_age)
monthly_hosp_rate_by_age[,(num_age_groups+1)] <- monthly_hosp_rate(0, 150)

rownames(monthly_inf_by_age) <- month.abb[1:num_months]
rownames(monthly_hosp_by_age) <- month.abb[1:num_months]
rownames(monthly_hosp_rate_by_age) <- month.abb[1:num_months]

age_cohorts <- vector(mode="character", length=num_age_groups)
for (i in 1:(num_age_groups)) {
  age_cohorts[i] <- paste(age_low[i], '-', age_high[i], sep='')
}
age_cohorts[num_age_groups] <-  paste(age_low[num_age_groups], '+', sep='')
age_cohorts[(num_age_groups+1)] <-  'All'
  
colnames(monthly_inf_by_age) <- age_cohorts
colnames(monthly_hosp_by_age) <- age_cohorts
colnames(monthly_hosp_rate_by_age) <- age_cohorts

## Plot infection and hospitalization for each age group by month of 2020.
library(ggplot2)
library(reshape2)
library(directlabels)

# add vector for months (by index 1, 2, ..., 11, 12)
monthly_inf_by_age$Months <- as.Date(start_date_loop)
monthly_hosp_by_age$Months <- as.Date(start_date_loop)
monthly_hosp_rate_by_age$Months <- as.Date(start_date_loop)

# melt data frames to plot with ggplot
inf_melt <- melt(monthly_inf_by_age, id.vars="Months")
hosp_melt <- melt(monthly_hosp_by_age, id.vars="Months")
hosp_rate_melt <- melt(monthly_hosp_rate_by_age, id.vars="Months")

# plot with ggplot
inf_plot <- ggplot(inf_melt, aes(Months,value, col=variable)) + geom_point() + geom_line() + geom_dl(aes(label = variable), method = list('last.qp', cex=.75, hjust = 1, vjust = 0)) + scale_x_date(date_breaks = "1 month", date_labels = "%b", limits = as.Date(c('2020-01-01','2020-12-01'))) + xlab('Month of 2020') + ylab('Monthly Cases (Per Capita)') + labs(col = "Cohorts") + theme(text = element_text(size = 10)) 
hosp_plot <- ggplot(hosp_melt, aes(Months,value, col=variable)) + geom_point() + geom_line() + geom_dl(aes(label = variable), method = list('last.qp', cex=.75, hjust = 1, vjust = 0)) + scale_x_date(date_breaks = "1 month", date_labels = "%b", limits = as.Date(c('2020-01-01','2020-12-01'))) + xlab('Month of 2020') + ylab('Monthly Hospitalizations (Per Capita)') + labs(col = "Cohorts") + theme(text = element_text(size = 10)) 
hosp_rate_plot <- ggplot(hosp_rate_melt, aes(Months,value, col=variable)) + geom_point() + geom_line() + geom_dl(aes(label = variable), method = list('last.qp', cex=.75, hjust = 1, vjust = 0)) + scale_x_date(date_breaks = "1 month", date_labels = "%b", limits = as.Date(c('2020-03-01','2020-12-01'))) + xlab('Month of 2020') + ylab('Case-Hospitalization Fraction') + labs(col = "Cohorts") + theme(text = element_text(size = 10)) + scale_y_continuous(labels = scales::percent, limits=c(0,0.5))

ggsave('Per_Capita_Cases.png', inf_plot)
ggsave('Per_Capita_Hospitalizations.png', hosp_plot)
ggsave('Case_Hospitalization_Fraction.png', hosp_rate_plot)

# k-means clustering analysis
inf_transpose <- t(monthly_inf_by_age[c(9, 10, 11, 12),-c(10,11,12)])
hosp_transpose <- t(monthly_hosp_by_age[c(9, 10, 11, 12),-c(10,11,12)])
hosp_fraction_transpose <- t(monthly_hosp_rate_by_age[c(9, 10, 11, 12),-c(10,11,12)])

kmean_inf <- kmeans(inf_transpose, 3)
kmean_hosp <- kmeans(hosp_transpose, 3)
kmean_hosp_fraction <- kmeans(hosp_fraction_transpose, 3)
