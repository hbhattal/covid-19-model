## Objective: calculate infection-hospitalization fraction in each age group.

# import 'raw_data' and 'date_pulled' into workspace 

rm(list = ls()[!ls() %in% c("raw_data", "date_pulled")]) 
library(date)

# set parameters
region <- 'BC'
offset <- 21 # how much of recent data to ignore
age_cutoff_1 <- 50
age_cutoff_2 <- 70

# ignore recent data
processed_data <- subset(raw_data, raw_data$reported_date <= (date_pulled-offset))

# subset data from HA of interest
if (region != 'BC') {
  processed_data <- subset(processed_data, processed_data$ha_res == region)
}

# subset hospitalized cases, 
hospitalized_data <- subset(processed_data, processed_data$ever_hospitalized == 'Yes')

# delete entries with incomplete data
hospitalized_data <- subset(hospitalized_data, hospitalized_data$admission_date != '')

# split into three Groups (by age)
hospitalized_data_Group1 <- subset(hospitalized_data, hospitalized_data$age_combined < age_cutoff_1)
hospitalized_data_Group2 <- subset(hospitalized_data, hospitalized_data$age_combined >= age_cutoff_1 & hospitalized_data$age_combined < age_cutoff_2)
hospitalized_data_Group3 <- subset(hospitalized_data, hospitalized_data$age_combined >= age_cutoff_2)

# calculate cumulative hospitalization case data (from first reported case 2020-01-26)
Dates <- as.character(seq(as.Date("2020-01-26"), as.Date(date_pulled-offset), by="days"))

daily_cases <- function(dataset) {
  daily <- numeric(length(Dates))
  for (i in 1:length(Dates)) {
    daily[i] <- length(grep(Dates[i], dataset$admission_date))
  }
  return(daily)
}

cumul_cases <- function(dataset) {
  daily <- daily_cases(dataset)
  cumul <- numeric(length(Dates))
  cumul[1] <- daily[1]
  for (i in 2:length(Dates)) {
    cumul[i] <- cumul[i-1] + daily[i]
  }
  return(cumul)
}

Group1 <- cumul_cases(hospitalized_data_Group1)
Group2 <- cumul_cases(hospitalized_data_Group2)
Group3 <- cumul_cases(hospitalized_data_Group3)
Combined <- Group1 + Group2 + Group3

cumul_cases_summ <- data.frame(Dates, Group1, Group2, Group3, Combined)
colnames(cumul_cases_summ) <- c('Dates', 'Group1', 'Group2', 'Group3','Combined')

### ----------------------------------------------------------------------------

# population data for each of eleven age cohorts (0-9, 10-19, ..., 90+) from 2016 Census
pop_total <- 5071336
pop_age <- c(457525, 492840, 590560, 607340, 617410, 709300, 611615, 347010, 172765, 40360 + 1325) 
pop_age_prop <- pop_age/sum(pop_age)

# function to calculate population in each group
pop_three_group <- function (age_cutoff_1, age_cutoff_2) {
  index1 = floor(age_cutoff_1/10);
  index2 = floor(age_cutoff_2/10);
  
  pop1 = sum(pop_age_prop[1:index1])*pop_total;
  pop2 = sum(pop_age_prop[(index1+1):index2])*pop_total;
  pop3 = sum(pop_age_prop[(index2+1):length(pop_age_prop)])*pop_total;
  
  return(c(pop1, pop2, pop3))
}

# seropositivity of all age groups (Saeed et al., 2021)
sero_mean <- 0.0056
sero_LL <- 0.0042
sero_UL <- 0.0069

# function to calculate estimated cumulative infections in each group 
sero_group <- function(population_size) {
  LL <- sero_LL*population_size
  mean <- sero_mean*population_size
  UL <- sero_UL*population_size
  
  return(c(LL, mean, UL))
}

# use functions to calculate estimate for cumulative infection in each group
pop_50_70 <- pop_three_group(age_cutoff_1, age_cutoff_2)

inf_group1 <- sero_group(pop_50_70[1])
inf_group2 <- sero_group(pop_50_70[2])
inf_group3 <- sero_group(pop_50_70[3])

# subset cumulative cases from May 9 and July 21, 2020 (with 14 day offset for median seroconversion)
sero_offset <- 14 # median days to seroconversion
cumul_cases_study <- subset(cumul_cases_summ, cumul_cases_summ$Dates >= (as.Date('2020-05-09')-sero_offset) & cumul_cases_summ$Dates <= (as.Date('2020-07-21')-sero_offset))
cumul_cases_group1 <- cumul_cases_study[,c(1,2)]
cumul_cases_group2 <- cumul_cases_study[,c(1,3)]
cumul_cases_group3 <- cumul_cases_study[,c(1,4)]

# function to calculate inf-hosp fraction
inf_hosp_fraction <- function(cumul_hosp_dataset, inf_data) {
  num_rows <- nrow(cumul_hosp_dataset)
  
  LL_inf_hosp_fraction <- numeric(num_rows)
  mean_inf_hosp_fraction <- numeric(num_rows)
  UL_inf_hosp_fraction <- numeric(num_rows)
  
  for (i in 1:num_rows) {
    LL_inf_hosp_fraction[i] <- cumul_hosp_dataset[i,2]/inf_data[1]
    mean_inf_hosp_fraction[i] <- cumul_hosp_dataset[i,2]/inf_data[2]
    UL_inf_hosp_fraction[i] <- cumul_hosp_dataset[i,2]/inf_data[3]
  }
  
  cumul_hosp_dataset[,3] <- LL_inf_hosp_fraction
  cumul_hosp_dataset[,4] <- mean_inf_hosp_fraction
  cumul_hosp_dataset[,5] <- UL_inf_hosp_fraction
  
  colnames(cumul_hosp_dataset) <- c('Dates', 'Cumul Hosp', 'UL Inf-Hosp', 'Mean Inf-Hosp', 'LL Inf-Hosp')
  return (cumul_hosp_dataset)
}

# calculate inf-hosp fractions (LL, mean, UL) for each group
inf_hosp_group1 <- inf_hosp_fraction(cumul_cases_group1, inf_group1)
inf_hosp_group2 <- inf_hosp_fraction(cumul_cases_group2, inf_group2)
inf_hosp_group3 <- inf_hosp_fraction(cumul_cases_group3, inf_group3)

# function to weight estimate by number of tests performed (Saeed et al., 2021)
weight <- function (Date_vector) {
  Date_vector <- as.Date(Date_vector)
  len_date_vector <- length(Date_vector)
  weight <- numeric(len_date_vector)
  
  for (i in 1:len_date_vector) {
    if ((as.Date('2020-05-09')-sero_offset) <= Date_vector[i] & Date_vector[i] <= (as.Date('2020-05-23')-sero_offset)) {
      weight[i] <- 12921
    }
    if ((as.Date('2020-05-24')-sero_offset) <= Date_vector[i] & Date_vector[i] <= (as.Date('2020-06-07')-sero_offset)) {
      weight[i] <- 16167
    }
    if ((as.Date('2020-06-08')-sero_offset) <= Date_vector[i] & Date_vector[i] <= (as.Date('2020-06-22')-sero_offset)) {
      weight[i] <- 22492
    }
    if ((as.Date('2020-06-23')-sero_offset) <= Date_vector[i] & Date_vector[i] <= (as.Date('2020-07-07')-sero_offset)) {
      weight[i] <- 18068
    }
    if ((as.Date('2020-07-08')-sero_offset) <= Date_vector[i] & Date_vector[i] <= (as.Date('2020-07-21')-sero_offset)) {
      weight[i] <- 4994
    }
  }
  
  return(weight)
}

weight_all <- weight(cumul_cases_study$Dates)
sum_weight_all <- sum(weight_all)

# function to calculate true inf-hosp fraction based on weight of number of test performed each week
true_inf_hosp <- function (inf_hosp_group) {
  avg_UL <- as.numeric(inf_hosp_group[,3]%*%weight_all)/sum_weight_all
  avg_mean <- as.numeric(inf_hosp_group[,4]%*%weight_all)/sum_weight_all
  avg_LL <- as.numeric(inf_hosp_group[,5]%*%weight_all)/sum_weight_all
  
  return(c(avg_UL, avg_mean, avg_LL))
}

true_inf_hosp_group1 <- true_inf_hosp(inf_hosp_group1)
true_inf_hosp_group2 <- true_inf_hosp(inf_hosp_group2)
true_inf_hosp_group3 <- true_inf_hosp(inf_hosp_group3)

