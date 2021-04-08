## Calculate infection ascertainment from infection- and case- hospitalization fraction

rm(list = ls())
library(ggplot2)
library(reshape2)
library(lubridate)


# infection-hosp fraction (UL, mean, LL) with seroconversion offset of 14 days
inf_hosp_fraction_group1 <- c(0.007070856, 0.005303142, 0.004303999)
inf_hosp_fraction_group2 <- c(0.02721318, 0.02040988, 0.01656454)
inf_hosp_fraction_group3 <- c(0.08711692, 0.06533769, 0.05302769)

# case-hosp fraction for each month
case_hosp_fraction_group1 <- read.csv('monthly_case_hosp_fraction_group1.csv') # pathway to read csv of monthly case-hospitalization fraction
case_hosp_fraction_group2 <- read.csv('monthly_case_hosp_fraction_group2.csv') # pathway to read csv of monthly case-hospitalization fraction
case_hosp_fraction_group3 <- read.csv('monthly_case_hosp_fraction_group3.csv') # pathway to read csv of monthly case-hospitalization fraction

# infection ascertainment
inf_asc_group1_LL <- inf_hosp_fraction_group1[1]/case_hosp_fraction_group1[,2]
inf_asc_group1_mean <- inf_hosp_fraction_group1[2]/case_hosp_fraction_group1[,2]
inf_asc_group1_UL <- inf_hosp_fraction_group1[3]/case_hosp_fraction_group1[,2]

inf_asc_group2_LL <- inf_hosp_fraction_group2[1]/case_hosp_fraction_group2[,2]
inf_asc_group2_mean <- inf_hosp_fraction_group2[2]/case_hosp_fraction_group2[,2]
inf_asc_group2_UL <- inf_hosp_fraction_group2[3]/case_hosp_fraction_group2[,2]

inf_asc_group3_LL <- inf_hosp_fraction_group3[1]/case_hosp_fraction_group3[,2]
inf_asc_group3_mean <- inf_hosp_fraction_group3[2]/case_hosp_fraction_group3[,2]
inf_asc_group3_UL <- inf_hosp_fraction_group3[3]/case_hosp_fraction_group3[,2]

inf_asc_summ_group1 <- data.frame(inf_asc_group1_LL, inf_asc_group1_mean, inf_asc_group1_UL)
inf_asc_summ_group2 <- data.frame(inf_asc_group2_LL, inf_asc_group2_mean, inf_asc_group2_UL)
inf_asc_summ_group3 <- data.frame(inf_asc_group3_LL, inf_asc_group3_mean, inf_asc_group3_UL)

months <- as.Date(case_hosp_fraction_group1$Months)

mean_summ <- data.frame(months, inf_asc_group1_mean, inf_asc_group2_mean, inf_asc_group3_mean)
colnames(mean_summ) <- c('Months', '0-49', '50-69', '70+')

# # test scaling the AF
# AF_offset <- 0
# mean_summ[,c(2,3,4)] <- mean_summ[,c(2,3,4)]+AF_offset

# plot with ggplot
inf_asc_melt <- melt(mean_summ, id.vars="Months")
inf_asc_melt$variable <- as.factor(inf_asc_melt$variable)
inf_asc_melt$UL <- inf_asc_melt$value*0.3333333
inf_asc_melt$LL <- inf_asc_melt$value*0.1884058

inf_asc_plot <- ggplot(inf_asc_melt, aes(Months,value, col=variable)) + geom_point() + geom_line() + scale_x_date(date_breaks = "1 month", date_labels = "%b %y", limits = as.Date(c('2020-03-01','2021-02-01'))) + xlab('Month') + ylab('Infection Ascertainment Fraction') + labs(col = "Age Group") + theme(text = element_text(size = 10)) + scale_y_continuous(labels = scales::percent)
print(inf_asc_plot)

ggsave("inf_asc_true.png", inf_asc_plot)

## Observed epidemic formatting
ObsEpi_processed <- read.csv('pathway_daily.csv') # pathway where csv file with daily reported case counts 

# ObsEpi_processed <- subset(ObsEpi, as.Date(ObsEpi$Dates) >= as.Date('2020-03-01') & as.Date(ObsEpi$Dates) <= as.Date('2021-02-28'))

ObsEpi_processed$Month_Yr <- format(as.Date(ObsEpi_processed$Dates), "%Y-%m")
mean_summ$Month_Yr <- format(as.Date(mean_summ$Months), "%Y-%m")

ObsEpi_group1 <- ObsEpi_processed[,c(1,2,6)]
ObsEpi_group2 <- ObsEpi_processed[,c(1,3,6)]
ObsEpi_group3 <- ObsEpi_processed[,c(1,4,6)]

mean_summ_group1 <- mean_summ[,c(1,2,5)]
mean_summ_group2 <- mean_summ[,c(1,3,5)]
mean_summ_group3 <- mean_summ[,c(1,4,5)]

ObsEpi_group1$factor <- numeric(nrow(ObsEpi_group1))
ObsEpi_group2$factor <- numeric(nrow(ObsEpi_group2))
ObsEpi_group3$factor <- numeric(nrow(ObsEpi_group3))


for (i in 1:nrow(ObsEpi_group1)) {
  for (j in 1:nrow(mean_summ_group1)) {
    if (mean_summ_group1$Month_Yr[j]==ObsEpi_group1$Month_Yr[i]) {
      ObsEpi_group1$factor[i] <- mean_summ_group1[j,2]
    }
  }
}

for (i in 1:nrow(ObsEpi_group2)) {
  for (j in 1:nrow(mean_summ_group2)) {
    if (mean_summ_group2$Month_Yr[j]==ObsEpi_group2$Month_Yr[i]) {
      ObsEpi_group2$factor[i] <- mean_summ_group2[j,2]
    }
  }
}

for (i in 1:nrow(ObsEpi_group3)) {
  for (j in 1:nrow(mean_summ_group3)) {
    if (mean_summ_group3$Month_Yr[j]==ObsEpi_group3$Month_Yr[i]) {
      ObsEpi_group3$factor[i] <- mean_summ_group3[j,2]
    }
  }
}

Dates_vector <- as.Date(ObsEpi_group1$Dates)

ObsEpi_group1$true <- ObsEpi_group1[,2] / ObsEpi_group1[,4]
ObsEpi_group2$true <- ObsEpi_group2[,2] / ObsEpi_group2[,4]
ObsEpi_group3$true <- ObsEpi_group3[,2] / ObsEpi_group3[,4]

true_epi <- data.frame(Dates_vector, ObsEpi_group1$true, ObsEpi_group2$true, ObsEpi_group3$true)


for (i in 1:nrow(true_epi)) {
  for (j in 1:ncol(true_epi)) {
    if (is.na(true_epi[i,j])) {
      true_epi[i,j] <- 0
    }
  }
}

true_epi$Combined <- true_epi[,2] + true_epi[,3] + true_epi[,4] 
colnames(true_epi) <- c('Dates', 'Group1', 'Group2', 'Group3', 'Combined')

write.csv(true_epi, 'TrueEpi.csv') # pathway to save csv


### Plot

# 2n+1 day moving avg function
moving_avg <- function(data_vector, n) {
  len <- length(data_vector)
  avg <- numeric(len)
  for (i in c(1:n, (len-n+1):len)) {
    avg[i] <- NA
  }
  for (i in (n+1):(len-n)) {
    avg[i] <- mean(data_vector[(i-n):(i+n)])
  }
  return(avg)
}

df_all <- data.frame(Dates_vector, ObsEpi_processed$Combined, moving_avg(ObsEpi_processed$Combined, 3), true_epi$Combined, moving_avg(true_epi$Combined, 3))
df <- data.frame(Dates_vector, moving_avg(ObsEpi_processed$Combined, 3), moving_avg(true_epi$Combined, 3))
colnames(df) <- c('Date', 'Observed Epidemic (7-Day Avg)', 'True Epidemic (7-Day Avg)')

# plot with ggplot
df_melt <- melt(df, id.vars="Date")
df_melt$UL <- numeric(nrow(df_melt))
df_melt$LL <- numeric(nrow(df_melt))

for (i in 1:nrow(df_melt)) {
  if(df_melt$variable[i] == 'True Epidemic (7-Day Avg)') {
    df_melt$UL[i] <- df_melt$value[i]*(1+0.1884058)
    df_melt$LL[i] <- df_melt$value[i]*(1-0.3333333)
  }
}


df_plot <- ggplot(df_melt, aes(Date,value, col=variable)) + 
  geom_line() + 
  scale_x_date(date_breaks = "1 month", date_labels = "%b %y", limits = as.Date(c('2020-03-01','2021-02-01'))) + 
  labs(x='Date', y='Daily Cases',col = '') + theme(text = element_text(size = 10)) +  
  geom_ribbon(aes(ymin=LL, ymax=UL), alpha=0.3, linetype = 0,) 


print(df_plot)
ggsave("daily_cases_true.png", df_plot) # pathway to save png



df$True_LL <- df$`True Epidemic (7-Day Avg)`*(1-0.3333333)
df$True_UL <- df$`True Epidemic (7-Day Avg)`*(1+0.1884058)

df_cumul <- data.frame(ObsEpi_processed$Dates, ObsEpi_processed$Combined, true_epi$Combined)
colnames(df_cumul) <- c('Date', 'Observed Epidemic (Cumulative)', 'True Epidemic (Cumulative)')
df_cumul$Date <- as.Date(df_cumul$Date)


for (i in 2:nrow(df_cumul)) {
  for (j in 2:ncol(df_cumul)) {
    df_cumul[i,j] <- df_cumul[i,j] + df_cumul[i-1,j]
  }
}



# plot with ggplot
df_cumul_melt <- melt(df_cumul, id.vars="Date")

df_cumul_melt$UL <- numeric(nrow(df_melt))
df_cumul_melt$LL <- numeric(nrow(df_melt))

for (i in 1:nrow(df_cumul_melt)) {
  if(df_cumul_melt$variable[i] == 'True Epidemic (Cumulative)') {
    df_cumul_melt$UL[i] <- df_cumul_melt$value[i]*(1+0.1884058)
    df_cumul_melt$LL[i] <- df_cumul_melt$value[i]*(1-0.3333333)
  }
}

df_cumul_melt$value <- df_cumul_melt$value/5071336
df_cumul_melt$UL <- df_cumul_melt$UL/5071336
df_cumul_melt$LL <- df_cumul_melt$LL/5071336

df_cumul_plot <- ggplot(df_cumul_melt, aes(Date,value, col=variable)) + 
  geom_line() + 
  scale_y_continuous(labels = scales::percent) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %y", limits = as.Date(c('2020-03-01','2021-02-01'))) + 
  labs(x='Date', y='Proportion Infected',col = '') + theme(text = element_text(size = 10)) +  
  geom_ribbon(aes(ymin=LL, ymax=UL), alpha=0.3, linetype = 0,) +
  geom_rect(aes(xmin=as.Date('2020-05-09')-14, xmax=as.Date('2020-07-21')-14, ymin=0.0042, ymax=0.0069), fill=NA, color='red',size=0.25)+
  geom_rect(aes(xmin=as.Date('2020-05-09')-14, xmax=as.Date('2020-07-21')-14, ymin=0.0056, ymax=0.0056), fill=NA, color='black', size=0.25)
  # geom_rect(aes(xmin=as.Date('2020-09-18')-14, xmax=as.Date('2020-09-29')-14, ymin=0.00099, ymax=0.01955), fill=NA, color='red',size=0.25)+
  # geom_rect(aes(xmin=as.Date('2020-09-18')-14, xmax=as.Date('2020-09-29')-14, ymin=0.00808, ymax=0.00808), fill=NA, color='black', size=0.25)


print(df_cumul_plot)
ggsave("cumul_cases_true_figure.png", df_cumul_plot) # pathway to save png

