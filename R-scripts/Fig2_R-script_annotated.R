library(readr)
library(tidyverse)
library(EpiEstim)
library(xlsx)
library(stats)
library(zoo)
library(DescTools)
set.seed(1234)
# Read in case data from JHU CSSE COVID-19 dataset
covid_cases<-read_csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv")

# Read in death data from JHU CSSE COVID-19 dataset
covid_deaths<-read_csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_US.csv")

# Create Connecticut (or other state or country) dataset with cases (new and cumulative)
long_cases<-gather(covid_cases, Date, Cumulative_Cases, '4/1/21':'7/12/21', factor_key=TRUE) #subset to the dates for which you want to calculate Rt
long_CT_cases<-long_cases[which(long_cases$Province_State=="Connecticut"),] #subsets the data to Connecticut cases (can be adjusted to other state or country)
long_CT_cases<-subset(long_CT_cases,Admin2!="Out of CT"&Admin2!= "Unassigned")
long_CT_cases<-aggregate(x=long_CT_cases$Cumulative_Cases,by=list(long_CT_cases$Date),FUN=sum)
names(long_CT_cases)<-c("Date","Cumulative_Cases")
long_CT_cases<-long_CT_cases %>%mutate(New_Cases = Cumulative_Cases - lag(Cumulative_Cases))

# Create Connecticut (or other state or country) dataset with deaths (new and cumulative)
long_deaths<-gather(covid_deaths, Date, Cumulative_Deaths, '4/1/21':'7/12/21', factor_key=TRUE)#subset to the dates for which you want to calculate Rt
long_CT_deaths<-long_deaths[which(long_deaths$Province_State=="Connecticut"),]#subsets the data to Connecticut deaths (can be adjusted to other state or country)
long_CT_deaths<-subset(long_CT_deaths,Admin2!="Out of CT"&Admin2!= "Unassigned")
long_CT_deaths<-aggregate(x=long_CT_deaths$Cumulative_Deaths,by=list(long_CT_deaths$Date),FUN=sum)
names(long_CT_deaths)<-c("Date","Cumulative_Deaths")
long_CT_deaths<-long_CT_deaths %>%mutate(New_Deaths = Cumulative_Deaths - lag(Cumulative_Deaths))

#Combine the case and death data into dataframe
covid_CT_wide<-cbind.data.frame(long_CT_cases$Date, long_CT_cases$Cumulative_Cases,
                                long_CT_cases$New_Cases,
                                long_CT_deaths$Cumulative_Deaths, long_CT_deaths$New_Deaths)
names(covid_CT_wide)<-c("Date", "Cumulative_Cases","New_Cases","Cumulative_Deaths","New_Deaths")
covid_CT_wide[is.na(covid_CT_wide)] <- 0
#Set any new cases or new deaths that are less than 0 to 0
covid_CT_wide$New_Cases<-ifelse(covid_CT_wide$New_Cases<0,0,covid_CT_wide$New_Cases)
covid_CT_wide$New_Deaths<-ifelse(covid_CT_wide$New_Deaths<0,0,covid_CT_wide$New_Deaths)

#This creates a variable called "end_of_week_date". This is not necessary to compute, but can be useful when you only have weekly data for the proportions of each variant
covid_CT_wide$Date<-as.character(covid_CT_wide$Date)
covid_CT_wide$Date<-as.Date(covid_CT_wide$Date,"%m/%d/%y")
end_of_week_date<-seq(as.Date("2021/4/7"), as.Date("2021/7/12"), by = "week") #change the second date to current date (last date of data)
weeklycases<-rollapply(covid_CT_wide$New_Cases, width=7, FUN=sum, by=7)#This gives you the cases per week in your state or country of interest
wide<-cbind.data.frame(end_of_week_date,weeklycases)

#Here you must read in the proportions of each variant. We calculated the proportions by: number of variant samples sequenced in lab per week divided by number of total samples sequence per week
wide$alpha_prop<-c() #read in "alpha" proportions.since we are using daily data you will copy the proportion 7 times per week (7 rows). 
wide$iota_prop<-c() #read in "iota" proportions. since we are using daily data you will copy the proportion 7 times per week (7 rows). 
wide$other_prop<-1- (wide$alpha_prop+wide$iota_prop) 

#Read in the variable that gives the total number of samples sequenced per week. This is the denominator of the variant proportions.
wide$n<-c() #number of samples sequenced per week. since we are using daily data you will copy n 7 times per week (7 rows). 

#This creates a new dataframe with the variables we need for our analysis. 
#Often you will have new case data that extends past your sequencing data. This subsets the dataframe so that it only contains days where you have sequencing data (i.e. variant proportions).
wide$New_Cases<-covid_CT_wide$New_Cases
wide$Date<-covid_CT_wide$Date

wide<-wide %>% 
  select(alpha_prop,
         iota_prop,
         other_prop,
         New_Cases,
         n,
         Date) %>%
  mutate(Date = as.Date(Date, "%m/%d/%y"))%>%
  drop_na

upperlimit<-nrow(wide)-3 #This is necessary because when we compute the 7 day rolling averages we don't have estimates for our first 3 days and last 3 days

#This calculates the number of cases of each variant by multiplying n (total number of samples sequenced per week) by the proportion of each variant per week 
wide$alpha<-wide$n*wide$alpha_prop
wide$iota<-wide$n*wide$iota_prop
wide$other<-wide$n*wide$other_prop


#This creates a 7 day rolling average
daily_7<-wide%>%dplyr::mutate(alpha_7 = zoo::rollmean(wide$alpha, k = 7, fill = NA),
                              iota_7 = zoo::rollmean(wide$iota, k = 7, fill = NA),
                              other_7 = zoo::rollmean(wide$other, k = 7, fill = NA),
                              n_7 = zoo::rollmean(n, k = 7, fill = NA))
daily_7$rolling_avg_alpha<-daily_7$alpha_7/daily_7$n_7
daily_7$rolling_avg_iota<-daily_7$iota_7/daily_7$n_7
daily_7$rolling_avg_other<-daily_7$other_7/daily_7$n_7
rolling_avg<-cbind.data.frame(daily_7$rolling_avg_alpha,daily_7$rolling_avg_iota,daily_7$rolling_avg_other)

alpha_rolling_avg_daily_cases<-daily_7$rolling_avg_alpha*daily_7$New_Cases
iota_rolling_avg_daily_cases<-daily_7$rolling_avg_iota*daily_7$New_Cases
other_rolling_avg_daily_cases<-daily_7$rolling_avg_other*daily_7$New_Cases
daily_cases_rolling_avg<-cbind.data.frame(alpha_rolling_avg_daily_cases,iota_rolling_avg_daily_cases,other_rolling_avg_daily_cases)


####ALPHA####
#This calculates the Jeffrey's interval (a binomial confidence interval) for the variant of interest. For our analysis we will only use the mean from the Jeffrey's interval but you could also use the lower and upper limits.
daily_alpha_confint<-BinomCI(x=daily_7$alpha_7[4:upperlimit], n=daily_7$n_7[4:upperlimit], conf.level = 0.95, sides = "two.sided",method = "jeffreys")
day<-seq(daily_7$Date[4], daily_7$Date[upperlimit],by="day")
alpha_jeffreys<-cbind.data.frame(day,daily_alpha_confint)

t_start<-seq(2,length(alpha_jeffreys$day)-21) #This uses a 3 week sliding interval. Other sliding interval lengths can be used.
t_end<-t_start+21 #This uses a 3 week sliding interval. Other sliding interval lengths can be used.

#This sets up our Rt estimate. Since we are uncertain about the serial interval, this allows it to explore multiple serial interval distributions.
config <- make_config(list(mean_si = 5.2, std_mean_si = 1,min_mean_si = 2.2, max_mean_si = 8.2,
                           std_si = 4, std_std_si = 0.5,min_std_si = 2.5, max_std_si = 5.5,
                           n1=500,n2=50,t_start=t_start, t_end=t_end)) 

#This multiplies the mean proportion of variant cases calculated in the Jeffrey's interval by the 7 day rolling average of New Cases
alpha_jeffreys_cases_mean<-alpha_jeffreys$est*daily_7$New_Cases[4:upperlimit]

#Here we calculate Rt
mean_Rt_alpha<-estimate_R(alpha_jeffreys_cases_mean, 
                          method="uncertain_si",
                          config = config)
#Here we smooth the Rt estimates using a smoothing spline
smooth_spline_alpha_mean<- with(mean_Rt_alpha$R, smooth.spline(mean_Rt_alpha$R$`t_end`, mean_Rt_alpha$R$`Mean(R)`, cv = TRUE))
smooth_spline_alpha_mean_df<-cbind.data.frame(daily_7$Date[26:upperlimit],smooth_spline_alpha_mean$x,smooth_spline_alpha_mean$y)



####IOTA####
#This calculates the Jeffrey's interval (a binomial confidence interval) for the variant of interest. For our analysis we will only use the mean from the Jeffrey's interval but you could also use the lower and upper limits.
daily_iota_confint<-BinomCI(x=daily_7$iota_7[4:upperlimit], n=daily_7$n_7[4:upperlimit], conf.level = 0.95, sides = "two.sided",method = "jeffreys")
day<-seq(daily_7$Date[4], daily_7$Date[upperlimit],by="day")
iota_jeffreys<-cbind.data.frame(day,daily_iota_confint)

t_start<-seq(2,length(iota_jeffreys$day)-21) #This uses a 3 week sliding interval. Other sliding interval lengths can be used.
t_end<-t_start+21 #This uses a 3 week sliding interval. Other sliding interval lengths can be used.

#This sets up our Rt estimate. Since we are uncertain about the serial interval, this allows it to explore multiple serial interval distributions.
config <- make_config(list(mean_si = 5.2, std_mean_si = 1,min_mean_si = 2.2, max_mean_si = 8.2,
                           std_si = 4, std_std_si = 0.5,min_std_si = 2.5, max_std_si = 5.5,
                           n1=500,n2=50,t_start=t_start, t_end=t_end)) 

#This multiplies the mean proportion of variant cases calculated in the Jeffrey's interval by the 7 day rolling average of New Cases
iota_jeffreys_cases_mean<-iota_jeffreys$est*daily_7$New_Cases[4:upperlimit]

#Here we calculate Rt
mean_Rt_iota<-estimate_R(iota_jeffreys_cases_mean, 
                          method="uncertain_si",
                          config = config)
#Here we smooth the Rt estimates using a smoothing spline
smooth_spline_iota_mean<- with(mean_Rt_iota$R, smooth.spline(mean_Rt_iota$R$`t_end`, mean_Rt_iota$R$`Mean(R)`, cv = TRUE))
smooth_spline_iota_mean_df<-cbind.data.frame(daily_7$Date[26:upperlimit],smooth_spline_iota_mean$x,smooth_spline_iota_mean$y)




####OTHER####
#This calculates the Jeffrey's interval (a binomial confidence interval) for the variant of interest. For our analysis we will only use the mean from the Jeffrey's interval but you could also use the lower and upper limits.
daily_other_confint<-BinomCI(x=daily_7$other_7[4:upperlimit], n=daily_7$n_7[4:upperlimit], conf.level = 0.95, sides = "two.sided",method = "jeffreys")
day<-seq(daily_7$Date[4], daily_7$Date[upperlimit],by="day")
other_jeffreys<-cbind.data.frame(day,daily_other_confint)

t_start<-seq(2,length(other_jeffreys$day)-21) #This uses a 3 week sliding interval. Other sliding interval lengths can be used.
t_end<-t_start+21 #This uses a 3 week sliding interval. Other sliding interval lengths can be used.

#This sets up our Rt estimate. Since we are uncertain about the serial interval, this allows it to explore multiple serial interval distributions.
config <- make_config(list(mean_si = 5.2, std_mean_si = 1,min_mean_si = 2.2, max_mean_si = 8.2,
                           std_si = 4, std_std_si = 0.5,min_std_si = 2.5, max_std_si = 5.5,
                           n1=500,n2=50,t_start=t_start, t_end=t_end)) 

#This multiplies the mean proportion of variant cases calculated in the Jeffrey's interval by the 7 day rolling average of New Cases
other_jeffreys_cases_mean<-other_jeffreys$est*daily_7$New_Cases[4:upperlimit]

#Here we calculate Rt
mean_Rt_other<-estimate_R(other_jeffreys_cases_mean, 
                          method="uncertain_si",
                          config = config)
#Here we smooth the Rt estimates using a smoothing spline
smooth_spline_other_mean<- with(mean_Rt_other$R, smooth.spline(mean_Rt_other$R$`t_end`, mean_Rt_other$R$`Mean(R)`, cv = TRUE))
smooth_spline_other_mean_df<-cbind.data.frame(daily_7$Date[26:upperlimit],smooth_spline_other_mean$x,smooth_spline_other_mean$y)


