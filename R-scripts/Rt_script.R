
##############################################################
## Calculate variant-specific effective reproduction number ##
##                                                          ##
## NOTE: this script is designed to be compatible with US   ##
##       data collated in the JHU CSSE database             ##
##                                                          ##
## Created by: Jessica Rothman                              ##
## Modified by: Mary Petrone                                ##
##############################################################

## (i) Prelimiary steps ##

# Load libraries
library(tidyverse)
library(EpiEstim)
library(zoo)
library(DescTools)

# set seed
set.seed(1234)

# set working directory
wd = "" # <-- update here
setwd(wd)

# read in function document
source("Rt_functions.R")

## (1) Tabulate and format incident cases ##
# NB use case counts collated by the JHU CSSE dashboard
filename = "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv"
study_start = '' # specify dates as m/d/yy
study_end = '' # specify dates as m/d/yy
region = "" # (optional) Any US state or territory
FIPS = c() # (optional) county-level FIPS code list

# The function get_incident returns dataframe with incident cases by week for the specified region or FIPS
weekly_cases = get_incident(fn = filename, start = study_start, end = study_end, region = NA, FIPS = FIPS)
head(weekly_cases) # visually inspect data

## (2) Estimate cases attributable to variants of interest ##
genome_counts_fn = "" # <-- update here
n_variants = 3 # number of variants you have data for

# calculate proportions of each variant
variant_props = get_variant_props(genome_counts_fn, n_variants)
tail(variant_props)

# Calculate number of cases attributed to each variant category
for(i in 1:n_variants){
  col_name = paste0("var",i,"_cases")
  cases = weekly_cases$Cases*variant_props[n_variants+2+i]
  
  weekly_cases[col_name] = cases
}

weekly_cases = cbind(weekly_cases, variant_props[,2:(n_variants+2)])

weekly_tot = format_data(weekly_cases) # need to rep each datapoint 7 times to mimic week
names(weekly_tot) = names(weekly_cases)
weekly_tot$Week = as.Date(weekly_tot$Week, "%Y-%m-%d")
tail(weekly_tot)

# Calculate a 7 day rolling frequency average for each variant
rolling_vars = c("Week", "Cases", paste0("var", 1:n_variants, "_rolling"), "total_rolling")
rolling_cols = ncol(weekly_tot) - seq(n_variants, 0, -1)

rolling_cases = as.data.frame(cbind(weekly_tot$Week, weekly_tot$Cases, sapply(X = weekly_tot[rolling_cols], FUN = zoo::rollmean, k = 7, fill = NA)))
names(rolling_cases) = rolling_vars
rolling_cases$Week = as.Date(rolling_cases$Week, "1970-01-01")

# Set temporal limits for jeffreys estimates
# This is necessary because when we compute the 7 day rolling averages we don't have estimates for our first 3 days and last 3 days
upperlimit = nrow(weekly_tot)-3
lowerlimit = 4

# set remaining parameters for estimates
rolling_df = rolling_cases
sliding = 21 # sliding window for Rt estimates
var_names = paste0("var", 1:n_variants, "_rolling")
count = "total_rolling"
date = "Week"
cases = "Cases"

jeff_list = Rt_list = smoothed_list = list()
for(i in 1:n_variants){
  
  print(paste0("Calculating estimates for ", var_names[i]))
  
  # get jeffreys intervals
  jeff_list[[i]] = get_jeffreys(rolling_df, upperlimit, lowerlimit, var = var_names[i], count, date)
  
  # estimate Rt
  Rt_list[[i]] = get_mean_Rt(jeffreys_mean = jeff_list[[i]], window = sliding)
  
  # smooth Rt estimates
  smoothed_list[[i]] = smooth_Rt(Rt_list[[i]], rolling_df, lowerlimit, date)
}

# collect mean Rt estimates for each variant
Rt_means = as.data.frame(matrix(NA, nrow = nrow(Rt_list[[1]]$R), ncol = n_variants))
smoothed = as.data.frame(matrix(NA, nrow = nrow(smoothed_list[[1]]), ncol = n_variants))
for(i in 1:n_variants){
  Rt_means[,i] = Rt_list[[i]]$R$`Mean(R)` 
  smoothed[,i] = smoothed_list[[i]]$`smooth_spline_mean$y`
}
dates = seq(rolling_cases$Week[lowerlimit + sliding], rolling_cases$Week[upperlimit] , 1)

Rt_means = cbind(dates, Rt_means)
head(Rt_means)

smoothed_df = cbind(dates, smoothed)

# write to file
Rt_out = "" # <-- update here
write.csv(Rt_means, Rt_out)

smooth_out = "" # <-- update here
write.csv(smoothed_df, smooth_out)

