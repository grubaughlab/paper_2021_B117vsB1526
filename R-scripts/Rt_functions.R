
#########################################################
## Functions for estimating variant-specific effective ##
##  reproduction numbers                               ##
## Created by: Mary Petrone                            ##
## Date: October 28, 2021                              ##
#########################################################


###################################################################################
## get_incident: reads a JHU CSSE timeseries (for US only!!) and returns         ##
##               dataframe of weekly incident cases for a desired region         ##
## variables: fn = path to JHU CSSE github repository, start = study start date  ##
##            end = study end date, region = US state or province,               ##
##            FIPS = county-level FIPS (single value or vector of values)        ##
###################################################################################

get_incident = function(fn, start, end, region = NA, FIPS = NA){
  
  # read in file
  df = read_csv(fn)
  
  # sum cases for entirety of study period by sub-region (e.g. county)
  long = gather(df, Date, Cumulative, all_of(start):all_of(end), factor_key=TRUE)
  
  # clean up data
  long_clean = long[!(grepl("Out of", long$Admin2)) & long$Admin2 != "Unassigned",]
  
  # filter dataset by region of interest
  long_clean_reg = long_clean[long_clean$Province_State == region | long_clean$FIPS %in% FIPS, c("Admin2", "Date", "Cumulative")]
  
  # tabulate cases for state by date
  sum_region = long_clean_reg %>% group_by(Date) %>% summarize(Cumulative = sum(Cumulative))
  
  incident_reg = sum_region %>% mutate(Incident = Cumulative - lag(Cumulative))
  incident_reg = incident_reg[!is.na(incident_reg$Incident),] # remove first data point
  incident_reg$Incident[incident_reg$Incident < 0] = 0 # remove reporting anomalies
  
  # create 'Week' variable
  incident_reg$Date = as.Date(incident_reg$Date, "%m/%d/%y")
  incident_reg$Week = cut.Date(incident_reg$Date, "week")
  
  # tabulate incident cases by week
  incident_weekly = incident_reg %>% group_by(Week) %>% summarize(Cases = sum(Incident))
  
  return(incident_weekly)
  
}



##############################################################################
## get_variant_props: reads a csv containing weekly variant counts and      ##
##                    returns dataframe with weekly variant proportions     ##
## variables: counts_fn = csv with weekly counts of individual variants,    ##
##            n_variants = number of variants to analzye                    ##
## notes: counts file must be a CSV with the first column containing the    ##
##        week of sample collection                                         ##
##############################################################################

get_variant_props = function(counts_fn, n_variants){
  
  # read in file with genome counts
  genome_counts = read.csv(counts_fn)
  genome_props = genome_counts
  
  # tabulate total number of genomes sequenced per week
  genome_props$total = rowSums(genome_props[,c(2:ncol(genome_props))])
  
  # calculate proportions of each genome
  for(i in 1:n_variants){
    col_name = paste0("var",i,"_prop")
    prop = genome_props[,c(i+1)]/genome_props$total
    genome_props[col_name] = prop
  }
  
  return(genome_props)
}



##############################################################################
## format_data: takes dataframe of weekly data and replicates each row      ##
##              7 times; returns formatted dataframe                        ##
## variables: df = dataframe to format                                      ##
##############################################################################

format_data = function(df){
  df_list = list()
  
  for(i in 1:nrow(df)){
    df_list[[i]] = as.data.frame(matrix(rep(df[i,], 7), nrow = 7, byrow = TRUE))
  }
  
  for(j in 1:length(df_list)){
    for(i in 1:ncol(df_list[[j]])){
      df_list[[j]][,i] = unlist(df_list[[j]][,i])
    }
  }
  
  new_df = df_list[[1]]
  for(i in 2:length(df_list)){
    new_df = rbind(new_df, df_list[[i]])
  }
  
  return(new_df)
}

#######################################################################################
## get_jeffreys: calculates jeffreys' intervals for individual variants              ##
## variables: df = dataframe of rolling means for individual variant counts,         ##
##            upperlimit = upper date bound for analysis,                            ##
##            lowerlimit = lower date bound for analysis,                            ##
##            var = column name in df containing variant-specific counts,            ##
##            count = column containing total number of genomes sequenced per week,  ##
##            date = name of column in df containing the week of collection          ##
## notes: For our analysis we will only use the mean from the Jeffrey's              ##
##        interval but you could also use the lower and upper limits.                ##
#######################################################################################

get_jeffreys = function(df, upperlimit, lowerlimit, var, count, date){

  print("Calculating Jeffrey's interval...")
  jeff_estim = BinomCI(x=df[c(lowerlimit:upperlimit), var], n = df[c(lowerlimit:upperlimit), count], 
                       conf.level = 0.95, sides = "two.sided", method = "jeffreys")

  print("Setting dates...")
  day = seq(df[lowerlimit, date], df[upperlimit, date], by = "day")
  jeffreys = cbind.data.frame(day, jeff_estim)
  
  #This multiplies the mean proportion of variant cases calculated in the Jeffrey's interval by the 7 day rolling average of New Cases
  jeffreys_mean_cases = jeffreys$est*df[c(lowerlimit:upperlimit), cases]
  
  return(jeffreys_mean_cases)
}



###########################################################################
## get_mean_Rt: calculates mean Rts for individual variants              ##
## variables: jeffreys_mean = output of get_jeffreys() function,         ##
##            window = size of sliding window (in days) to use           ##
###########################################################################

get_mean_Rt = function(jeffreys_mean, window){

  t_start = seq(2,length(jeffreys_mean) - (window - 1))
  t_end = t_start + (window-1)
  
  #This sets up our Rt estimate. Since we are uncertain about the serial interval, 
  # this allows it to explore multiple serial interval distributions.
  print("Configuring estimates...")
  config = make_config(list(mean_si = 5.2, std_mean_si = 1,min_mean_si = 2.2, max_mean_si = 8.2,
                            std_si = 4, std_std_si = 0.5,min_std_si = 2.5, max_std_si = 5.5,
                            n1 = 500, n2 = 50, t_start = t_start, t_end = t_end)) 
  
  # Generate Rt estimates
  print("Estimating Rt...")
  mean_Rt = estimate_R(incid = jeffreys_mean, method="uncertain_si", config = config)

  return(mean_Rt)
}


#########################################################################################
## smooth_Rt: smooths mean Rt estimates for individual variants                        ##
## variables: Rt_df = output of mean_Rt, rolling_df = dataframe containg rolling       ##
##              means of weekly variant frequencies,                                   ##
##              lowerlimit = lower date bound for analysis,                            ##
##              date = name of column in rolling_df containing the week of collection  ##
#########################################################################################

smooth_Rt = function(Rt_df, rolling_df, lowerlimit, date){
  # Smooth the Rt estimates using a smoothing spline
  print("Smoothing...")
  smooth_spline_mean = with(Rt_df$R, smooth.spline(Rt_df$R$`t_end`, Rt_df$R$`Mean(R)`, cv = TRUE))
  
  smooth_dates_start = rolling_df$Week[c((lowerlimit + Rt_df$R$`t_end`[1]-1))]
  smooth_dates_end = smooth_dates_start + (nrow(Rt_df$R)-1)
  smooth_dates = seq(smooth_dates_start, smooth_dates_end, 1)
  
  smooth_spline_mean_df = cbind.data.frame(smooth_dates, smooth_spline_mean$x, smooth_spline_mean$y)
  
  return(smooth_spline_mean_df)
}







