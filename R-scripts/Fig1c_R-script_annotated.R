
###########################
## Figure 1c            ##
## Get daily frequencies ##
###########################

library(tidyverse)

# Read in genomic metadata with (1) collection date, (2) location (e.g. country, state, county), (3) PANGO lineage
# Column headings used here: date, state, county, lineage

wd = "" # set working directly
setwd(wd)

genome_meta = "" # CSV file with metadata
genomes = read.csv(genome_meta)

date_format = "" # date format from CSV (e.g. %m/%d/%y)
genomes$date = as.Date(genomes$date, format)

# If applicable, select genomes from specific geographic location
# For this paper, we first assessed frequency in the state of Connecticut (CT)
genomes_CT = genomes[genomes$state == "Connecticut",]

# If applicable, create a general lineage classification
# For this paper, we created 3 categories: B.1.526*, B.1.1.7, and other
genomes_CT$gen.lin = genomes_CT$lineage

# Note: we did not sequence any B.1.526.3 genomes in Connecticut during our study period
genomes_CT$gen.lin[genomes_CT$gen.lin == "B.1.526" | genomes_CT$gen.lin == "B.1.526.1" | 
                     genomes_CT$gen.lin == "B.1.526.2"] = "B.1.526*"
genomes_CT$gen.lin[genomes_CT$gen.lin != "B.1.1.7" & genomes_CT$gen.lin != "B.1.526*"] = "other"

# Tabulate general lineages by date (can be daily or weekly)
daily_B117 = genomes_CT %>% group_by(date) %>% summarize(B117 = sum(gen.lin == "B.1.1.7"))
daily_B1526 = genomes_CT %>% group_by(date) %>% summarize(B1526 = sum(gen.lin == "B.1.526*"))
daily_other = genomes_CT %>% group_by(date) %>% summarize(other = sum(gen.lin == "other"))

# combine tabulations into a single dataframe and write to file
daily_tot = cbind(daily_B117, daily_B1526, daily_other)

daily_out = "" # name of CSV file where data will be written
write.csv(daily_tot, daily_out)


### If applicable, repeat process for counties ###
genomes_county = genomes_CT[genomes_CT$county == "New Haven",]
daily_B117 = genomes_county %>% group_by(date) %>% summarize(B117 = sum(gen.lin == "B.1.1.7"))
daily_B1526 = genomes_county %>% group_by(date) %>% summarize(B1526 = sum(gen.lin == "B.1.526*"))
daily_other = genomes_county %>% group_by(date) %>% summarize(other = sum(gen.lin == "other"))

daily_county = cbind(daily_B117, daily_B1526, daily_other)

daily_county_out = "" # name of CSV file where data will be written
write.csv(daily_county, daily_county_out)

## For county-level frequencies, we used the 'daily_county' output 
## as the input for our 7-day rolling average estimates (Fig. 1c)
