# Filters MASTREE+ data to units of interest and calculates summary statistics and masting metrics
# pooled version combines fruits/individual and seeds/individual
#generates different versions of the dataset with different duration thresholds of time series 3-6 years
# Written by Jessie Foest, with modifications by Esther Dale and spatial synchrony function by Andrew Tanentzap, Feb 2021

# Data / functions ----
library(tidyverse)
library(here)
library(forecast)
library(geosphere)
library(maditr)

## ACF function ====
any_acf <- function(ts, lag = 1, partial = FALSE) {
  
  # calculate the correlations
  acf <-  Acf(ts, plot = FALSE)
  pacf <- Pacf(ts, plot = FALSE)
  # ACF (partial = FALSE)
  ifelse(
    partial == FALSE,
    # extract the acf value
    acf[[1]][{{lag}}+1] # +1 deals with the fact that Acf starts at lag 0
    ,
    # PACF (partial = TRUE)
    # extract the pacf value
    pacf[[1]][{{lag}}]
  )
}

## Spatial Synchrony Function ====
spat_syn <- function(data_tibble, species, scale_km, global, min_num_yrs) {
  if(global==T){scale_km <- 1e6} # set to exceed circumfrence of Earth
  # subset the table
  mastree <- data_tibble
  species_sub <- which(mastree$Species == species)
  
  # extract the unique coordinates in the table and calculate a distance matrix
  unique_site_coords <- distinct(mastree[species_sub,c('Longitude','Latitude')])
  site_d_mat <- distm(unique_site_coords,fun=distGeo)/1000  # covert to km distance
  
  # for each unique set of coordinates, find other sites within scale_km distance
  # keep only lower triangle of distance matrix to avoid duplication, i.e. mirroring coordinate combinations
  if(nrow(site_d_mat)>1){
    site_pairs <- sapply(1:(nrow(site_d_mat)-1), function(x) { which(site_d_mat[(x+1):nrow(site_d_mat),x] <= scale_km)+x}) # add 1 because focal site x is incorrectly taken as zero-point
    # find coordinates that are within scale_km of another set of coordinates
    sites_with_matches <- which(sapply(site_pairs, length) > 0)
  } else { sites_with_matches <- 1; site_pairs <- NA }
  
  all_cors <- unlist(sapply(sites_with_matches, function(y){
    # subset rows from each unqiue site that has other sites (matches) within scale_km 
    # first find time series with the same focal coordinates (more than 1 "site" can exist at a coordinate) 
    site_x1 <- mastree[species_sub,][apply(mastree[species_sub,c('Longitude','Latitude')],1,identical,unlist(unique_site_coords[y,])),]
    # then find time series that are within scale_km of the focal coordinates      
    site_y1 <- merge(mastree[species_sub,],unique_site_coords[site_pairs[[y]],],by=c('Longitude','Latitude'))
    # bind everything together     
    all_sites <- bind_rows(site_x1,site_y1)
    # combine time series into columns from which easier to calculate pairwise correlations 
    all_sites_wide <- as.data.frame(all_sites %>% dcast(Year ~ Alpha_Number + Site_number + Variable_number + Unit, value.var = "Value"))
    # keep just those time series exceeding minimum length defined by min_num_yrs
    all_sites_wide <- subset(all_sites_wide,select=colSums(is.na(all_sites_wide)==F) >= min_num_yrs)
    if(ncol(all_sites_wide)>2){
      # calculate spatial synchrony with all other sites/indices      
      cor_mat_sites <- suppressWarnings(cor(subset(all_sites_wide,select=-c(Year)),method='spearman',use='pairwise.complete.obs'))
      # h/t https://stackoverflow.com/questions/43283103/speeding-up-count-of-pairwise-observations-in-r
      # limit pairwise comparisons to those with min_num_yrs 
      NAmat = matrix(0, nrow = nrow(subset(all_sites_wide,select=-c(Year))), ncol = ncol(subset(all_sites_wide,select=-c(Year))))
      NAmat[ !is.na(subset(all_sites_wide,select=-c(Year))) ] = 1
      filter2 = (tcrossprod(t(NAmat)) < min_num_yrs)
      cor_mat_sites[filter2] <- NA
      # keep just the lower triangle
      cor_mat_sites <- cor_mat_sites[lower.tri(cor_mat_sites)]
      return(cor_mat_sites[is.na(cor_mat_sites)==F])
    }
  })) 
  if(is.null(all_cors)){ output <- cbind(NA,NA,0) } else(    
    output <- cbind(mean(all_cors),sd(all_cors),length(all_cors)) )
  colnames(output) <- c('mean','sd','N')
  rownames(output) <- species
  return(output)
}  

## Load mastree Data ----
load("mastree.RData") # raw mastree data
#mastree+ data can be accessed through the MASTREE+ app: https://mastreeplus.shinyapps.io/mastreeplus/

# Data subsetting ----
#loop through different time series duration filters
ts.lengths <- 3:6
for(i in 1:length(ts.lengths)){
sample <- mastree %>% filter(VarType == "C", 
                            !(Variable %in% c("flower", "pollen", "pollen_sediment", "dendro")),# These variables are excluded
                            !grepl("spp", Species),
                            !(Species %in% c("Mixed species", "MIXED", "MIXED SPECIES"))) %>% # exclude records where the exact species isn't recorded, but noted instead as spp.
  # Remove series that have fewer than the ith threshold of continuous records
  filter(Length >= ts.lengths[i]) # remove series with shorter duration than the threshold.

## Detect multiple variables within a paper (= duplication within dataset of same raw data)====                        
sample  %>% 
  group_by(Alpha_Number,Species, Segment, Site_number, Latitude, Longitude) %>%
  mutate(mult = n_distinct(Variable), multunit = n_distinct(Unit)) %>%
  filter(mult >1 | multunit >1) %>% select(-Year, -Value) %>% unique() 

## remove case with >1 variable per species+location+sampling effort. Note: Paper duplicates removed in other script ====
dataset <- sample %>% filter(!(Alpha_Number == 2505 & Variable_number == 1),
                             !(Alpha_Number == 1500 & Variable_number == 2),
                             !(Alpha_Number == 1472 & Variable_number !=  1))

pooled_dataset <- dataset %>% mutate(Unit = case_when(Unit == "fruits/individual" ~ "seeds/individual", 
                                                      TRUE ~ Unit))

# Basic indices ----

basic_indices <- dataset %>%
  group_by(Alpha_Number, Segment, Site_number, Variable_number, Species) %>% 
  mutate( 
    time_series = cur_group_id(),
    # basic metrics Standard deviation, Mean, CV
    SD = sd(Value, na.rm = TRUE), 
    M = mean(Value, na.rm = TRUE), 
    CV = SD/M) %>% 
  ungroup()

# Add AR indices ----
ar_indices <- basic_indices %>% 
  # put all data in a time series friendly format 
  select(time_series, Year, Value) %>% 
  arrange(time_series, Year) %>%   
  group_by(time_series) %>% 
  # ACF correlation at lag 1,2,3
  mutate(
    AR1 = tibble(Value) %>% any_acf(., lag = 1) # Autocorrelation at lag 1. Also: tibble ensures sig_acf is applied to each group separately
  ) %>% 
  # rejoin indices columns and select the vars of interest
  left_join(basic_indices, by = c("time_series", "Year", "Value")) 

# Spatial indices ----

# initial calculation
spat_10km <- t(sapply(unique(dataset$Species), function(z){spat_syn(data_tibble=dataset,species= z, scale_km=10,global=F,min_num_yrs=3) }))
spat_100km <- t(sapply(unique(dataset$Species), function(z){spat_syn(data_tibble=dataset,species= z, scale_km=100,global=F,min_num_yrs=3) }))
spat_all <- t(sapply(unique(dataset$Species), function(z){spat_syn(data_tibble=dataset,species=z,global=T,min_num_yrs=3) }))

# reformat and combine
spat_10km <- data.frame(Species = row.names(spat_10km), spat_10km) %>% as_tibble() %>% mutate(km = "10") %>% select(Species, mean_spat = X1, sd_spat = X2, N_spat = X3, km)
spat_100km <- data.frame(Species = row.names(spat_100km), spat_100km) %>% as_tibble() %>% mutate(km = "100") %>% select(Species, mean_spat = X1, sd_spat = X2, N_spat = X3, km)
spat_all <- data.frame(Species = row.names(spat_all), spat_all) %>% as_tibble() %>% mutate(km = "global") %>% select(Species, mean_spat = X1, sd_spat = X2, N_spat = X3, km)

combined_spat <- bind_rows(spat_10km, spat_100km, spat_all)

# summary table with spatial information added ----

# summarize the tibble and add spatial info
sum_all_indices <- ar_indices %>% ungroup() %>% 
  group_by(Species, Variable, Unit) %>% 
  summarise(
    mean_AR1 = mean(AR1),
    sd_AR1 = sd(AR1),
    mean_CV = mean(CV),
    sd_CV = sd(CV),
    n_timeseries = n_distinct(time_series),
    mean_duration=mean(Length),
    total_duration=sum(Length),
    mean_latitude=mean(Latitude),
    range_latitude=(max(Latitude)-min(Latitude)),
    mean_latitude_absolute= mean(abs(Latitude))
    ) %>%  
  #join the spatial info
  left_join(combined_spat, by = "Species") %>% 
  pivot_wider(names_from = km, values_from = c(mean_spat, sd_spat, N_spat))

filename.i <- paste("Summary_All_indices_", ts.lengths[i], "_years.Rdata", sep="")
save(sum_all_indices, file =filename.i)

# repeat for pooled object ----

# Basic indices ----

basic_indices_pooled <- pooled_dataset %>%
  group_by(Alpha_Number, Segment, Site_number, Variable_number, Species) %>% 
  mutate( 
    time_series = cur_group_id(),
    # basic metrics Standard deviation, Mean, CV
    SD = sd(Value, na.rm = TRUE), 
    M = mean(Value, na.rm = TRUE), 
    CV = SD/M) %>% 
  ungroup()

# Add AR indices ----
ar_indices_pooled <- basic_indices_pooled %>% 
  # put all data in a time series friendly format 
  select(time_series, Year, Value) %>% 
  arrange(time_series, Year) %>%   
  group_by(time_series) %>% 
  # ACF correlation at lag 1,2,3
  mutate(
    AR1 = tibble(Value) %>% any_acf(., lag = 1) # Autocorrelation at lag 1. Also: tibble ensures sig_acf is applied to each group separately
  ) %>% 
  # rejoin indices columns and select the vars of interest
  left_join(basic_indices_pooled, by = c("time_series", "Year", "Value")) 

# Spatial indices ----

# initial calculation
spat_10km_p <- t(sapply(unique(pooled_dataset$Species), function(z){spat_syn(data_tibble=pooled_dataset,species= z, scale_km=10,global=F,min_num_yrs=3) }))
spat_100km_p <- t(sapply(unique(pooled_dataset$Species), function(z){spat_syn(data_tibble=pooled_dataset,species= z, scale_km=100,global=F,min_num_yrs=3) }))
spat_all_p <- t(sapply(unique(pooled_dataset$Species), function(z){spat_syn(data_tibble=pooled_dataset,species=z,global=T,min_num_yrs=3) }))

# reformat and combine
spat_10km_p <- data.frame(Species = row.names(spat_10km_p), spat_10km_p) %>% as_tibble() %>% mutate(km = "10") %>% select(Species, mean_spat = X1, sd_spat = X2, N_spat = X3, km)
spat_100km_p <- data.frame(Species = row.names(spat_100km_p), spat_100km_p) %>% as_tibble() %>% mutate(km = "100") %>% select(Species, mean_spat = X1, sd_spat = X2, N_spat = X3, km)
spat_all_p <- data.frame(Species = row.names(spat_all_p), spat_all_p) %>% as_tibble() %>% mutate(km = "global") %>% select(Species, mean_spat = X1, sd_spat = X2, N_spat = X3, km)

combined_spat_p <- bind_rows(spat_10km_p, spat_100km_p, spat_all_p)

# summary table with spatial information added ----

# summarize the tibble and add spatial info
sum_all_indices_pooled <- ar_indices_pooled %>% ungroup() %>% 
  group_by(Species, Variable, Unit) %>% 
  summarise(
    mean_AR1 = mean(AR1),
    sd_AR1 = sd(AR1),
    mean_CV = mean(CV),
    sd_CV = sd(CV),
    n_timeseries = n_distinct(time_series),
    mean_duration=mean(Length),
    total_duration=sum(Length),
    mean_latitude=mean(Latitude),
    range_latitude=(max(Latitude)-min(Latitude)),
    mean_latitude_absolute= mean(abs(Latitude))
  ) %>%  
  #join the spatial info
  left_join(combined_spat_p, by = "Species") %>% 
  pivot_wider(names_from = km, values_from = c(mean_spat, sd_spat, N_spat))

save(sum_all_indices_pooled, file =gsub(".Rdata", "_pooled.Rdata", filename.i))

  }
