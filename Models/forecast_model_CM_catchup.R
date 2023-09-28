
# ---
# Camille Minaudo
# title: "NEON forecast challenge - GLEON2022"
# date: "Nov 2022"
# https://github.com/camilleminaudo/neon4cast
# ---


# cat("\014")
# rm(list = ls())

library(tidyverse)
library(neon4cast)
library(lubridate)
library(rMR)

scriptpath <- dirname(rstudioapi::getSourceEditorContext()$path) # path of current script file
path_repo_root <- dirname(scriptpath)
source(paste(path_repo_root,'generate_GLEON_lm_lag_1day.R', sep = "/"))
source(paste(path_repo_root,'download_noaa.R', sep = "/"))

target_url <- "https://data.ecoforecast.org/neon4cast-targets/aquatics/aquatics-targets.csv.gz"

sites <- read_csv("https://raw.githubusercontent.com/eco4cast/neon4cast-targets/main/NEON_Field_Site_Metadata_20220412.csv") |> 
  dplyr::filter(aquatics == 1) |> 
  dplyr::select(field_site_id) |> 
  dplyr::pull()

# generate ensembles
if (dir.exists('./Forecasts') != T) {
  dir.create('./Forecasts', recursive = T)
}




# check for any missing forecasts
message("==== Checking for missed GLEON_lm_lag_1day forecasts ====")
challenge_model_name <- 'GLEON_lm_lag_1day'


# Dates of forecasts 
today <- paste(Sys.Date() - days(2), '00:00:00')
this_year <- data.frame(date = as.character(paste0(seq.Date(as_date('2023-01-01'), to = as_date(today), by = 'day'), ' 00:00:00')),
                        exists = NA)

# what forecasts have already been submitted?
challenge_s3_region <- "data"
challenge_s3_endpoint <- "ecoforecast.org"

# is that file present in the bucket?
for (i in 1:nrow(this_year)) {
  forecast_file <- paste0('aquatics-', as_date(this_year$date[i]), '-', challenge_model_name, '.csv.gz')
  
  this_year$exists[i] <- suppressMessages(aws.s3::object_exists(object = file.path("raw", 'aquatics', forecast_file),
                                                                bucket = "neon4cast-forecasts",
                                                                region = challenge_s3_region,
                                                                base_url = challenge_s3_endpoint))
}

# which dates do you need to generate forecasts for?
missed_dates <- this_year |> 
  filter(exists == F) |> 
  pull(date) |> 
  as_date()

if (length(missed_dates) != 0) {
  message('there are new forecasts to be submitted!')
  
  for (i in 1:length(missed_dates)) {
    message(paste("Running forecast",i,"out of",length(missed_dates),sep = " "))
    forecast_date <- missed_dates[i]
    
    # Download the NOAA data needed
    noaa_data <- download_noaa(sites = sites,
                               forecast_date = forecast_date)
    
    noaa_future = noaa_data$future
    noaa_past = noaa_data$past
    
    
    forecast_file <- generate_GLEON_lm_lag_1day(forecast_date = forecast_date,
                                                forecast_name = 'GLEON_lm_lag_1day', 
                                                target_url = target_url,
                                                sites = sites,
                                                noaa_future = noaa_future,
                                                noaa_past = noaa_past)
    
    
    if (!is.na(forecast_file)) {
      neon4cast::submit(file.path('./Forecasts/', forecast_file), ask = F)
      message('submitting new GLEON_lm_lag_1day forecast')
    }
    
  }
} else {
  message('no new forecasts')  
}






