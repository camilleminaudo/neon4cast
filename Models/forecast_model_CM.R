
# ---
  # Camille Minaudo
  # title: "NEON forecast challenge - GLEON2022"
  # date: "Nov 2022"
# ---
  

cat("\014")
rm(list = ls())

library(tidyverse)
library(neon4cast)
library(lubridate)
library(rMR)

forecast_date <- Sys.Date()
noaa_date <- Sys.Date() - days(2)  #Need to use yesterday's NOAA forecast because today's is not available yet

#Step 0: Define team name and team members 

model_id <- "GLEON_lm_lag_1day"

#Step 1: Download latest target data and site description data

target <- readr::read_csv("https://data.ecoforecast.org/neon4cast-targets/aquatics/aquatics-targets.csv.gz", guess_max = 1e6)

site_data <- readr::read_csv("https://raw.githubusercontent.com/eco4cast/neon4cast-targets/main/NEON_Field_Site_Metadata_20220412.csv") |> 
  dplyr::filter(aquatics == 1)
site_data <- site_data %>%
  filter(field_site_subtype == 'Lake')

target <- target %>%
  dplyr::filter(site_id %in% site_data$field_site_id)

# selecting only water temperature and oxygen as target variables to predict
target <- target %>%
  dplyr::filter(variable == 'temperature' | variable == 'oxygen')


#Step 2: Get drivers

df_past <- neon4cast::noaa_stage3()

df_future <- neon4cast::noaa_stage2()

sites <- unique(target$site_id)

#Step 3.0: Generate forecasts for each site

forecast <- NULL

# sites <- sites[1:5]


# ------------ Some needed functions ---------------

ipredict <-function(model,newdata) {
  
  i=1
  newdata$temperature[i] <- predict(model,newdata=newdata[i,])[i]
  
  for (i in seq(2,dim(newdata)[1])) {
    if (is.na(newdata$temperature[i])) {
      newdata$lag[i] <- newdata$temperature[i-1]
      newdata$temperature[i] <-predict(model,newdata=newdata[seq(1,i),])[i]
    }
  }
  P = as.numeric(unlist(newdata$temperature))
  return(P)
}

#a function to calculate lags
lagmatrix <- function(x,max.lag){embed(c(rep(NA,max.lag),x),max.lag)}
lag <- function(x,lag) {
  out<-lagmatrix(x,lag+1)[,lag]
  return(out[1:length(out)-1])
}


# ------------ actual forecast ---------------

 for(i in 1:length(sites)){

  message(paste0("Running site: ", sites[i]))

  # Get site information for elevation
  site_info <- site_data %>% dplyr::filter(field_site_id == sites[i])

  print("loading past data")
  noaa_past <- df_past |>
    dplyr::filter(site_id == sites[i],
                  variable == "air_temperature") |>
    dplyr::rename(ensemble = parameter) %>%
    dplyr::select(datetime, prediction, ensemble) |>
    dplyr::collect()
  
  print("...done!")

  
  print("loading future data")
  
  variables <- c("air_temperature")
  
  noaa_future <- df_future |> 
    dplyr::filter(reference_datetime == noaa_date,
                  datetime >= forecast_date,
                  site_id %in% sites[i],
                  variable %in% variables) |> 
    dplyr::collect()
  print("...done!")

  # Aggregate (to day) and convert units of drivers
  print("aggregate to daily")
  
  noaa_past_mean <- noaa_past %>%
    mutate(date = as_date(datetime)) %>%
    group_by(date) %>%
    summarize(air_temperature = mean(prediction, na.rm = TRUE), .groups = "drop") %>%
    rename(datetime = date) %>%
    mutate(air_temperature = air_temperature - 273.15)

  
  noaa_future_daily <- noaa_future |> 
    mutate(datetime = as_date(datetime)) |> 
    # mean daily forecasts at each site per ensemble
    group_by(datetime, site_id, parameter, variable) |> 
    summarize(prediction = mean(prediction)) |>
    pivot_wider(names_from = variable, values_from = prediction) |>
    # convert to Celsius
    mutate(air_temperature = air_temperature - 273.15) |> 
    select(datetime, site_id, air_temperature, parameter)
  print("...done!")
  
  noaa_future_site <- noaa_future_daily |> 
    dplyr::filter(site_id == sites[i])
  
  # noaa_future_site <- noaa_future %>%
  #   mutate(datetime = as_date(datetime)) %>%
  #   group_by(datetime, ensemble) %>%
  #   summarize(air_temperature = mean(prediction), .groups = "drop") |>
  #   mutate(air_temperature = air_temperature - 273.15) |>
  #   select(datetime, air_temperature, ensemble)

  #Merge in past NOAA data into the targets file, matching by date.
  print("Merge in past NOAA data into the targets file")
  site_target <- target |>
    select(datetime, site_id, variable, observation) |>
    dplyr::filter(variable %in% c("temperature", "oxygen"),
           site_id == sites[i]) |>
    pivot_wider(names_from = "variable", values_from = "observation") |>
    left_join(noaa_past_mean, by = c("datetime"))
  
  print("...done!")
  
  #Check that temperature and oxygen are available at site
  if("temperature" %in% names(site_target) & "oxygen" %in% names(site_target)){

    if(length(which(!is.na(site_target$air_temperature) & !is.na(site_target$temperature))) > 0){
      
      print("Build model")
      
      #Fit linear model based on past data: water temperature = m * air temperature + b
      # fit <- lm(site_target$temperature~site_target$air_temperature)
      
      # linear fit with lag of 1 day : using today's obs to predict tomorrow
      site_target$lag <- I(lag(site_target$temperature,1))
      fit <- lm(data = site_target, temperature ~ air_temperature+lag)
      

      #use linear regression to forecast water temperature for each ensemble member
      print("Predict temperature")
      # forecasted_temperature <- fit$coefficients[1] + fit$coefficients[2] * noaa_future_site$air_temperature
      
      # prediction: using today's obs to predict tomorrow
      noaa_future_site$temperature <- noaa_future_site$lag <- NA * noaa_future_site$air_temperature
      
      # first, using model to predict temperature until supposed day of observation
      temperature_extended <- ipredict(model = fit, newdata = site_target)
      site_target$temperature[is.na(site_target$temperature)] <- temperature_extended[is.na(site_target$temperature)]
      
      noaa_future_site$lag[1] <- tail(site_target$temperature[!is.na(site_target$temperature)], 1)
      
      forecasted_temperature <- ipredict(model = fit, newdata = noaa_future_site)
      

      #use forecasted temperature to predict oyxgen by assuming that oxygen is saturated.
      print("Predict oxygen")
      forecasted_oxygen <- rMR::Eq.Ox.conc(forecasted_temperature, elevation.m = ,site_info$field_mean_elevation_m,
                                           bar.press = NULL,
                                           bar.units = NULL,
                                           out.DO.meas = "mg/L",
                                           salinity = 0,
                                           salinity.units = "pp.thou")

      temperature <- tibble(datetime = noaa_future_site$datetime,
                            site_id = sites[i],
                            ensemble = noaa_future_site$ensemble,
                            prediction = forecasted_temperature,
                            variable = "temperature")

      oxygen <- tibble(datetime = noaa_future_site$datetime,
                       site_id = sites[i],
                       ensemble = noaa_future_site$ensemble,
                       prediction = forecasted_oxygen,
                       variable = "oxygen")


      #Build site level dataframe.  Note we are not forecasting chla
      forecast <- dplyr::bind_rows(forecast, temperature, oxygen)
      print("<=== done for this site ===>")
    }
  }
}


# ------------ save and submit  ---------------


my_forecast_EFI <- forecast %>%
  mutate(model_id = model_id,
         reference_datetime = as_date(min(datetime)) - days(2),
         family = 'ensemble',
         parameter = as.character(variable)) %>%
  select(model_id, datetime, reference_datetime, site_id, family, parameter, variable, prediction)


# forecast <- forecast |>
#   mutate(reference_datetime = forecast_date,
#          family = "ensemble",
#          model_id = model_id) |>
#   rename(parameter = ensemble) |>
#   select(model_id, datetime, reference_datetime, site_id, family, parameter, variable, prediction)

#Visualize forecast.  Is it reasonable?

ggplot(my_forecast_EFI, aes(datetime, prediction, colour = site_id))+geom_path()+theme_bw()+facet_wrap(variable~., scales = "free_y")

#Forecast output file name in standards requires for Challenge.
# csv.gz means that it will be compressed
file_date <- my_forecast_EFI$reference_datetime[1]
forecast_file <- paste0("aquatics","-",file_date,"-",model_id,".csv.gz")

#Write csv to disk
# setwd("~/my_forecasts")

write_csv(my_forecast_EFI, forecast_file)

neon4cast::forecast_output_validator(forecast_file)

# Step 4: Submit forecast!

neon4cast::submit(forecast_file = forecast_file, metadata = NULL, ask = FALSE)
