



generate_GLEON_lm_lag_1day <- function(forecast_date, 
                                       forecast_name = 'GLEON_lm_lag_1day',
                                       target_url,
                                       sites,
                                       noaa_future, #dateframe containing the relevant noaa data
                                       noaa_past) { #dateframe containing the relevant noaa data
  
  # First check that the NOAA data looks right
  if (min(noaa_future$datetime) != forecast_date) {
    return(NA)
    stop('Error: NOAA forecast is not correct!')
  }
  
  targets <- readr::read_csv(target_url, guess_max = 10000)
  
  target <- targets |> 
    select(datetime, site_id, variable, observation) |> 
    filter(variable == 'temperature' | variable == 'oxygen')
  
  
  forecast <- NULL
  
  
  site_data <- readr::read_csv("https://raw.githubusercontent.com/eco4cast/neon4cast-targets/main/NEON_Field_Site_Metadata_20220412.csv") |>
    dplyr::filter(aquatics == 1)
  site_data <- site_data %>%
    dplyr::filter(field_site_subtype == 'Lake')
  
  
  
  
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
  get_lag <- function(x,mylag) {
    out<-lagmatrix(x,mylag+1)[,mylag]
    return(out[1:length(out)-1])
  }
  
  
  # ------------ actual forecast ---------------
  
  for(i in 1:length(sites)){
    
    message(paste0("Running site: ", sites[i]))
    
    # Get site information for elevation
    site_info <- site_data %>% dplyr::filter(field_site_id == sites[i])
    
    print("loading past data")
    my_noaa_past <- noaa_past |>
      dplyr::filter(site_id == sites[i]) 
    
    
    my_noaa_future <- noaa_future |>
      dplyr::filter(site_id == sites[i]) 
    
    
    
    #Merge in past NOAA data into the targets file, matching by date.
    print("Merge in past NOAA data into the targets file")
    site_target <- target |>
      select(datetime, site_id, variable, observation) |>
      dplyr::filter(variable %in% c("temperature", "oxygen"),
                    site_id == sites[i]) |>
      pivot_wider(names_from = "variable", values_from = "observation") |>
      left_join(my_noaa_past, by = c("datetime"))
    
    print("...done!")
    
    #Check that temperature and oxygen are available at site
    if("temperature" %in% names(site_target) & "oxygen" %in% names(site_target)){
      
      if(length(which(!is.na(site_target$air_temperature) & !is.na(site_target$temperature))) > 0){
        
        print("Build model")
        
        #Fit linear model based on past data: water temperature = m * air temperature + b
        # fit <- lm(site_target$temperature~site_target$air_temperature)
        
        # linear fit with lag of 1 day : using today's obs to predict tomorrow
        site_target$lag <- I(get_lag(site_target$temperature,1))
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
        if(dim(site_info)[1]>0){
          print("Predict oxygen")
          forecasted_oxygen <- rMR::Eq.Ox.conc(forecasted_temperature, elevation.m = site_info$field_mean_elevation_m,
                                               bar.press = NULL,
                                               bar.units = NULL,
                                               out.DO.meas = "mg/L",
                                               salinity = 0,
                                               salinity.units = "pp.thou")
        } else {
          print("Cannot predict oxygen because elevation data is missing from site_info file")
          forecasted_oxygen <- NA*forecasted_temperature
        }
        
        
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
  
  my_forecast_EFI <- forecast %>%
    mutate(model_id = forecast_name,
           reference_datetime = as_date(min(datetime)) - days(2),
           family = 'ensemble',
           parameter = as.character(variable)) %>%
    select(model_id, datetime, reference_datetime, site_id, family, parameter, variable, prediction)
  
  
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
  
}
















generate_climatology <- function(forecast_date, 
                                 forecast_name = 'climatology',
                                 target_url,
                                 sites) {
  
  
  targets <- readr::read_csv(target_url, guess_max = 10000)
  
  # calculates a doy average for each target variable in each site
  target_clim <- targets %>%  
    mutate(doy = yday(datetime),
           year = year(datetime)) %>% 
    filter(year < year(forecast_date),
           site_id %in% sites) |> 
    group_by(doy, site_id, variable) %>% 
    summarise(mean = mean(observation, na.rm = TRUE),
              sd = sd(observation, na.rm = TRUE),
              .groups = "drop") %>% 
    mutate(mean = ifelse(is.nan(mean), NA, mean)) 
  
  # subset to the DOY needed
  curr_month <- month(forecast_date)
  if(curr_month < 10){
    curr_month <- paste0("0", curr_month)
  }
  
  curr_year <- year(forecast_date)
  start_date <- forecast_date + days(1)
  
  forecast_dates <- seq(start_date, as_date(start_date + days(34)), "1 day")
  forecast_doy <- yday(forecast_dates)
  
  forecast_dates_df <- tibble(datetime = forecast_dates,
                              doy = forecast_doy)
  # generate forecast
  forecast <- target_clim %>%
    mutate(doy = as.integer(doy)) %>% 
    filter(doy %in% forecast_doy) %>% 
    full_join(forecast_dates_df, by = 'doy') %>%
    arrange(site_id, datetime)
  
  
  subseted_site_names <- unique(forecast$site_id)
  site_vector <- NULL
  for(i in 1:length(subseted_site_names)){
    site_vector <- c(site_vector, rep(subseted_site_names[i], length(forecast_dates)))
  }
  
  forecast_tibble1 <- tibble(datetime = rep(forecast_dates, length(subseted_site_names)),
                             site_id = site_vector,
                             variable = "temperature")
  
  forecast_tibble2 <- tibble(datetime = rep(forecast_dates, length(subseted_site_names)),
                             site_id = site_vector,
                             variable = "oxygen")
  
  forecast_tibble3 <- tibble(datetime = rep(forecast_dates, length(subseted_site_names)),
                             site_id = site_vector,
                             variable = "chla")
  
  forecast_tibble <- bind_rows(forecast_tibble1, forecast_tibble2, forecast_tibble3)
  
  foreast <- right_join(forecast, forecast_tibble)
  
  combined <- forecast %>% 
    select(datetime, site_id, variable, mean, sd) %>% 
    group_by(site_id, variable) %>% 
    # remove rows where all in group are NA
    filter(all(!is.na(mean))) %>%
    # retain rows where group size >= 2, to allow interpolation
    filter(n() >= 2) %>%
    mutate(mu = imputeTS::na_interpolation(mean),
           sigma = median(sd, na.rm = TRUE)) %>%
    pivot_longer(c("mu", "sigma"),names_to = "parameter", values_to = "prediction") |> 
    mutate(family = "normal") |> 
    ungroup() |> 
    mutate(reference_datetime = lubridate::as_date(min(datetime)) - lubridate::days(1),
           model_id = forecast_name) |> 
    select(model_id, datetime, reference_datetime, site_id, family, parameter, variable, prediction)
  
  # write forecast
  file_date <- combined$reference_datetime[1]
  
  file_name <- paste0('aquatics-', forecast_date, '-', forecast_name, '.csv.gz')
  
  readr::write_csv(combined, file.path('Forecasts', file_name))  
  
  message(forecast_name, ' generated for ', forecast_date)
  
  
  valid <- neon4cast::forecast_output_validator(file.path('./Forecasts', file_name))
  
  
  if (!valid) {
    file.remove(file.path('./Forecasts/', file_name))
    message('forecast not valid')
  } else {
    return(file_name)
  }
  
}



generate_persistenceRW <- function(forecast_date, 
                                   forecast_name = 'persistenceRW',
                                   target_url,
                                   sites) {
  source('R/fablePersistenceRW.R')
  # 1.Read in the targets data
  targets <- readr::read_csv(target_url, guess_max = 10000) |> 
    mutate(observation = ifelse(observation == 0 & variable == "chla", 0.00001, observation))
  
  
  # 2. Make the targets into a tsibble with explicit gaps
  targets_ts <- targets %>%
    filter(datetime < forecast_date) |> 
    as_tsibble(key = c('variable', 'site_id'), index = 'datetime') %>%
    # add NA values up to today (index)
    tsibble::fill_gaps(.end = forecast_date)
  
  
  # 3. Run through each via map
  # Requires a dataframe that has each of the variable in the RW_forecast function
  site_var_combinations <- expand.grid(site = unique(targets$site_id),
                                       var = unique(targets$variable)) %>%
    # assign the transformation depending on the variable. chla and oxygen get a log(x) transformation
    mutate(transformation = ifelse(var %in% c('chla', 'oxygen'), 'log', 'none')) %>%
    mutate(boot_number = 200,
           ref_date = forecast_date,
           h = 35,
           bootstrap = T, 
           verbose = T)
  
  # runs the RW forecast for each combination of variable and site_id
  RW_forecasts <- purrr::pmap_dfr(site_var_combinations, RW_daily_forecast) 
  # convert the output into EFI standard
  RW_forecasts_EFI <- RW_forecasts %>%
    rename(parameter = .rep,
           prediction = .sim) %>%
    # For the EFI challenge we only want the forecast for future
    filter(datetime > forecast_date) %>%
    group_by(site_id, variable) %>%
    mutate(reference_datetime = min(datetime) - lubridate::days(1),
           family = "ensemble",
           model_id = forecast_name) %>%
    select(model_id, datetime, reference_datetime, site_id, family, parameter, variable, prediction) 
  
  # write forecast
  file_date <- RW_forecasts_EFI$reference_datetime[1]
  
  file_name <- paste0('aquatics-', file_date, '-', forecast_name, '.csv.gz')
  
  readr::write_csv(RW_forecasts_EFI, file.path('Forecasts', file_name))  
  
  message(forecast_name, ' generated for ', forecast_date)
  
  
  valid <- neon4cast::forecast_output_validator(file.path('./Forecasts', file_name))
  
  
  if (!valid) {
    file.remove(file.path('./Forecasts/', file_name))
    message('forecast not valid')
  } else {
    return(file_name)
  }
  
}