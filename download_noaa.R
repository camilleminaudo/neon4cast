# Functionise the download noaa script
# Returns a list of two dataframes, one containing past data and one containing future data
download_noaa <- function(sites,
                          forecast_date) {
  
  
  # Past stacked weather
  df_past <- neon4cast::noaa_stage3()
  # Forecasts
  
  # New forecast only available at 5am UTC the next day
  noaa_date <- forecast_date - days(1)

  if (class(try(neon4cast::noaa_stage2(start_date = noaa_date), silent = T))[1] == 'try-error') {
    return(NA)
    stop('No NOAA forecast available for ', noaa_date, ' skipping')
    
  } else {
    df_future <- neon4cast::noaa_stage2(start_date = noaa_date)
  }
  
  # Only need the air temperature from the test_site
  noaa_past <- df_past |> 
    dplyr::filter(site_id %in% sites,
                  datetime >= ymd('2017-01-01'),
                  datetime < forecast_date,
                  variable == "air_temperature") |>  
    dplyr::collect()
  message('stage 3 data downloaded')
  
  
  # aggregate the past to mean daily values
  noaa_past_agg <- noaa_past |> 
    mutate(datetime = as_date(datetime)) |> 
    group_by(datetime, site_id) |> 
    summarize(air_temperature = mean(prediction, na.rm = TRUE), .groups = "drop") |> 
    rename(datetime = datetime) |> 
    # convert air temp to C
    mutate(air_temperature = air_temperature - 273.15)
  
  
  noaa_future <- df_future |> 
    dplyr::filter(datetime >= forecast_date,
                  site_id %in% sites,
                  variable == "air_temperature") |> 
    dplyr::collect()
  
  # Aggregate for each ensemble for future
  noaa_future_agg <- noaa_future |> 
    mutate(datetime = as_date(datetime)) |> 
    group_by(datetime, site_id, parameter) |> 
    summarize(air_temperature = mean(prediction)) |> 
    mutate(air_temperature = air_temperature - 273.15) |> 
    select(datetime, site_id, air_temperature, parameter)
  
  message('stage 2 data downloaded')
  
 return(list(future = noaa_future_agg, 
             past = noaa_past_agg)) 
}
