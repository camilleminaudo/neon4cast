


Sys.unsetenv("AWS_ACCESS_KEY_ID")
Sys.unsetenv("AWS_SECRET_ACCESS_KEY")
Sys.unsetenv("AWS_DEFAULT_REGION")
Sys.unsetenv("AWS_S3_ENDPOINT")
Sys.setenv(AWS_EC2_METADATA_DISABLED="TRUE")
library(arrow)
library(tidyverse)

s3 <- s3_bucket(paste0("neon4cast-scores/parquet/aquatics"), endpoint_override= "data.ecoforecast.org")

your_model_id <- 'GLEON_lm_lag_1day'

scores <- open_dataset(s3) |>
  filter(model_id == your_model_id) |> 
  collect()

my_scores <- scores |> 
  group_by(model_id, reference_datetime, site_id, datetime) |> 
  summarize(crps_score = mean(crps, na.rm = TRUE)) 


ggplot(my_scores, aes(datetime, crps_score, colour = site_id ))+geom_point()+theme_bw()+facet_wrap(site_id~.)
