installed_packages <- installed.packages()
if (is.element('rMR', installed_packages) == F) {
  install.packages('rMR')
}

# my scripts
source("./Models/forecast_model_CM.R")
message('My dynamic model successfully submitted!')



