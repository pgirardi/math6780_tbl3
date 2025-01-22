##' Download Targets
##' @return data.frame in long format with days as rows, and time, site_id, variable, and observed as columns
download_targets <- function(){
  readr::read_csv("https://data.ecoforecast.org/neon4cast-targets/aquatics/aquatics-targets.csv.gz", guess_max = 1e6)
}

##' Download Site metadata
##' @return metadata dataframe
download_site_meta <- function(){
  site_data <- readr::read_csv("https://raw.githubusercontent.com/eco4cast/neon4cast-targets/main/NEON_Field_Site_Metadata_20220412.csv") 
  site_data %>% filter(as.integer(aquatics) == 1)
}


##' append historical meteorological data into target file
##' @param target targets dataframe
##' @return updated targets dataframe with added weather data
merge_met_past <- function(target){
  
  ## connect to data
  df_past <- neon4cast::noaa_stage3()
  
  ## filter for site and variable
  sites <- unique(target$site_id)
  
  ## temporary hack to remove a site that's mid-behaving
  sites = sites[!(sites=="POSE")] 
  target = target |> filter(site_id %in% sites)  
  
  ## grab air temperature from the historical forecast
  noaa_past <- df_past |> 
    dplyr::filter(site_id %in% sites,
                  variable == "air_temperature") |> 
    dplyr::collect()
  
  ## aggregate to daily
  noaa_past_mean = noaa_past |> 
    mutate(datetime = as.Date(datetime)) |>
    group_by(datetime, site_id) |> 
    summarise(air_temperature = mean(prediction),.groups = "drop")
  
  ## Aggregate (to day) and convert units of drivers
  target <- target %>% 
    group_by(datetime, site_id,variable) %>%
    summarize(obs2 = mean(observation, na.rm = TRUE), .groups = "drop") %>%
    mutate(obs3 = ifelse(is.nan(obs2),NA,obs2)) %>%
    select(datetime, site_id, variable, obs3) %>%
    rename(observation = obs3) %>%
    filter(variable %in% c("temperature", "oxygen")) %>% 
    tidyr::pivot_wider(names_from = "variable", values_from = "observation")
  
  ## Merge in past NOAA data into the targets file, matching by date.
  target <- left_join(target, noaa_past_mean, by = c("datetime","site_id"))
  
}

##' Download NOAA GEFS weather forecast
##' @param forecast_date start date of forecast
##' @return dataframe
download_met_forecast <- function(forecast_date){
  noaa_date <- forecast_date - lubridate::days(1)  #Need to use yesterday's NOAA forecast because today's is not available yet
  
  ## connect to data
  df_future <- neon4cast::noaa_stage2(start_date = as.character(noaa_date))
  
  ## filter available forecasts by date and variable
  met_future <- df_future |> 
    dplyr::filter(datetime >= lubridate::as_datetime(forecast_date), 
                  variable == "air_temperature") |> 
    dplyr::collect()
  
  ## aggregate to daily
  met_future <- met_future %>% 
    mutate(datetime = lubridate::as_date(datetime)) %>% 
    group_by(datetime, site_id, parameter) |> 
    summarize(air_temperature = mean(prediction), .groups = "drop") |> 
    #    mutate(air_temperature = air_temperature - 273.15) |> 
    select(datetime, site_id, air_temperature, parameter)
  
  return(met_future)
}
```

### OWNER

1. In RStudio, click File > New File > R Script
2. Copy and Paste the above function into this file
3. Save the file as "01_download_data.R"
4. From the Git tab, click the box next to the file you just created. This is equivalent to _git add_
5. Click Commit, enter a log message, and click Commit. This is equivalent to _git commit_
6. To push the change up to Github click on the green up arrow. This is equivalent to _git push_

## Task 3: Collaborator adds model calibration and forecast

With the first function complete, let's now imagine that a **COLLABORATOR** has been tasked with adding the next two functions. To do so they must first fork and clone the repository

### COLLABORATOR

1. Go to Github and navigate to the project repository within the OWNER's workspace.
2. Click Fork, which will make a copy of the repository to your own workspace.
3. Copy the URL to your own version and follow the instructions above for cloning the repository in RStudio.
4. Open a new file, enter the code below, and then save the file as "02_calibrate_forecast.R"
```{r}
##' Calibrate aquatic forecast model
calibrate_forecast <- function(target){
  fit <- list() 
  sites <- unique(target$site_id)
  for(i in 1:length(sites)){
    site_target <- target |> 
      filter(site_id == sites[i])
    
    if(length(which(!is.na(site_target$air_temperature) & !is.na(site_target$temperature))) > 0){
      #Fit linear model based on past data: water temperature = m * air temperature + b
      fit[[i]] <- lm(temperature~air_temperature,data = site_target)
    }
  }
  names(fit) <- sites
  return(fit)
}
```
5. Open another new file, enter the code below, and then save the file as "03_run_forecast.R"
```{r}
##' run aquatic forecast into the future
##' @param model site-specific list of forecast models
##' @param met_forecast weather forecast dataframe
##' @param site_data dataframe of site metadata
##' @return dataframe in EFI standard format
run_forecast <- function(model,met_forecast,site_data){
  
  forecast <- NULL
  sites <- names(model)
  
  for(i in 1:length(sites)){
    
    # Get site information for elevation
    site_info <- site_data %>% filter(field_site_id == sites[i]) 
    
    met_future_site <- met_future |> 
      filter(site_id == sites[i])
    
    if(!is.null(model[[i]])){
      
      #use model to forecast water temperature for each ensemble member
      forecasted_temperature <- predict(model[[i]],met_future_site)
      
      #use forecasted temperature to predict oyxgen by assuming that oxygen is saturated.  
      forecasted_oxygen <- rMR::Eq.Ox.conc(forecasted_temperature, 
                                           elevation.m = site_info$field_mean_elevation_m, 
                                           bar.press = NULL, 
                                           bar.units = NULL,
                                           out.DO.meas = "mg/L",
                                           salinity = 0, 
                                           salinity.units = "pp.thou")
      ## organize outputs
      temperature <- tibble(datetime = met_future_site$datetime,
                            site_id = sites[i],
                            parameter = met_future_site$parameter,
                            prediction = forecasted_temperature,
                            variable = "temperature")
      
      oxygen <- tibble(datetime = met_future_site$datetime,
                       site_id = sites[i],
                       parameter = met_future_site$parameter,
                       prediction = forecasted_oxygen,
                       variable = "oxygen")
      
      
      #Build site level dataframe.
      forecast <- dplyr::bind_rows(forecast, temperature, oxygen)
      
    }
    
  }
  
  ## reorganize into EFI standard
  forecast <- forecast |> 
    mutate(reference_datetime = forecast_date) |>
    select(datetime, reference_datetime, site_id, variable, parameter, prediction)
  
  return(forecast)
}