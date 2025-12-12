library(ncdf4)
library(reshape2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(Metrics)

wd <- "~/Documents/isi_cal/"
input <- "ISIMIP_Local_Lakes/LocalLakes"

dir <- list.dirs(path=paste0(wd,input))
lakes <- unlist(strsplit(dir, "/home/dmercado/Documents/isi_cal/ISIMIP_Local_Lakes/LocalLakes/"))
lakes <- lakes[lakes != ""]
lakes <- lakes[2:length(lakes)]

setwd(paste0(wd))

#c <- 1
years_range <- 1901:2019
escenario <- "historical"
model <- "GSWP3-W5E5"

# NSE
nse <- function(obs, pred) {
  numerator <- sum((obs - pred)^2)
  denominator <- sum((obs - mean(obs))^2)
  return(1 - (numerator / denominator))
}
#BIAS
bias <- function(obs, sim){return(mean(sim - obs, na.rm = TRUE))}
# RMSE
rmse <- function(obs, pred) {return(sqrt(mean((obs - pred)^2)))}


for (cali in c("calibrated", "uncalibrated")){
  metrics <- c();
  for (c in 1:length(lakes)){
    
    #observed temperature  
    dir_temp <- dir[c+1]
    temp_obs <- read.csv(paste0(dir_temp, "/", lakes[c],"_temp_daily.csv"))
    temp_obs <- temp_obs[c("TIMESTAMP", "DEPTH", "WTEMP")]
    temp_obs$TIMESTAMP <- paste(as.Date(as.character(temp_obs$TIMESTAMP), format="%Y%m%d"), "00:00:00")
    temp_obs$DEPTH <- temp_obs$DEPTH*(-1)
    
    temp_obs <- temp_obs[temp_obs$TIMESTAMP> as.Date("1995-01-01"),]
    
    temp_obs <- na.omit(temp_obs)
    
    # Convert depths to positive
    temp_obs <- temp_obs %>%
      mutate(DEPTH_POS = abs(DEPTH))
    
    temp_obs_binned <- temp_obs %>%
      mutate(depth_bin = floor(DEPTH_POS)) %>%
      group_by(TIMESTAMP, depth_bin) %>%
      summarise(mean_temp = mean(WTEMP, na.rm = TRUE), .groups = "drop")
    temp_obs_binned$depth_bin <- -1*temp_obs_binned$depth_bin
    
    #write.table(temp_obs, file = paste0(wd, "calibrated/temp.obs"),
    #            quote = F, row.names = F, col.names = F)
    
    # Load NetCDF file
    ncfile <- nc_open(paste0(cali, "/run/", escenario,"/",model,"/out_",tolower(model),"_",escenario,"_", tolower(lakes[c]),"_",years_range[1],"_",years_range[length(years_range)],".nc"))
    
    # Explore variable names
    print(ncfile)
    
    # Extract variables
    time <- ncvar_get(ncfile, "time")  # typically in days
    time <- time/86400 #converto to days
    depth <- ncvar_get(ncfile, "z")    # depth levels
    temp <- ncvar_get(ncfile, "temp")  # 2D matrix: depth x time
    
    # Convert time if needed (e.g., days since a reference date)
    time_units <- ncatt_get(ncfile, "time", "units")$value
    print(time_units)  # e.g., "days since 2000-01-01"
    
    # Convert to POSIXct
    origin <- sub("seconds since ", "", time_units)
    time_dates <- as.Date(time, origin = origin)
    
    # Close NetCDF file
    nc_close(ncfile)
    
    # Make a dataframe
    temp_df <- melt(temp)
    colnames(temp_df) <- c("DepthIndex", "TimeIndex", "Temperature")
    
    # Add actual time and depth
    temp_df$Time <- time_dates[temp_df$TimeIndex]
    temp_df$Depth <- depth[temp_df$DepthIndex]
    
    temp_df <- temp_df[c("Temperature", "Time", "Depth")]
    
    temp_df <- temp_df %>%
      mutate(DEPTH_POS = abs(Depth))
    
    temp_df_binned <- temp_df %>%
      mutate(depth_bin_df = floor(DEPTH_POS)) %>%
      group_by(Time, depth_bin_df) %>%
      summarise(mean_temp = mean(Temperature, na.rm = TRUE), .groups = "drop")
    temp_df_binned$depth_bin_df <- -1*temp_df_binned$depth_bin_df
    
    temp_df_binned$Time <- as.Date(temp_df_binned$Time)
    temp_obs_binned$TIMESTAMP <- as.Date(temp_obs_binned$TIMESTAMP)
    
    # Rename columns to have the same names
    temp_df_binned <- temp_df_binned %>%
      rename(date = Time, depth = depth_bin_df, temp_model = mean_temp)
    
    temp_obs_binned <- temp_obs_binned %>%
      rename(date = TIMESTAMP, depth = depth_bin, temp_obs = mean_temp)
    
    # Merge by date and depth
    merged_df <- full_join(temp_df_binned, temp_obs_binned, by = c("date", "depth"))
    
    merged_df <- merged_df[merged_df$date> as.Date("1995-01-01"),]
    
    matched_data <- na.omit(merged_df)
    
    metrics <- c(metrics, rmse(matched_data$temp_obs, matched_data$temp_model)) 
    metrics <- c(metrics, nse(matched_data$temp_obs, matched_data$temp_model))
    metrics <- c(metrics, cor(matched_data$temp_obs, matched_data$temp_model))
    metrics <- c(metrics, bias(matched_data$temp_obs, matched_data$temp_model))
  }
  
  save_metrics <- data.frame(model="GOTM", lake=sort(rep(lakes, 4)),
                             cali=cali, 
                             metric=rep(c("rmse", "nse", "r", "bias"),length(lakes)),
                             value=metrics)
  write.csv(save_metrics, file=paste0(cali,"/",cali,"_performance.csv"),
            quote = F, row.names = F)
  
}
