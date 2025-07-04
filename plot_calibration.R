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

rmse_total <- list()
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
  ncfile <- nc_open(paste0("calibrated/run/", escenario,"/out_",tolower(model),"_",escenario,"_", tolower(lakes[c]),"_",years_range[1],"_",years_range[length(years_range)],".nc"))
  
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
  
  # For each row in temp_obs, find the closest depth in temp_df for the same date
  rmse_by_depth <- matched_data %>%
    group_by(depth) %>%
    summarise(
      RMSE = rmse(temp_model, temp_obs),
      .groups = "drop"
    )
  
  rmse_total[[c]] <- rmse_by_depth
  
  matched_with_rmse <- matched_data %>%
    left_join(rmse_by_depth, by = "depth")
  
  if(lakes[c]=="Tahoe"){
    p_save <- ggplot(matched_with_rmse, aes(x = temp_obs, y = temp_model, color = factor(depth))) +
      geom_point(alpha = 0.7) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
      #facet_wrap(~ Depth_obs) +
      labs(
        title = paste0("Observed vs Simulated Water Temperature: ", lakes[c]),
        #subtitle = paste0("RMSE by Depth:\n", 
        #                  paste0("Depth ", rmse_by_depth$depth, ": ", 
        #                         round(rmse_by_depth$RMSE, 2), collapse = " | ")),
        subtitle = paste0("Mean RMSE: ", round(mean(rmse_by_depth$RMSE),2)),
        x = "Observed Temperature (°C)",
        y = "Simulated Temperature (°C)",
        color = "Observed Depth (m)"
      ) +
      theme_minimal()+
      theme(legend.position="none")
  }else{
    p_save <- ggplot(matched_with_rmse, aes(x = temp_obs, y = temp_model, color = factor(depth))) +
      geom_point(alpha = 0.7) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
      #facet_wrap(~ Depth_obs) +
      labs(
        title = paste0("Observed vs Simulated Water Temperature: ", lakes[c]),
        #subtitle = paste0("RMSE by Depth:\n", 
        #                  paste0("Depth ", rmse_by_depth$depth, ": ", 
        #                         round(rmse_by_depth$RMSE, 2), collapse = " | ")),
        subtitle = paste0("Mean RMSE: ", round(mean(rmse_by_depth$RMSE),2)),
        x = "Observed Temperature (°C)",
        y = "Simulated Temperature (°C)",
        color = "Observed Depth (m)"
      ) +
      theme_minimal()
  }
    
  pdf(paste0("run/calibration_plot/", tolower(lakes[c]),".pdf"), width = 10, height = 6)
  print(p_save)
  dev.off()
  
}

save(rmse_total, file=paste0("run/calibration_plot/rmse_lakes.RData"))



#OLD VERSION:
for (c in 60:length(lakes)){
  
  #observed temperature  
  dir_temp <- dir[c+1]
  temp_obs <- read.csv(paste0(dir_temp, "/", lakes[c],"_temp_daily.csv"))
  temp_obs <- temp_obs[c("TIMESTAMP", "DEPTH", "WTEMP")]
  temp_obs$TIMESTAMP <- paste(as.Date(as.character(temp_obs$TIMESTAMP), format="%Y%m%d"), "00:00:00")
  temp_obs$DEPTH <- temp_obs$DEPTH*(-1)
  
  temp_obs <- temp_obs[temp_obs$TIMESTAMP> as.Date("1995-01-01"),]
  
  temp_obs <- na.omit(temp_obs)
  
  #write.table(temp_obs, file = paste0(wd, "calibrated/temp.obs"),
  #            quote = F, row.names = F, col.names = F)
  
  # Load NetCDF file
  ncfile <- nc_open(paste0("calibrated/run/", escenario,"/out_",tolower(model),"_",escenario,"_", tolower(lakes[c]),"_",years_range[1],"_",years_range[length(years_range)],".nc"))
  
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
  
  temp_df$Time <- as.Date(temp_df$Time)
  temp_obs$TIMESTAMP <- as.Date(temp_obs$TIMESTAMP)
  
  # For each row in temp_obs, find the closest depth in temp_df for the same date
  matched_data <- temp_obs %>%
    rowwise() %>%
    mutate(
      match_row = list(
        temp_df %>%
          filter(Time == TIMESTAMP) %>%
          mutate(depth_diff = abs(Depth - DEPTH)) %>%
          slice_min(depth_diff, n = 1)
      )
    ) %>%
    unnest(match_row) %>%
    # Clean and rename columns
    transmute(
      Date = TIMESTAMP,
      Depth_df = Depth,
      Depth_obs = DEPTH,
      Temperature = Temperature,
      WTEMP = WTEMP
    )
  
  rmse_by_depth <- matched_data %>%
    group_by(Depth_obs) %>%
    summarise(
      RMSE = rmse(WTEMP, Temperature),
      .groups = "drop"
    )
  
  rmse_total[[c]] <- rmse_by_depth
  
  matched_with_rmse <- matched_data %>%
    left_join(rmse_by_depth, by = "Depth_obs")
  
  p_save <- ggplot(matched_with_rmse, aes(x = Temperature, y = WTEMP, color = factor(Depth_obs))) +
    geom_point(alpha = 0.7) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    #facet_wrap(~ Depth_obs) +
    labs(
      title = "Simulated vs Observed Temperature",
      subtitle = paste0("RMSE by Depth:\n", 
                        paste0("Depth ", rmse_by_depth$Depth_obs, ": ", round(rmse_by_depth$RMSE, 2), collapse = " | ")),
      x = "Simulated Temperature (°C)",
      y = "Observed Temperature (°C)",
      color = "Observed Depth (m)"
    ) +
    theme_minimal()
  
  pdf(paste0("run/calibration_plot/", tolower(lakes[c]),".pdf"), width = 10, height = 6)
  print(p_save)
  dev.off()
  
}
