library(ncdf4)
library(reshape2)
library(ggplot2)

wd <- "~/Documents/isi_cal/"
input <- "ISIMIP_Local_Lakes/LocalLakes"

dir <- list.dirs(path=paste0(wd,input))
lakes <- unlist(strsplit(dir, "/home/dmercado/Documents/isi_cal/ISIMIP_Local_Lakes/LocalLakes/"))
lakes <- lakes[lakes != ""]
lakes <- lakes[2:length(lakes)]

setwd(paste0(wd,"/run/"))

#c <- 1

for (c in 1:length(lakes)){
  # Load NetCDF file
  ncfile <- nc_open(paste0("1st_run/output_", tolower(lakes[c]),".nc"))
  
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
  
  #Plot
  p_save <- ggplot(temp_df, aes(x = Time, y = Depth, fill = Temperature)) +
    geom_tile() +
    #scale_y_reverse() +  # Depth increases downward
    scale_fill_viridis_c(name = "Temp (°C)") +
    labs(title = "Lake Water Temperature Profile",
         x = "Date", y = "Depth (m)") +
    theme_minimal()
  
  pdf(paste0("1st_plot/1stemp_", tolower(lakes[c]),".pdf"), width = 10, height = 6)
  print(p_save)
  dev.off()
  
}

