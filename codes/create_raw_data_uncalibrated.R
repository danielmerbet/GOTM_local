library(ncdf4)
library(reshape)
library(dplyr)
library(data.table)

wd <- "~/Documents/isi_cal/"
input <- "ISIMIP_Local_Lakes/LocalLakes"

dir <- list.dirs(path=paste0(wd,input))
lakes <- unlist(strsplit(dir, "/home/dmercado/Documents/isi_cal/ISIMIP_Local_Lakes/LocalLakes/"))
lakes <- lakes[lakes != ""]
lakes <- lakes[2:length(lakes)]

setwd(wd)

phase <- "isimip3b"

for (escenario in c("picontrol")){ #"historical", "ssp126", "ssp370", "ssp585"
  print(escenario)
  if (escenario=="historical"){
    years_range <- 1850:2014
  }else if (escenario=="picontrol"){
    years_range <- 1850:2100
  }else{
    years_range <- 2015:2100
  }
  
  for (model in c("GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL")){
    print(model)
    for (c in 1:73){
      # Load NetCDF file
      if (lakes[c] %in% c("Kivu", "Tahoe")){
        next
      }
      
      ncfile <- nc_open(paste0("uncalibrated/run/",phase,"/", escenario, "/",model,"/out_",tolower(model),"_",escenario,"_", tolower(lakes[c]),"_",years_range[1],"_",years_range[length(years_range)],".nc"))
      
      # Explore variable names
      #print(ncfile)
      
      # Extract variables
      time <- ncvar_get(ncfile, "time")  # typically in days
      time <- time/86400 #converto to days
      depth <- ncvar_get(ncfile, "z")    # depth levels
      temp <- ncvar_get(ncfile, "temp")  # 2D matrix: depth x time
      qs <- ncvar_get(ncfile, "qh")  # 2D matrix: depth x time
      ql <- ncvar_get(ncfile, "qe")  # 2D matrix: depth x time
      ice <- ncvar_get(ncfile, "Hice")  # 2D matrix: depth x time
      
      # Convert time if needed (e.g., days since a reference date)
      time_units <- ncatt_get(ncfile, "time", "units")$value
      #print(time_units)  # e.g., "days since 2000-01-01"
      
      # Convert to POSIXct
      origin <- sub("seconds since ", "", time_units)
      time_dates <- as.Date(time, origin = origin)
      
      # Close NetCDF file
      nc_close(ncfile)
      
      # Make a dataframe
      ######################################################################
      ######################################################################
      #TEMPERATURE
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
      beep()
      temp_df_binned$Time <- as.Date(temp_df_binned$Time)
      colnames(temp_df_binned) <- c("datetime", "depth", "wtemp")
      
      temp_df_binned <- temp_df_binned[temp_df_binned$datetime> as.Date(paste0(years_range[1],"-01-01")),]
      
      outfile <- paste0("uncalibrated/raw_data/", escenario, "/", model, "/temp_",tolower(lakes[c]), ".csv")
      
      fwrite(temp_df_binned, file = outfile, quote = FALSE, row.names = FALSE)
      
      #write.csv(temp_df_binned, file = outfile,
      #          quote = F, row.names = F)
      ######################################################################
      ######################################################################
      #SENSIBLE HEAT
      qs_df <- melt(qs)
      colnames(qs_df) <- c("TimeIndex", "q_sens")
      
      # Add actual time
      qs_df$Time <- time_dates[qs_df$TimeIndex]
      qs_df <- qs_df[c("Time", "q_sens")]
      colnames(qs_df) <- c("datetime", "q_sens")
      qs_df <- qs_df[qs_df$datetime> as.Date(paste0(years_range[1],"-01-01")),]
      write.csv(qs_df, file = paste0("uncalibrated/raw_data/", escenario, "/", model, "/qs_",tolower(lakes[c]), ".csv"),
                quote = F, row.names = F)
      ######################################################################
      ######################################################################
      #LATENT HEAT
      ql_df <- melt(ql)
      colnames(ql_df) <- c("TimeIndex", "q_lat")
      
      # Add actual time
      ql_df$Time <- time_dates[ql_df$TimeIndex]
      ql_df <- ql_df[c("Time", "q_lat")]
      colnames(ql_df) <- c("datetime", "q_lat")
      ql_df <- ql_df[ql_df$datetime> as.Date(paste0(years_range[1],"-01-01")),]
      write.csv(ql_df, file = paste0("uncalibrated/raw_data/", escenario, "/", model, "/ql_",tolower(lakes[c]), ".csv"),
                quote = F, row.names = F)
      ######################################################################
      ######################################################################
      #ICE THICKNESS
      ice_df <- melt(ice)
      colnames(ice_df) <- c("TimeIndex", "ice_height")
      
      # Add actual time
      ice_df$Time <- time_dates[ice_df$TimeIndex]
      ice_df <- ice_df[c("Time", "ice_height")]
      colnames(ice_df) <- c("datetime", "ice_height")
      ice_df <- ice_df[ice_df$datetime> as.Date(paste0(years_range[1],"-01-01")),]
      
      write.csv(ice_df, file = paste0("uncalibrated/raw_data/", escenario, "/", model, "/ice_",tolower(lakes[c]), ".csv"),
                quote = F, row.names = F)
      
      beep()          
      print(lakes[c])
    }
  }
}





