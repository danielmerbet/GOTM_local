library(ggplot2)
library(dplyr)
library(rLakeAnalyzer)
library(tidyr)
library(RSQLite)
library(lubridate)

wd <- "~/Documents/isi_cal/"
input <- "ISIMIP_Local_Lakes/LocalLakes"

dir <- list.dirs(path=paste0(wd,input))

lakes <- unlist(strsplit(dir, "/home/dmercado/Documents/isi_cal/ISIMIP_Local_Lakes/LocalLakes/"))
lakes <- lakes[lakes != ""]
lakes <- lakes[2:length(lakes)]

#lake characteristics
lake_char <- read.csv(paste0(wd,"LakeCharacteristics.csv"))

#cc function
calc_cc_dewp <- function(date, airt, relh, swr, lat, lon, elev) {
  # Load required libraries
  library(lubridate)
  
  # Convert date to POSIXct
  date <- as.POSIXct(strptime(date, format = "%Y-%m-%d %H:%M:%S"))
  
  # Extend to next day minus one hour
  end_date <- (tail(date, 1) + days(1)) - hours(1)
  date <- seq.POSIXt(from = head(date, 1), to = end_date, by = "hour")
  
  yday <- as.numeric(format(date, "%j"))  # Day of year
  hour <- as.numeric(format(date, "%H"))
  hour[hour == 0] <- 24
  
  # Local standard meridian
  stdmer <- seq(-180, 165, by = 15)
  Lsm <- stdmer[which.min(abs(lon - stdmer))]  # Closest 15-degree meridian
  
  Hsc <- 1390         # Solar constant (W/m2)
  cd <- 0.06          # Dust coefficient
  Rg <- 0.045         # Reflectivity of ground (extended mixed forest)
  pi_val <- pi
  
  theta <- lat * pi_val / 180   # Latitude in radians
  
  # Earth-Sun distance
  r <- 1 + 0.017 * cos((2 * pi_val / 365) * (186 - yday))
  
  # Declination of the sun
  d <- 23.45 * pi_val / 180 * cos((2 * pi_val / 365) * (172 - yday))
  
  # Time correction factor
  dts <- (1 / 15) * (Lsm - lon)
  
  value <- sin(theta) * sin(d) / (cos(theta) * cos(d))
  value[value > 1] <- 1
  value[value < -1] <- -1
  
  tss <- (12 / pi_val) * acos(-value) + dts + 12   # Sunset time
  tsu <- -tss + 2 * dts + 24                       # Sunrise time
  
  gamma <- rep(0, length(tss))
  gamma[hour > tsu & hour < tss] <- 1
  
  # Calculate hb and he
  hb <- he <- numeric(length(hour))
  hb1 <- he1 <- (pi_val / 12) * (hour - dts)
  
  hb1[hour <= 12] <- hb1[hour <= 12] + pi_val
  hb1[hour > 12] <- hb1[hour > 12] - pi_val
  hb <- hb1
  hb[hb > 2 * pi_val] <- hb[hb > 2 * pi_val] - 2 * pi_val
  hb[hb < 0] <- hb[hb < 0] + 2 * pi_val
  
  he1[hour <= 12] <- he1[hour <= 12] + pi_val
  he1[hour > 12] <- he1[hour > 12] - pi_val
  he <- he1
  he[he > 2 * pi_val] <- he[he > 2 * pi_val] - 2 * pi_val
  he[he < 0] <- he[he < 0] + 2 * pi_val
  
  # Extraterrestrial radiation
  Ho <- Hsc / r^2 * (sin(theta) * sin(d) +
                       (12 / pi_val) * cos(theta) * cos(d) * (sin(he) - sin(hb))) * gamma
  
  # Radiation scattering and absorption
  w <- (he + hb) / 2
  alpha1 <- abs(sin(theta) * sin(d) + cos(theta) * cos(d) * cos(w))
  alpha <- atan(alpha1 / sqrt(1 - alpha1^2))  # Solar altitude
  
  # Optical air mass
  theta_am1 <- ((288 - 0.0065 * elev) / 288)^5.256
  theta_am2 <- sin(alpha) + 0.15 * ((alpha * 180 / pi_val) + 3.855)^(-1.253)
  theta_am <- theta_am1 / theta_am2
  
  # Dewpoint temperature
  dewt_daily <- 243.04 * (log(relh / 100) + (17.625 * airt / (243.04 + airt))) /
    (17.625 - log(relh / 100) - (17.625 * airt / (243.04 + airt)))
  
  # Expand daily dewpoint to hourly
  dewt_hourly <- rep(dewt_daily, each = 24)
  
  # Precipitable water content
  Pwc <- 0.85 * exp(0.11 + 0.0614 * dewt_hourly)
  
  # Atmospheric transmission coefficients
  a2 <- exp(-(0.465 + 0.134 * Pwc) * (0.179 + 0.421 * exp(-0.721 * theta_am)) * theta_am)
  a1 <- exp(-(0.465 + 0.134 * Pwc) * (0.129 + 0.171 * exp(-0.88 * theta_am)) * theta_am)
  at <- (a2 + 0.5 * (1 - a1 - cd)) / (1 - 0.5 * Rg * (1 - a1 - cd))
  
  Ho <- at * Ho
  Ho[Ho < 0] <- 1
  
  # Daily mean Ho
  Ho_daily <- rowMeans(matrix(Ho, ncol = 24, byrow = TRUE), na.rm = TRUE)
  
  ccsim <- rep(NA_real_, length(Ho_daily))
  for (i in seq_along(Ho_daily)) {
    if (Ho_daily[i] > swr[i]) {
      ccsim[i] <- sqrt((1 - (swr[i] / Ho_daily[i])) / 0.65)
    }
  }
  
  ccsim[ccsim > 1] <- 1
  
  # Handle edge NAs
  if (is.na(ccsim[1])) ccsim[1] <- 0.5
  if (is.na(ccsim[length(ccsim)])) ccsim[length(ccsim)] <- 0.5
  
  # Interpolate NA values
  if (any(is.na(ccsim))) {
    ccsim <- approx(seq_along(ccsim), ccsim, xout = seq_along(ccsim))$y
  }
  
  ccsim[ccsim > 1] <- 1
  
  return(list(ccsim = ccsim, dewt_daily = dewt_daily))
}

phase <- "isimip3b"
for (escenario in c("ssp126", "ssp370", "ssp585")){ #"historical", "picontrol", 
  #escenario <- "picontrol"
  if (escenario=="historical"){
    years_range <- 1850:2014
  }else if (escenario=="picontrol"){
    years_range <- 1850:2100
  }else{
    years_range <- 2015:2100
  }
  
  for (model in c("GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL")){
    #if (escenario=="historical" & model=="GFDL-ESM4") next
    for (c in 1:73){
      #c <- c+1
      #if (c==1){next}
      
      dir_temp <- dir[c+1]
      
      #observed temperature  
      #temp_obs <- read.csv(paste0(dir_temp, "/", lakes[c],"_temp_daily.csv"))
      #temp_obs <- temp_obs[c("TIMESTAMP", "DEPTH", "WTEMP")]
      #temp_obs$TIMESTAMP <- paste(as.Date(as.character(temp_obs$TIMESTAMP), format="%Y%m%d"), "00:00:00")
      #temp_obs$DEPTH <- temp_obs$DEPTH*(-1)
      
      #temp_obs <- temp_obs[temp_obs$TIMESTAMP> as.Date("1995-01-01"),]
      
      #temp_obs <- na.omit(temp_obs)
      
      #write.table(temp_obs, file = paste0(wd, "calibrated/temp.obs"),
      #            quote = F, row.names = F, col.names = F)
      
      #meteorology
      if (escenario=="picontrol"){
        clim <- paste0(wd, "ISIMIP_Local_Lakes/",phase,"/",escenario,"/",model,"/",tolower(model),"_",escenario,"_bias-adjusted_", tolower(lakes[c]),"_daily_1601_",years_range[length(years_range)],".csv")
      }else{
        clim <- paste0(wd, "ISIMIP_Local_Lakes/",phase,"/",escenario,"/",model,"/",tolower(model),"_",escenario,"_bias-adjusted_", tolower(lakes[c]),"_daily_",years_range[1],"_",years_range[length(years_range)],".csv")
      }
      clim_temp <- read.csv(clim)
      
      #select lake characteristics
      row_char <- which(lake_char$Lake.Name.Folder==tolower(lakes[c]))
      
      #v, u, ps, tas, dew, cc, rsds, pr
      cc_dew <- calc_cc_dewp(paste(clim_temp$time, "00:00:00"),
                             clim_temp$tas-273.15, 
                             clim_temp$hurs, clim_temp$rsds,
                             lake_char$latitude.dec.deg[row_char], 
                             lake_char$longitude.dec.deg[row_char],
                             lake_char$elevation.m[row_char])
      
      clim_save <- data.frame(date=paste(clim_temp$time, "12:00:00"), 
                              v=clim_temp$sfcwind, 
                              u=rep(0, nrow(clim_temp)),
                              ps=clim_temp$ps/100, tas=clim_temp$tas-273.15,
                              dew=cc_dew$dewt_daily,
                              cc=cc_dew$ccsim, swr=clim_temp$rsds, 
                              pr=clim_temp$pr*86400/(1000*24*60*60))
      
      if(escenario=="picontrol"){
        clim_save[clim_save$date >= as.Date(paste0(years_range[1],"-01-01")), ]
      }
      
      #add 20 years spin-up
      clim_save$date <- as.Date(clim_save$date)
      spinup_data <- clim_save[clim_save$date >= as.Date(paste0(years_range[1],"-01-01")) &
                                 clim_save$date <= as.Date(paste0(years_range[1]+19, "-12-31")), ]
      spinup_dates <- seq(as.Date(paste0(years_range[1]-20,"-01-01")),
                          as.Date(paste0(years_range[1]-1,"-12-31")),by=1)
      spinup_data <- spinup_data[1:length(spinup_dates),]
      spinup_data$date <- spinup_dates
      clim_spinup <- rbind(spinup_data, clim_save)
      
      clim_spinup$date <- paste0(clim_spinup$date, " 12:00:00")
      
      write.table(clim_spinup, file = paste0(wd, "calibrated/meteo_file.dat"),
                  quote = F, row.names = F, col.names = F)
      
      #hypsometry
      hyp <- paste0(wd, "ISIMIP_Local_Lakes/LocalLakes/", lakes[c],"/", lakes[c], "_hypsometry.csv")
      hyp_temp <- read.csv(hyp)
      
      col1 <- c(nrow(hyp_temp), 0, cumsum(rev(diff(hyp_temp$DEPTH))))
      col2 <- c(3, rev(hyp_temp$BATHYMETRY_AREA))
      
      hyp_save <- data.frame(col1, col2)
      
      write.table(hyp_save, file = paste0(wd, "calibrated/hypsograph.dat"),
                  quote = F, row.names = F, col.names = F)
      
      #set gotm.yaml
      gotm_ymal <- readLines(paste0(wd, "calibrated/gotm.yaml"))
      gotm_ymal[3] <- paste0("   name: ", tolower(lakes[c]))
      gotm_ymal[4] <- paste0("   latitude: ", lake_char$latitude.dec.deg[row_char])
      gotm_ymal[5] <- paste0("   longitude: ", lake_char$longitude.dec.deg[row_char])
      gotm_ymal[6] <- paste0("   depth: ", lake_char$max.depth.m[row_char])
      gotm_ymal[10] <- paste0(paste0("   start: ", clim_spinup$date[1]))
      gotm_ymal[11] <- paste0(paste0("   stop: ", clim_spinup$date[length(clim_spinup$date)]))
      gotm_ymal[14] <- paste0("   nlev: ", round(lake_char$max.depth.m[row_char]*2))
      if (clim_save$tas[1]<0){
        gotm_ymal[26] <- paste0("      t_1: ", 0)
        gotm_ymal[28] <- paste0("      t_2: ", 0)
      }else{
        gotm_ymal[26] <- paste0("      t_1: ", round(clim_save$tas[1], 2))
        gotm_ymal[28] <- paste0("      t_2: ", round(clim_save$tas[1],2 ))
      }
      
      #add calibrated parameters
      setwd(paste0(wd, "run/"))
      db <- dbConnect(RSQLite::SQLite(), dbname = paste0("parsac_db/parsac_",tolower(lakes[c]),".db") )
      #db <- dbConnect(RSQLite::SQLite(), dbname = paste0("parsac.db") )
      myresults <- dbReadTable(db,"results")
      param <- myresults$parameters[dim(myresults)[1]]
      param <- as.numeric(unlist(strsplit(param, ";")))
      
      gotm_ymal[104] <- paste0("         scale_factor: ", round(param[1], 4))
      gotm_ymal[70] <- paste0("         scale_factor: ", round(param[2], 4))
      gotm_ymal[164] <- paste0("      k_min: ", param[3])
      gotm_ymal[134] <- paste0("      constant_value: ", round(param[4], 4))
      gotm_ymal[139] <- paste0("      constant_value: ", round(param[5], 4))
      gotm_ymal[144] <- paste0("      constant_value: ", round(param[6], 4))
      
      setwd(wd)
      
      writeLines(gotm_ymal, con = paste0(wd, "calibrated/gotm.yaml"))
      
      setwd(paste0(wd, "calibrated/"))
      
      system("./gotm")
      
      system(paste0("mv output.nc run/",phase,"/", escenario,"/",model, "/out_",tolower(model),"_",escenario,"_", tolower(lakes[c]),"_",years_range[1],"_",years_range[length(years_range)],".nc" ))
      
      #system("/home/dmercado/miniconda3/bin/parsac calibration run parsac.xml --maxiter 500 -n 10")
      
      #system(paste0("mv parsac.db parsac_db/parsac_", tolower(lakes[c]), ".db" ))
      
      print(paste0("finished ", lakes[c]))
    }
  }
}
