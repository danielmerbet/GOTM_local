library(marelac)
library(rLakeAnalyzer)
library(tidyverse)
library(LakeEnsemblR)

# helper function to extract depth with density difference of 0.1
md_fun <- function(dens_diff, depth){
  md <- suppressWarnings(tryCatch({
    approx(dens_diff, depth, 0.1)$y},
    error = function(e) {NA}))
  return(md)
}

homedir <- "~/Documents/isi_cal/"
calibration <- "calibrated"
scenario <- "ssp126" #"historical"
model <- "gotm"
hypsodir <- paste0(homedir, "/ISIMIP_Local_Lakes/hypso/")
lake_info <- read.csv(paste0(homedir,"LakeCharacteristics.csv"))

for (calibration in c("calibrated")){ #"calibrated", 
  for (scenario in c("picontrol")){ #"historical", , "ssp370", "ssp585","ssp126"
    for (gcm in c("GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL")){

      readin <- paste0(homedir, calibration, "/raw_data/", scenario, "/", gcm, "/")
      writeout <- paste0(homedir, calibration, "/annual/", scenario, "/", gcm, "/")
      
      if (file.exists(paste0(writeout, "GOTM_", scenario, "_", tolower(gcm), "_", calibration, ".csv" ))){
        next
      }
      
      setwd(readin)
      for (i in 1:73){
        lake_name <- strsplit(strsplit(list.files(path=readin, pattern="qs")[i], ".csv")[[1]], "qs_")[[1]][2]
        lake_info_temp <- lake_info[lake_info$Lake.Name.Folder==lake_name,]
        
        #open files
        temp <- read.csv(paste0("temp_", lake_name, ".csv"))
        temp <- temp |> mutate(year = year(datetime)) 
        temp$datetime <- as.POSIXct(temp$datetime, tz = "UTC")
        ice <- read.csv(paste0("ice_", lake_name, ".csv"))
        ice <- ice |> mutate(year = year(datetime)) 
        ql <- read.csv(paste0("ql_", lake_name, ".csv"))
        ql <- ql |> mutate(year = year(datetime)) 
        ql$datetime <- as.POSIXct(ql$datetime, tz = "UTC")
        qs <- read.csv(paste0("qs_", lake_name, ".csv"))
        qs <- qs |> mutate(year = year(datetime)) 
        qs$datetime <- as.POSIXct(qs$datetime, tz = "UTC")
        
        #open hypsometry
        hypso <- read.csv(paste0(hypsodir, str_to_title(lake_name), "_hypsometry.csv"))
        
        #calculate annual variables
        bottemp_max <- temp |> mutate(doy = yday(datetime)) |>
          filter(!is.na(wtemp) & !is.na(depth)) |>
          group_by(year, doy) |>
          filter(depth == max(depth)) |> group_by(year) |>
          reframe(bottemp_max = max(wtemp))
        
        bottemp_mean <- temp |> mutate(doy = yday(datetime)) |>
          filter(!is.na(wtemp) & !is.na(depth)) |>
          group_by(year, doy) |>
          filter(depth == max(depth)) |> group_by(year) |>
          reframe(bottemp_mean = mean(wtemp))
        
        bottemp_min <- temp |> mutate(doy = yday(datetime)) |>
          filter(!is.na(wtemp) & !is.na(depth)) |>
          group_by(year, doy) |>
          filter(depth == max(depth)) |> group_by(year) |>
          reframe(bottemp_min = min(wtemp))
        
        surftemp_max <- temp |> mutate(doy = yday(datetime)) |>
          filter(!is.na(wtemp) & !is.na(depth)) |>
          group_by(year, doy) |>
          filter(depth == 0) |> group_by(year) |>
          reframe(surftemp_max = max(wtemp))
        
        surftemp_mean <- temp |> mutate(doy = yday(datetime)) |>
          filter(!is.na(wtemp) & !is.na(depth)) |>
          group_by(year, doy) |>
          filter(depth == 0) |> group_by(year) |>
          reframe(surftemp_mean = mean(wtemp))
        
        surftemp_min <- temp |> mutate(doy = yday(datetime)) |>
          filter(!is.na(wtemp) & !is.na(depth)) |>
          group_by(year, doy) |>
          filter(depth == 0) |> group_by(year) |>
          reframe(surftemp_min = min(wtemp))
        
        heat_mean <- temp |> mutate(doy = yday(datetime)) |>
          filter(!is.na(wtemp) & !is.na(depth)) |> 
          group_by(year, doy) |>
          reframe(heat = internal.energy(wtemp, depth, hypso$BATHYMETRY_AREA,
                                         hypso$DEPTH)) |>
          group_by(year) |>
          reframe(heat_mean = mean(heat))
        
        ice$datetime <- as.POSIXct(ice$datetime, tz = "UTC")
        ice_start <- ice |> mutate(year = year(datetime) + ifelse(month(datetime) %in% 10:12, 1, 0)) |>
          mutate(doy = as.numeric(datetime - as.POSIXct(paste0(year-1,"-10-01 00:00:00"),
                                                        tz = "UTC"))/(24*60*60) + 1) |>
          filter(!is.na(ice_height)) |> group_by(year) |>
          reframe(ice_start = doy[min(which(ice_height>0))])
        
        ice_end <- ice |> mutate(year = year(datetime) + ifelse(month(datetime) %in% 10:12, 1, 0)) |>
          mutate(doy = as.numeric(datetime - as.POSIXct(paste0(year-1,"-10-01 00:00:00"),
                                                        tz = "UTC"))/(24*60*60) + 1) |>
          filter(!is.na(ice_height)) |> group_by(year) |>
          reframe(ice_end = doy[max(which(ice_height>0))])
        
        ice_sum <- ice |> mutate(year = year(datetime) + ifelse(month(datetime) %in% 10:12, 1, 0)) |>
          mutate(doy = as.numeric(datetime - as.POSIXct(paste0(year-1,"-10-01 00:00:00"),
                                                        tz = "UTC"))/(24*60*60) + 1) |>
          filter(!is.na(ice_height)) |> group_by(year) |>
          reframe(ice_sum = sum(ice_height>0))
        
        latentheatf_mean <- ql |> group_by(year) |>
          reframe(latentheatf_mean = mean(q_lat, na.rm = TRUE))
        
        mixeddepth_mean <- temp |> mutate(doy = yday(datetime),
                                          dens = sw_dens(S = 0, t = wtemp, method = "Chen")) |>
          filter(!is.na(dens) & !is.na(depth)) |> group_by(year, doy) |>
          reframe(dens = dens,
                  depth = depth,
                  dens_diff = dens - dens[depth == 0],) |> group_by(year, doy) |>
          reframe(mixeddepth = md_fun(dens_diff, depth)) |> group_by(year) |>
          reframe(mixeddepth_mean = mean(mixeddepth, na.rm = TRUE))
        
        sensheatf_mean <- qs |> group_by(year) |>
          reframe(sensheatf_mean = mean(q_sens, na.rm = TRUE))
        
        lat <- lake_info_temp$latitude.dec.deg
        strat_start <- temp |>
          mutate(lat = lat,
                 year = ifelse(lat > 0, year(datetime),
                               year(datetime) + ifelse(month(datetime) %in% 10:12, 1, 0)),
                 doy = ifelse(lat > 0, yday(datetime),
                              as.numeric(datetime - as.POSIXct(paste0(year-1,"-10-01 00:00:00"),
                                                               tz = "UTC"))/(24*60*60)),
                 dens = sw_dens(S = 0, t = wtemp, method = "Chen")) |>
          filter(!is.na(dens) & !is.na(depth)) |> group_by(year, doy) |>
          reframe(dens_diff = dens[depth == max(depth)] - dens[depth == 0]) |>
          group_by(year) |>
          reframe(strat_start = doy[min(which(dens_diff > 0.1))])
        
        strat_end <- temp |>
          mutate(lat = lat,
                 year = ifelse(lat > 0, year(datetime),
                               year(datetime) + ifelse(month(datetime) %in% 10:12, 1, 0)),
                 doy = ifelse(lat > 0, yday(datetime),
                              as.numeric(datetime - as.POSIXct(paste0(year-1,"-10-01 00:00:00"),
                                                               tz = "UTC"))/(24*60*60)),
                 dens = sw_dens(S = 0, t = wtemp, method = "Chen")) |>
          filter(!is.na(dens) & !is.na(depth)) |> group_by(year, doy) |>
          reframe(dens_diff = dens[depth == max(depth)] - dens[depth == 0]) |>
          group_by(year) |>
          reframe(strat_end = doy[max(which(dens_diff > 0.1))])
        
        strat_sum <- temp |>
          mutate(lat = lat,
                 year = ifelse(lat > 0, year(datetime),
                               year(datetime) + ifelse(month(datetime) %in% 10:12, 1, 0)),
                 doy = ifelse(lat > 0, yday(datetime),
                              as.numeric(datetime - as.POSIXct(paste0(year-1,"-10-01 00:00:00"),
                                                               tz = "UTC"))/(24*60*60)),
                 dens = sw_dens(S = 0, t = wtemp, method = "Chen")) |>
          filter(!is.na(dens) & !is.na(depth)) |> group_by(year, doy) |>
          reframe(dens_diff = dens[depth == max(depth)] - dens[depth == 0]) |>
          group_by(year) |>
          reframe(strat_sum = sum(dens_diff > 0.1))
        
        stratstrength_mean <- temp |> mutate(doy = yday(datetime),
                                             dens = sw_dens(S = 0, t = wtemp, method = "Chen")) |>
          filter(!is.na(dens) & !is.na(depth)) |> group_by(year, doy) |>
          reframe(dens_diff = dens[depth == max(depth)] - dens[depth == 0]) |>
          group_by(year) |>
          reframe(stratstrength_mean = mean(dens_diff[dens_diff > 0.1]))
        
        year_min <- min(ice_start$year) + 1 #+ spinup
        year_max <- max(ice_start$year) - 1
        
        suppressMessages(out <- full_join(bottemp_max, bottemp_mean) |>
                           full_join(bottemp_min) |>
                           full_join(surftemp_max) |> full_join(surftemp_mean) |>
                           full_join(surftemp_min) |> full_join(heat_mean) |>
                           full_join(ice_start) |> full_join(ice_end) |>
                           full_join(ice_sum) |> full_join(latentheatf_mean) |>
                           full_join(mixeddepth_mean) |> full_join(sensheatf_mean) |>
                           full_join(strat_start) |> full_join(strat_end) |>
                           full_join(strat_sum) |> full_join(stratstrength_mean) |>
                           pivot_longer(cols = c(-year)) |>
                           mutate(gcm = gcm, scenario = scenario, cali = calibration, 
                                  lake = lake_name, model = model) |>
                           dplyr::select(year, model, scenario, gcm, lake, cali, name, value) |>
                           filter(year %in% year_min:year_max) |>
                           mutate(value = round(value, 5)))
        
        if (i==1){
          out_total <- out
        }else{
          out_total <- rbind(out_total, out)
        }
        print(paste0("lake ", lake_name, " finished"))
      }
      
      out_total$model <- "GOTM"
      out_total$gcm <- tolower(gcm)
      
      write.csv(out_total, file=paste0(writeout, "GOTM_", scenario, "_", tolower(gcm), "_", calibration, ".csv" ),
                quote = F, row.names = F)
      print(paste(gcm, scenario, calibration, "finished"))
    }
  }
}



