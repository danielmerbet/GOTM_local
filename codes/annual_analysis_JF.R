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

# function to calculate annual variables for lake calibration project and
# write them to a csv file
#
# arguments:
# - ncdf: the LakeEnsemblR ncdf file with the model output
# - hypos: hypsograph of the lake with columns Depth_meter and Area_meterSquared
# - lat: latitude of the lake
# - lake: lake name
# - gcm: gcm name
# - scen: scenario name
# - cali: calibrated or uncalibrated
# - out_f: output folder
# - spinup: spinup period for the lake

calc_vars <- function(ncdf, hypso, lat, lake, gcm, scen, cali, out_f = ".",
                      spinup) {
  
  temp <- load_var(ncdf, var = "temp", print = FALSE)
  
  temp <- reshape2::melt(temp, id.vars = colnames(temp[[1]])) |>
    mutate(year = year(datetime)) |> rename(model = "L1") |>
    pivot_longer(c(-datetime, -model)) |>
    mutate(depth = as.numeric(str_remove(name, "wtr_"))) |>
    rename(wtemp = 'value') |> select(model, datetime, depth, wtemp) |>
    filter(model != "Obs")
  
  ice <- load_var(ncdf, var = "ice_height", print = FALSE)
  
  ice <- reshape2::melt(ice, id.vars = colnames(ice[[1]])) |>
    rename(model = "L1") |> mutate(year = year(time))|>
    filter(model != "Obs")
  
  qh <- load_var(ncdf, var = "q_lat", print = FALSE)
  
  qh <- reshape2::melt(qh, id.vars = colnames(qh[[1]])) |>
    rename(model = "L1") |> mutate(year = year(time))|>
    filter(model != "Obs")
  
  qs <- load_var(ncdf, var = "q_sens", print = FALSE)
  
  qs <- reshape2::melt(qs, id.vars = colnames(qs[[1]])) |>
    rename(model = "L1") |> mutate(year = year(time))|>
    filter(model != "Obs")
  
  bottemp_max <- temp |> mutate(year = year(datetime),
                                doy = yday(datetime)) |>
    filter(!is.na(wtemp) & !is.na(depth)) |>
    group_by(year, doy, model) |>
    filter(depth == max(depth)) |> group_by(year, model) |>
    reframe(bottemp_max = max(wtemp))
  
  bottemp_mean <- temp |> mutate(year = year(datetime),
                                 doy = yday(datetime)) |>
    filter(!is.na(wtemp) & !is.na(depth)) |>
    group_by(year, doy, model) |>
    filter(depth == max(depth)) |> group_by(year, model) |>
    reframe(bottemp_mean = mean(wtemp))
  
  bottemp_min <- temp |> mutate(year = year(datetime),
                                doy = yday(datetime)) |>
    filter(!is.na(wtemp) & !is.na(depth)) |>
    group_by(year, doy, model) |>
    filter(depth == max(depth)) |> group_by(year, model) |>
    reframe(bottemp_min = min(wtemp))
  
  surftemp_max <- temp |> mutate(year = year(datetime),
                                 doy = yday(datetime)) |>
    filter(!is.na(wtemp) & !is.na(depth)) |>
    group_by(year, doy, model) |>
    filter(depth == 0) |> group_by(year, model) |>
    reframe(surftemp_max = max(wtemp))
  
  surftemp_mean <- temp |> mutate(year = year(datetime),
                                  doy = yday(datetime)) |>
    filter(!is.na(wtemp) & !is.na(depth)) |>
    group_by(year, doy, model) |>
    filter(depth == 0) |> group_by(year, model) |>
    reframe(surftemp_mean = mean(wtemp))
  
  surftemp_min <- temp |> mutate(year = year(datetime),
                                 doy = yday(datetime)) |>
    filter(!is.na(wtemp) & !is.na(depth)) |>
    group_by(year, doy, model) |>
    filter(depth == 0) |> group_by(year, model) |>
    reframe(surftemp_min = min(wtemp))
  
  heat_mean <- temp |> filter(!is.na(wtemp) & !is.na(depth)) |> mutate(year = year(datetime),
                                                                       doy = yday(datetime)) |>
    group_by(year, doy, model) |>
    reframe(heat = internal.energy(wtemp, depth, hypso$Area_meterSquared,
                                   hypso$Depth_meter)) |>
    group_by(year, model) |>
    reframe(heat_mean = mean(heat))
  
  ice_start <- ice |> mutate(year = year(time) + ifelse(month(time) %in% 10:12, 1, 0)) |>
    mutate(doy = as.numeric(time - as.POSIXct(paste0(year-1,"-10-01 00:00:00"),
                                              tz = "UTC"))/(24*60*60) + 1) |>
    filter(!is.na(ice_height)) |> group_by(year, model) |>
    reframe(ice_start = doy[min(which(ice_height>0))])
  
  ice_end <- ice |> mutate(year = year(time) + ifelse(month(time) %in% 10:12, 1, 0)) |>
    mutate(doy = as.numeric(time - as.POSIXct(paste0(year-1,"-10-01 00:00:00"),
                                              tz = "UTC"))/(24*60*60) + 1) |>
    filter(!is.na(ice_height)) |> group_by(year, model) |>
    reframe(ice_end = doy[max(which(ice_height>0))])
  
  ice_sum <- ice |> mutate(year = year(time) + ifelse(month(time) %in% 10:12, 1, 0)) |>
    mutate(doy = as.numeric(time - as.POSIXct(paste0(year-1,"-10-01 00:00:00"),
                                              tz = "UTC"))/(24*60*60) + 1) |>
    filter(!is.na(ice_height)) |> group_by(year, model) |>
    reframe(ice_sum = sum(ice_height>0))
  
  latentheatf_mean <- qh |> mutate(year = year(time)) |> group_by(year, model) |>
    reframe(latentheatf_mean = mean(q_lat, na.rm = TRUE))
  
  mixeddepth_mean <- temp |> mutate(year = year(datetime),
                                    doy = yday(datetime),
                                    dens = sw_dens(S = 0, t = wtemp, method = "Chen")) |>
    filter(!is.na(dens) & !is.na(depth)) |> group_by(year, doy, model) |>
    reframe(dens = dens,
            depth = depth,
            dens_diff = dens - dens[depth == 0],) |> group_by(year, doy, model) |>
    reframe(mixeddepth = md_fun(dens_diff, depth)) |> group_by(year, model) |>
    reframe(mixeddepth_mean = mean(mixeddepth, na.rm = TRUE))
  
  sensheatf_mean <- qs |> mutate(year = year(time)) |> group_by(year, model) |>
    reframe(sensheatf_mean = mean(q_sens, na.rm = TRUE))
  
  strat_start <- temp |>
    mutate(lat = lat,
           year = ifelse(lat > 0, year(datetime),
                         year(datetime) + ifelse(month(datetime) %in% 10:12, 1, 0)),
           doy = ifelse(lat > 0, yday(datetime),
                        as.numeric(datetime - as.POSIXct(paste0(year-1,"-10-01 00:00:00"),
                                                         tz = "UTC"))/(24*60*60)),
           dens = sw_dens(S = 0, t = wtemp, method = "Chen")) |>
    filter(!is.na(dens) & !is.na(depth)) |> group_by(year, doy, model) |>
    reframe(dens_diff = dens[depth == max(depth)] - dens[depth == 0]) |>
    group_by(year, model) |>
    reframe(strat_start = doy[min(which(dens_diff > 0.1))])
  
  strat_end <- temp |>
    mutate(lat = lat,
           year = ifelse(lat > 0, year(datetime),
                         year(datetime) + ifelse(month(datetime) %in% 10:12, 1, 0)),
           doy = ifelse(lat > 0, yday(datetime),
                        as.numeric(datetime - as.POSIXct(paste0(year-1,"-10-01 00:00:00"),
                                                         tz = "UTC"))/(24*60*60)),
           dens = sw_dens(S = 0, t = wtemp, method = "Chen")) |>
    filter(!is.na(dens) & !is.na(depth)) |> group_by(year, doy, model) |>
    reframe(dens_diff = dens[depth == max(depth)] - dens[depth == 0]) |>
    group_by(year, model) |>
    reframe(strat_end = doy[max(which(dens_diff > 0.1))])
  
  strat_sum <- temp |>
    mutate(lat = lat,
           year = ifelse(lat > 0, year(datetime),
                         year(datetime) + ifelse(month(datetime) %in% 10:12, 1, 0)),
           doy = ifelse(lat > 0, yday(datetime),
                        as.numeric(datetime - as.POSIXct(paste0(year-1,"-10-01 00:00:00"),
                                                         tz = "UTC"))/(24*60*60)),
           dens = sw_dens(S = 0, t = wtemp, method = "Chen")) |>
    filter(!is.na(dens) & !is.na(depth)) |> group_by(year, doy, model) |>
    reframe(dens_diff = dens[depth == max(depth)] - dens[depth == 0]) |>
    group_by(year, model) |>
    reframe(strat_sum = sum(dens_diff > 0.1))
  
  stratstrength_mean <- temp |> mutate(year = year(datetime),
                                       doy = yday(datetime),
                                       dens = sw_dens(S = 0, t = wtemp, method = "Chen")) |>
    filter(!is.na(dens) & !is.na(depth)) |> group_by(year, doy, model) |>
    reframe(dens_diff = dens[depth == max(depth)] - dens[depth == 0]) |>
    group_by(year, model) |>
    reframe(stratstrength_mean = mean(dens_diff[dens_diff > 0.1]))
  
  year_min <- min(ice_start$year) + 1 + spinup
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
                     pivot_longer(cols = c(-year, -model)) |>
                     mutate(gcm = gcm, scenario = scen, cali = cali, lake = lake) |>
                     select(year, model, scenario, gcm, lake, cali, name, value) |>
                     filter(year %in% year_min:year_max) |>
                     mutate(value = round(value, 5), model = paste0(model, "-LER")))
  
  
  ret <- lapply(unique(out$model), function(m){
    out |> filter(model == m) |>
      write.csv(file = file.path(out_f, paste0(m, "_", lake, "_",
                                               scen, "_", gcm, "_",
                                               cali, ".csv")),
                quote = FALSE, row.names = FALSE)
  })
  
  gc()
  return(out)
}
