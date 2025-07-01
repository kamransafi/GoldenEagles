### Code to take classified soaring flight, full data at a coarse scale, and emigration dates
### and estimate the effects of soaring flight and daily movement on excursion and emigration timing
### Hester Bronnvik
### hbronnvik@ab.mpg.de
### 25.07.2025

library(geodata)
library(survival)
library(ggfortify)
library(survminer)
library(geosphere)
library(sf)
library(tidyverse)
options(warn = 1)

wind_speed <- function(u,v) {
  return(sqrt(u*u+v*v))
}
# plotting
theme_set(theme_classic()+theme(axis.text = element_text(color = "black", size = 10), 
                                text = element_text(size = 12),
                                panel.background = element_rect(fill = "#F5F5F5"),
                                strip.background = element_rect(fill = "#F5F5F5", color = "#F5F5F5"),
                                strip.text = element_text(face = "bold")))
colfunc <- colorRampPalette(c("#6F101E", "#8F1D28", "#BF413E", "#FACEB6", "white", 
                              "#A4CCE3", "#3D84B4", "#1D5497", "#0A3162"))



### I need the emigration dates and the excursion dates

# dispersal and tagging dates
dates <- readRDS("/home/hbronnvik/Documents/chapter4/dispersal_dates_20250417.rds") %>% 
  mutate(dispersal_date = as.POSIXct(ifelse(did_disperse == F, NA, dispersal_date), tz = "UTC"),
         pfdp = as.numeric(difftime(dispersal_date, deploy_on_timestamp, units = "days")))
# timing of excursions
ex_info <- readRDS("/home/hbronnvik/Documents/chapter4/excursions_250127.rds")
# find the first excursions
first_ex <- ex_info %>% 
  # the first excursion on which a bird traveled at least 10 km
  # and maybe if it is only a couple points it was an error
  filter(excursion_dist >= 10 & excursion_observations >= 3) %>% 
  group_by(individual.id) %>% 
  slice(1) %>% 
  ungroup() %>% 
  left_join(dates %>% dplyr::select(individual.id, deploy_on_timestamp, dispersal_date, did_disperse)) %>% 
  mutate(days_to_leave = round(as.numeric(difftime(timestamp, deploy_on_timestamp, units = "days")))) %>% 
  dplyr::select(individual.id, individual.local.identifier, excursion_id, excursion_observations,
                excursion_duration, excursion_dist, excursion_dist_sum, days_to_leave, did_disperse)


### I need the daily scale movements

# a list of 89 birds and their movements without bursts
locs <- lapply(list.files("/home/hbronnvik/Documents/GE_data/no_burst_in_out", full.names = T), function(l){
  ind <- readRDS(l)
  # remove any >= 1 km inaccuracies
  ind <- ind[which(ind$eobs.horizontal.accuracy.estimate <= 1000),]
  ind <- ind[, c("location.long", "location.lat", "timestamp", "individual.id", "height.above.ellipsoid", "outside")]
  return(ind)
})
# find the 99.99% cutoff for removing highly unlikely outlier heights
cutoff <- locs %>% 
  data.table::rbindlist(fill = T, use.names = T) %>% 
  dplyr::select(height.above.ellipsoid) %>% 
  deframe() %>% 
  quantile(0.9999, na.rm = T) %>% 
  as.numeric()
# calculate individual metrics
daily_move <- lapply(locs, function(x){
  # deal with tags that transmit post-mortem (e.g. 1600584021)
  x <- x[x$timestamp <= x$timestamp[1]+days(365),]
  # the emigration date
  emig <- dates$dispersal_date[dates$individual.id == unique(x$individual.id)]
  # if there is none, the final date
  if(is.na(emig)){
    emig <- max(x$timestamp)
  }
  # add the tagging date
  x$tagging_date <- dates$deploy_on_timestamp[dates$individual.id == unique(x$individual.id)]
  # add a day column
  x$date <- date(x$timestamp)
  # get the time before emigration
  ind <- x[x$date <= date(emig),]
  # group by day and get the maximum height
  ind_heights <- ind %>% 
    drop_na(location.long) %>% 
    # remove any outliers using the setting above
    filter(height.above.ellipsoid < cutoff) %>% 
    group_by(date) %>% 
    summarize(max_height = max(height.above.ellipsoid))
  # now get geodesic distances between midday locations
  ind_dists <- ind %>% 
    # for low-battery tags, we often get one fix at 19:00, but we do want only midday locations
    filter(hour(timestamp) %in% c(12, 13)) %>% 
    group_by(date) %>% 
    # take the hour closest to noon
    mutate(lag = abs(difftime(as.POSIXct(paste(unique(date), "12:00:00"), tz = "UTC"),
                              timestamp, units = "secs"))) %>% 
    slice_min(lag) %>% 
    # break ties
    slice(1) %>% 
    ungroup() %>% 
    # estimate how much change there was between noon locations (in m)
    mutate(interday_dist = distVincentyEllipsoid(cbind(lag(location.long), lag(location.lat)),
                                 cbind(location.long, location.lat))) %>% 
    #clean up
    dplyr::select(individual.id, tagging_date, timestamp, outside, date, interday_dist) %>% 
    rename(midday = timestamp)
  # combine and save
  ind <- ind_heights %>% 
    left_join(ind_dists, by = join_by(date))
  # gives a height and distance for every day of the PFDP
  return(ind)
})
daily_move <- data.table::rbindlist(daily_move)

# see what we have
daily_move %>% 
  group_by(individual.id) %>% 
  summarize(obs = n()) %>% 
  arrange(obs)

daily_move %>% 
  mutate(dst = as.numeric(difftime(midday, tagging_date, units = "days")),
         log_dist = log(interday_dist)) %>% 
  pivot_longer(cols = c(log_dist, max_height)) %>% 
  ggplot(aes(dst, value)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~name, scales = "free_y")

### I need the fine-scale movements
# bursts with thermal soaring split into events and with wind estimates
winds <- lapply(list.files("/home/hbronnvik/Documents/GE_data/wind_estimates_25", full.names = T), function(l){
  ind <- readRDS(l)
  ind <- ind[, c("location.long", "location.lat", "timestamp", "individual.id", 
                 "height_above_ellipsoid", "thermal_event", "soarClust", "thermalClust",
                 "thermal_duration", "vspeed_thermal", "turn_var_thermal", "windX", 
                 "windY", "airSpeed", "vert_speed", "ThermalStrength", "CircRadius", "CircTime")]
  return(ind)
})
# now I need to get per-thermal metrics and average those to the day
soar_metrics <- lapply(winds, function(x){
  # deal with tags that transmit post-mortem (e.g. 1600584021)
  x <- x[x$timestamp <= x$timestamp[1]+days(365),]
  # the emigration date
  emig <- dates$dispersal_date[dates$individual.id == unique(x$individual.id)]
  # if there is none, the final date
  if(is.na(emig)){
    emig <- max(x$timestamp)
  }
  # add the tagging date
  x$tagging_date <- dates$deploy_on_timestamp[dates$individual.id == unique(x$individual.id)]
  # add a day column
  x$date <- date(x$timestamp)
  # get the time before emigration
  ind <- x[x$date <= date(emig)+days(10),]
  # if there are PFDP data, crunch them
  if(nrow(ind) > 0){
    # per thermal measures of flight
    flight <- ind %>% 
      mutate(wind_sp = wind_speed(windX, windY)) %>% 
      group_by(thermal_event) %>% 
      summarize(date = unique(date(timestamp)),
                therm_circ_rad = mean(CircRadius, na.rm = T),
                vspeed_thermal = unique(vspeed_thermal),
                therm_wind_sp = mean(wind_sp, na.rm = T))
    # compress those into days
    daily_flight <- flight %>% 
      group_by(date)%>% 
      summarize(med_circ_rad = median(therm_circ_rad),
                med_vspeed = median(vspeed_thermal),
                med_ws = median(therm_wind_sp)) %>% 
      mutate(individual.id = unique(x$individual.id),
             tagging_date = unique(x$tagging_date),
             dispersal_date = emig, 
             dst = as.numeric(difftime(date, tagging_date, units = "days")))
    return(daily_flight)
  }
})
soar_metrics <- data.table::rbindlist(soar_metrics)
# 43 birds left 
length(unique(soar_metrics$individual.id))

soar_metrics %>% 
  pivot_longer(cols = c(med_circ_rad, med_vspeed, med_ws)) %>% 
  ggplot(aes(dst, value)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~name, scales = "free_y")

# daily_move <- lapply(split(daily_move, daily_move$individual.id), function(d){
#   has_noon <- ifelse(12 %in% unique(hour(d$midday)), T, F)
#   if(has_noon == T){
#     d <- d %>% 
#       filter(hour(midday) == 12)
#   }
#   return(d)
# }) %>% reduce(rbind)

### I need to format those with status and interval
coarse_intervals <- lapply(split(daily_move, daily_move$individual.id), function(x){
  emig <- dates$dispersal_date[dates$individual.id == unique(x$individual.id)]
  excu <- first_ex$days_to_leave[first_ex$individual.id == unique(x$individual.id)]

  x <- x %>% 
    # the start and end of each interval
    mutate(dst = as.numeric(difftime(date, date(tagging_date), units = "days")),
           t1 = dst-1, t2 = dst)
  emig_x <- x %>% 
    mutate(# the status (returns NA if the bird did not emigrate)
      emigration_status = ifelse(date == date(emig), 1, 0))
  # add excursion status if there was an excursion
  if(length(excu) > 0){
    excu_x <- x %>% 
      filter(dst <= excu) %>% 
      # the other status
      mutate(excursion_status = ifelse(dst == excu, 1, 0)) %>% 
      dplyr::select(date, excursion_status)
  }else{
    excu_x <- x %>% 
      # the other status
      mutate(excursion_status = NA) %>% 
      dplyr::select(date, excursion_status)
  }
  x <- emig_x %>% 
    left_join(excu_x, by = join_by(date))
  return(x)
})
coarse_intervals <- data.table::rbindlist(coarse_intervals)
# scale the predictors
coarse_intervals <- coarse_intervals %>% 
  mutate(interday_dist_z = scale(interday_dist)[,1],
         max_height_z = scale(max_height)[,1])

# the fine-scale metrics have the complication of missing data, I need to allow for a ten day gap
soar_intervals <- lapply(split(soar_metrics, soar_metrics$individual.id), function(x){
  emig <- dates$dispersal_date[dates$individual.id == unique(x$individual.id)]
  excu <- first_ex$days_to_leave[first_ex$individual.id == unique(x$individual.id)]
  
  x <- x %>% 
    # the start and end of each interval
    mutate(dst = as.numeric(difftime(date, date(tagging_date), units = "days")),
           t1 = dst-1, t2 = dst)
  emig_x <- x %>% 
    # we carefully kept 10 days past emigration above in case, now trim them off
    filter(date <= emig) %>% 
    mutate(# the status (returns NA if the bird did not emigrate)
      emigration_status = ifelse(date == date(emig), 1, 0)) %>% 
    dplyr::select(date, emigration_status)
  # if we need to deal with a gap during emigration
  if(length(unique(emig_x$emigration_status)) == 1){
    # take the first observation within 10 days of excursion
    after <- x %>% 
      filter(date > emig) %>% 
      slice(1) %>% 
      # set that to the event occurrence
      mutate(emigration_status = 1)
    # add that to the status
    if(nrow(after) > 0){
      emig_x <- x %>% 
        filter(date <= emig) %>% 
        mutate(emigration_status = 0) %>% 
        rbind(after) %>% 
        dplyr::select(date, emigration_status)
    }else{
      # or remove the status along with the 10 day buffer if that buffer was not big enough
      emig_x <- x %>% 
        filter(date <= emig) %>% 
        mutate(emigration_status = NA) %>% 
        dplyr::select(date, emigration_status)
    }
  }
  
  # add excursion status if there was an excursion
  if(length(excu) > 0){
    excu_x <- x %>% 
      filter(dst <= excu) %>% 
      # the other status
      mutate(excursion_status = ifelse(dst == excu, 1, 0)) %>% 
      dplyr::select(date, excursion_status)
    
    # if we need to deal with a gap during that excursion
    if(length(unique(excu_x$excursion_status)) == 1){
      # take the first observation within 10 days of excursion
      after <- x %>% 
        filter(between(dst, excu, excu+10)) %>% 
        slice(1) %>% 
        # set that to the event occurrence
        mutate(excursion_status = 1)
      # add that to the status
      excu_x <- x %>% 
        filter(dst <= excu) %>% 
        mutate(excursion_status = 0) %>% 
        rbind(after) %>% 
        dplyr::select(date, excursion_status)
    }
    
  }else{
    excu_x <- x %>% 
      # the other status
      mutate(excursion_status = NA) %>% 
      dplyr::select(date, excursion_status)
  }
  x <- x %>% 
    left_join(emig_x, by = join_by(date)) %>% 
    left_join(excu_x, by = join_by(date))
  return(x)
})
soar_intervals <- data.table::rbindlist(soar_intervals)

### I need four models of those
# two for the daily metrics
fit_c1 <- coxph(Surv(t1, t2, excursion_status) ~ interday_dist_z+max_height_z,
                data = coarse_intervals)
fit_c2 <- coxph(Surv(t1, t2, emigration_status) ~ interday_dist_z+max_height_z,
                data = coarse_intervals)
# and two for the soaring metrics
fit_f1 <- coxph(Surv(t1, t2, excursion_status) ~ med_circ_rad+med_vspeed+med_ws,
                data = soar_intervals)
fit_f2 <- coxph(Surv(t1, t2, emigration_status) ~ med_circ_rad+med_vspeed+med_ws,
      data = soar_intervals)

# list to process
fits <- list(fit_c1, fit_c2, fit_f1, fit_f2)
names(fits) <- c("fit_c1", "fit_c2", "fit_f1", "fit_f2")

# create a data frame with HRs and CIs
hr_table <- lapply(1:length(fits), function(f){
  x <- fits[[f]]
  data.frame(variable = names(coef(x)),
             HR = round(exp(coef(x)), 3),
             lower = round(exp(confint(x)[, 1]), 3),
             upper = round(exp(confint(x)[, 2]), 3),
             p_value = round(summary(x)$coefficients[, "Pr(>|z|)"], 3),
             model = names(fits)[[f]])
}) %>% reduce(rbind)
print(hr_table)

sd(coarse_intervals$interday_dist, na.rm = T)
sd(coarse_intervals$max_height, na.rm = T)

### I need to make some plots to illustrate these

curve_plts <- lapply(1:length(fits), function(f){
  x <- fits[[f]]
  if(grepl("_c", names(fits)[f])){
    # for coarse scale data (_c), the data are coarse_intervals
    df <- coarse_intervals
  }else{
    df <- soar_intervals
  }
  if(grepl("1", names(fits)[f])){
    # for excursions (event 1)
    col <- "#2966A2"
    samp_n <- sum(df$excursion_status, na.rm = T)
    event <- "excursion"
    if(grepl("_c", names(fits)[f])){
      st <- "Daily movement"
    }else{
      st <- "Soaring performance"
    }
  }else{
    col <- "#A12A30"
    samp_n <- sum(df$emigration_status, na.rm = T)
    event <- "emigration"
    st <- ""
  }
  ## plot the outputs
  # as a curve
  p <- ggsurvplot(surv_fit(x, data = df),
                  palette = col,
                  xlab = "Days since tagging", 
                  ylab = paste0("Probability of\n", event," (n = ", samp_n, ")"),
                  subtitle = st) 
  p <- p$plot + 
    geom_hline(yintercept = 0.5, lty = 2) +
    geom_vline(xintercept = read.table(textConnection(capture.output(survfit(x))),skip=2,header=TRUE)$median,
               lty = 2) +
    scale_y_reverse(expand = c(0,0), limits = c(1.1, -0.2), breaks = c(1, 0.75, 0.5, 0.25, 0),
                    labels = c("0%", "25%", "50%", "75%", "100%")) +
    scale_x_continuous(limits = c(0, 310)) +
    theme_classic() +
    theme(axis.text = element_text(color = "black", size = 10),
          text = element_text(size = 12),
          panel.background = element_rect(fill = "#F5F5F5"),
          strip.background = element_rect(fill = "#F5F5F5", color = "#F5F5F5"),
          strip.text = element_text(face = "bold"),
          legend.position = "none",
          plot.margin=unit(c(0,0,0,1), 'cm'),
          plot.subtitle = element_text(face = "bold"))
  # the coefficients
  surv_coef <- x %>%  
    broom::tidy(conf.int = TRUE, exponentiate = TRUE) %>%  
    dplyr::select(term, estimate, starts_with("conf")) %>% 
    mutate(term = sub("interday_dist", "Daily distance",
                      sub("max_height", "Daily max. height",
                          sub("med_circ_rad", "Median circling radius", 
                              sub("med_ws", "Median wind speed",
                                  sub("med_vspeed", "Median vertical speed",
                                      sub("_z", "", term)))))))
  print(surv_coef)
  
  coefs <- ggplot(surv_coef, aes(estimate, fct_rev(term))) +
    geom_vline(lty = 2, xintercept = 1) +
    geom_pointrange(aes(xmin = conf.low, xmax = conf.high), color = col) +
    # geom_text(aes(label = round(estimate, 4)), vjust=-1) +
    labs(x = "Estimate and 95% hazard ratio CI", y = "")
  res <- ggpubr::ggarrange(p, coefs, ncol = 2)
  return(res)
})

# png(filename = "/home/hbronnvik/Documents/chapter4/independence_25/surv_interval.png",
#     height = 8.3, width = 11, units = "in", res = 300)
ggarrange(curve_plts[[3]], curve_plts[[4]],
          curve_plts[[1]], curve_plts[[2]], 
          nrow = 4, labels = "AUTO")
# dev.off()

# check the model fits
res_plts <- lapply(1:length(fits), function(f){
  x <- fits[[f]]
  if(grepl("_c", names(fits)[f])){
    # for coarse scale data (_c), the data are coarse_intervals
    df <- coarse_intervals
  }else{
    df <- soar_intervals
  }
  if(grepl("1", names(fits)[f])){
    # for excursions (event 1)
    col <- "#2966A2"
    samp_n <- sum(df$excursion_status, na.rm = T)
    event <- "excursion"
    if(grepl("_c", names(fits)[f])){
      st <- "Daily movement"
    }else{
      st <- "Soaring performance"
    }
  }else{
    col <- "#A12A30"
    samp_n <- sum(df$emigration_status, na.rm = T)
    event <- "emigration"
    st <- ""
  }
  len <- ifelse(grepl("_c", names(fits)[f]), 2, 3)
  
  # p-values
  zph <- cox.zph(x)
  print(zph)
  individual_ps <- zph$table[-nrow(zph$table), "p"] %>% 
    as.data.frame() %>% 
    rownames_to_column() %>% 
    rename(name = 1, p_value = 2)
  # individual_ps <- setNames(zph$table[-nrow(zph$table), "p"], 
  #                      rownames(zph$table)[-nrow(zph$table)])
  # extract schoenfeld residuals and plot
  sch_res <- residuals(x, type = "schoenfeld") %>% 
    as.data.frame() %>% 
    mutate(times = round(as.numeric(sub("X", "", row.names(.))))) %>% 
    pivot_longer(-times, values_to = "residual") %>% 
    left_join(individual_ps) %>% 
    mutate(lab = ifelse(name == "interday_dist_z", "Beta(t) for daily distance",
                              ifelse(name == "max_height_z", "Beta(t) for max. daily height",
                                     ifelse(name == "med_circ_rad", "Beta(t) for circling radius",
                                            ifelse(name == "med_vspeed", "Beta(t) for vertical speed",
                                                   ifelse(name == "med_ws", "Beta(t) for wind speed",
                                                          NA))))),
           lab = paste0(lab, "\nSchoenfeld Individual Test p: ", round(p_value, 5)))
  
  
  schoenfeld_res <- ggplot(sch_res, aes(x = times, y = residual)) +
    geom_point() +
    geom_smooth(method = "loess", color = col) +
    labs(x = "Event time", y = "Schoenfeld residual") +
    facet_wrap(~lab, nrow = length(unique(sch_res$name)))

  return(schoenfeld_res)
})


# png(filename = "/home/hbronnvik/Documents/chapter4/independence_25/surv_interval_residuals.png",
#     width = 8.3, height = 11, units = "in", res = 300)
ggarrange(res_plts[[3]]+labs(title = "Soaring performance"), 
          res_plts[[4]]+labs(title = ""),
          res_plts[[1]]+labs(title = "Daily movement"), 
          res_plts[[2]]+labs(title = ""),
          labels = "AUTO")
# dev.off()

# facet labels
fac_labs <- c("A) Days until first excursion", 
              "B) Excursion distance (km)", 
              "C) Excursion duration (hours)",
              "D) Days until emigration")
names(fac_labs) <- c("days_to_leave", "excursion_dist", "excursion_duration", "pfdp")
# png(filename = paste0("/home/hbronnvik/Documents/chapter4/figures/independence/ex_descriptions.png"),
#     height = 8.3, width = 11, units = "in", res = 300)
first_ex %>% 
  dplyr::select(excursion_duration, excursion_dist, days_to_leave) %>% 
  pivot_longer(cols = 1:3) %>% 
  rbind(dates %>% 
          filter(did_disperse == T) %>% 
          dplyr::select(pfdp) %>% 
          rename(value = pfdp) %>% 
          mutate(name = "pfdp")) %>% 
  ggplot(aes(value)) +
  geom_histogram(color = "black", bins = 30) +
  scale_y_continuous(expand = c(0, 0), n.breaks = 7) +
  labs(x = "", y = "Count") +
  facet_wrap(~name, scales = "free", labeller = labeller(name = fac_labs))
# dev.off()

# GPS data with bursts removed
polys <- lapply(list.files("/home/hbronnvik/Documents/chapter4/polygons", full.names = T), readRDS) %>% reduce(rbind)

library(OpenStreetMap)
rel_poly <- polys %>% 
  filter(id == 1201987685)
rel_poly <- rel_poly %>% 
  rbind(rel_poly[1,])

ind <- list.files("/home/hbronnvik/Documents/GE_data/no_burst_in_out", 
                  full.names = T, pattern = "1201987685") %>% 
  readRDS() %>% 
  filter(between(location.lat, min(rel_poly$y)-0.1, max(rel_poly$y)+0.1) &
           between(location.long, min(rel_poly$x)-0.1, max(rel_poly$x)+0.1)) %>% 
  # add on emigration timing
  mutate(emigration = dates$dispersal_date[dates$individual.id == unique(individual.id)],
         relative_day = round(as.numeric(difftime(date(timestamp), emigration, units = "days")))) %>% 
  filter(between(relative_day, -200, 201))

sa_map <- openmap(c(max(rel_poly$y)+0.1, max(rel_poly$x)+0.1), 
                  c(min(rel_poly$y)-0.1, min(rel_poly$x)-0.1), 
                  zoom = 10, type = "esri-terrain", mergeTiles = T)
# reproject onto WGS84
sa_map <- openproj(sa_map)
# make the map

# png(filename = paste0("/home/hbronnvik/Documents/chapter4/figures/independence/example_poly_Stürfis.png"),
#     height = 8.3, width = 11, units = "in", res = 300)
autoplot.OpenStreetMap(sa_map) +
  geom_path(data = rel_poly, aes(x, y), color = "black", linewidth = 1) +
  geom_point(data = rel_poly, aes(x, y), color = "black", size = 0.75) +
  geom_point(data = ind, aes(location.long, location.lat, color = relative_day)) +
  scale_color_gradientn("Days to/since\nemigration", colors = rev(colfunc(200))) +
  # coord_sf(xlim = c(min(rel_poly$x)-0.1, max(rel_poly$x)+0.1),
  #          ylim = c(min(rel_poly$y)-0.1, max(rel_poly$y)+0.1)) +
  theme_void()
# dev.off()

daily_data <- lapply(c("max_height", "interday_dist"), function(p){
  pane_lab <- ifelse(p == "max_height", "Maximum height\nabove ellipsoid (m)",
                     ifelse(p == "interday_dist", "log distance\nbetween days (m)",
                                   p))
  pl <- coarse_intervals %>%
    mutate(interday_dist = log(interday_dist)) %>% 
    dplyr::select(max_height, interday_dist, dst, date) %>% 
    pivot_longer(cols = c(max_height, interday_dist)) %>% 
    mutate(month = month(date),
           month = factor(month, levels = c(6:12, 1:5))) %>% 
    filter(name == p) %>% 
    ggplot(aes(dst, value, color = month)) +
    geom_point() +
    geom_smooth(color = "grey50") +
    scale_x_continuous(breaks = seq(0, 300, 20), limits = c(0, 280)) +
    scale_color_manual("Month", values = c("#90272D", "#A42D32", "#C95A53", "#F4C1AB", 
                                           "#FDEDE4", "#DDECF4", "#9AC5DE", "#4F91BC",
                                           "#2B69A4", "#174A88", "#0A3162")) +
    labs(x = "Days since tagging", y = pane_lab)
  return(pl)
})
dd <- ggpubr::ggarrange(daily_data[[1]], daily_data[[2]], nrow = 2, labels = c("D", "E"), 
                  common.legend = T, legend = "right")

# visualize (as above with the daily-scale patterns)
pls <- lapply(c("med_circ_rad", "med_vspeed", "med_ws"), function(p){
  pane_lab <- ifelse(p == "med_circ_rad", "Circling\nradius (m)",
                     ifelse(p == "med_vspeed", "Vertical\nspeed (m/s)",
                            ifelse(p == "med_ws", "Wind\nspeed (m/s)", 
                                   p)))
  pl <- soar_intervals %>%
    dplyr::select(med_circ_rad, med_vspeed, med_ws, dst, date) %>% 
    pivot_longer(cols = med_circ_rad:med_ws) %>% 
    mutate(month = month(date),
           month = factor(month, levels = c(6:12, 1:5))) %>% 
    filter(name == p) %>% 
    ggplot(aes(dst, value, color = month)) +
    geom_point() +
    geom_smooth(color = "grey50") +
    scale_x_continuous(breaks = seq(0, 300, 20), limits = c(0, 280)) +
    scale_color_manual("Month", values = c("#C95A53", "#F4C1AB", "#FDEDE4", 
                                           "#DDECF4", "#9AC5DE", "#4F91BC",
                                           "#2B69A4", "#174A88", "#0A3162")) +
    labs(x = "Days since tagging", y = pane_lab) +
    theme(legend.position = "none")
  return(pl)
})
fd <- ggpubr::ggarrange(pls[[1]], pls[[2]], pls[[3]], nrow = 3, labels = "AUTO")
# png(filename = paste0("/home/hbronnvik/Documents/chapter4/independence_25/input_daily_interval_data.png"),
#     height = 8.3, width = 11, units = "in", res = 300)
ggpubr::ggarrange(fd, dd, ncol = 2)
# dev.off()

# plot all metrics together
# png(filename = paste0("/home/hbronnvik/Documents/chapter4/figures/independence/both_total.png"),
#     height = 8.3, width = 11, units = "in", res = 300)
ggpubr::ggarrange(ggpubr::ggarrange(pls[[1]], pls[[2]], pls[[3]], nrow = 3, 
                                    labels = c("A", "C", "E")),
                  ggpubr::ggarrange(idd, mgs2, mgs, nrow = 3, 
                                    labels = c("B", "D", "F"), 
                                    common.legend = T, legend = "right"),
                  nrow = 1, ncol = 2)
# dev.off()