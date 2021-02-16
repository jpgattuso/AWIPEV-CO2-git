# Extract data from d_hour
Sys.setlocale("LC_ALL", "en_US.UTF-8")
Sys.setenv(TZ='UTC') # on utilise UTC
rm(list = ls())
if (!require("tidyverse")) install.packages("tidyverse")
library(tidyverse)

#define who is the user and define path
if (Sys.getenv("LOGNAME") == "gattuso") path = "../../pCloud\ Sync/Documents/experiments/exp168_awipev-CO2/"
if (Sys.getenv("LOGNAME") == "samir") path = "../../pCloud\ Sync/exp168_awipev-CO2/"

d <- readRDS(file= paste0(path, "fb_data/d_all.rds")) %>% 
  dplyr::select(datetime, temp_insitu_11m, sal_fb) %>% 
  dplyr::mutate(week = lubridate::isoweek(datetime)) %>% 
  dplyr::filter(week >= 22 & week <= 35 )

weekly_means <- d %>% 
  dplyr::group_by(week) %>% 
  dplyr::summarize(mean_temp_insitu_11m = mean(temp_insitu_11m, na.rm = TRUE),
                   sd_temp_insitu_11m = sd(temp_insitu_11m, na.rm = TRUE),
                   mean_sal_fb = mean(sal_fb, na.rm = TRUE),
                   sd_sal_fb = sd(sal_fb, na.rm = TRUE))
write.csv(file=paste0(path, "fb_data/sal_temp.csv"), weekly_means)
          