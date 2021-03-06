## Figures AWIPEV-CO2 for ESSD paper

# set up ####
Sys.setlocale("LC_ALL", "en_US.UTF-8")
Sys.setenv(TZ='UTC') # on utilise UTC
rm(list = ls())
if (!require("tidyverse")) install.packages("tidyverse")
library(tidyverse)
if (!require("seacarb")) install.packages("seacarb")
library("seacarb")
if (!require("gridExtra")) install.packages("gridExtra")
library("gridExtra")
if (!require("reshape2")) install.packages("reshape2")
library("reshape2")
if (!require("lubridate")) install.packages("lubridate")
library("lubridate")
if (!require("lmtest")) install.packages("lmtest")
library("lmtest")
if (!require("grid")) install.packages("grid")
library(grid)
if (!require("viridis")) install.packages("viridis")
library("viridis")
if (!require("lmodel2")) install.packages("lmodel2")
library("lmodel2")
if (!require("xts")) install.packages("xts")
library("xts")
if (!require("scales")) install.packages("scales")
library("scales")
if (!require("cowplot")) install.packages("cowplot")
library(cowplot)
if (!require("htmlwidgets")) install.packages("htmlwidgets")
library("htmlwidgets")
if (!require("hms")) install.packages("hms")
library("hms")
if (!require("ragg")) install.packages("ragg")
library("ragg")

#################### define who is the user and define path
if (Sys.getenv("LOGNAME") == "gattuso") path = "../../pCloud\ Sync/Documents/experiments/exp168_awipev-CO2/"
if (Sys.getenv("LOGNAME") == "samir") path = "../../pCloud\ Sync/exp168_awipev-CO2/"

#################### Function RMSE root mean square error
rmse <- function(error) { sqrt(mean(error^2, na.rm=TRUE)) }
#rmse(fit$residuals)

#################### Function RMSE for model 2 regression root mean square error
lm2 <- function(x, y) { 
  fit <- lmodel2(y  ~ x , nperm = 99)
  intercept <- fit$regression.results[2,2]
  slope <- fit$regression.results[2,3]
  predicted <- intercept + x * slope
  error <- predicted - x
  my_list <- list(fit = fit, rmse = sqrt(mean(error^2, na.rm=TRUE)))
  return(my_list)
  }

mytheme <- function(size_labs = 6, face_font="plain", ...) {
  theme_bw() +
    theme(axis.text.x = element_blank(),
          #axis.text.x = element_text(face=face_font, size=size_labs, color="black"),
          axis.title.x = element_text(face=face_font, size=size_labs, margin = margin(0,0,0,0,"pt")),
          axis.text.y = element_text(face=face_font, color="black", size=size_labs),
          axis.title.y = element_text(face=face_font, size=size_labs),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(size=0.1),
          title = element_text(face=face_font, size=size_labs),
          axis.ticks.length = unit(1.1, "mm"),
          panel.grid.major = element_line(size = 0.05, color="grey50", linetype="dashed"),
          panel.grid.minor = element_blank(),
          #panel.grid.minor = element_line(size = rel(0.5), color = "grey50", linetype="dotted"),
          aspect.ratio = 1 / 3,
          plot.margin = margin(t = -0.3, r = -0.3, b = 0, l = -0.3, unit = "cm")
    )
}

# Read data ####
d <-
  readRDS(file = paste0(
    path = path,
    file = "fb_data/pangaea.rds"
  ))

# Read Zeppelin atmospheric CO2 data ####
# Atmospheric CO2 Zeppelin
# downloaded 2020-08-19 from https://gaw.kishou.go.jp/search/file/0054-6001-1001-01-01-9999
# value:units : ppm
# value:long_name : measured_mole_fraction_of_trace_gas_in_dry_air
atm_co2 <-
  read_delim(file = paste0(path, "atmospheric_co2_zeppelin/co2_zep_surface-insitu_54_9999-9999_hourly.txt"),
             delim = " ",
             comment = "#",
             #skip = 182,
             col_names = TRUE,
             na = "-999.999",
  ) %>%
  dplyr::mutate(co2 = value,
                datetime = as_datetime(paste(paste(year, month, day, sep="-"), paste(hour, minute, second, sep="-"), sep = " "))
  ) %>%
  dplyr::select(datetime, co2)

d <- dplyr::left_join(d, atm_co2, by = "datetime")

# time series ####
xlim <- c(as.POSIXct("2015-07-01 00:00:00"), as.POSIXct("2020-07-01 00:00:00"))
# calculate omeaga arag
z <- d %>% 
  dplyr::filter(!(is.na(pco2) | is.na(at_calc) | is.na(s_mix) |is.na(t_11m) | is.na(pressure)))
oa_tbf <- as_tibble(
  carb(24, z$pco2, z$at_calc*1e-6, S = z$s_mix, T = z$t_11m, 
       P = z$pressure/10, Pt = 0, Sit = 0, k1k2 = "l", kf = "dg", 
       ks = "d", pHscale = "T", b = "u74", warn = "n")
) %>% 
  dplyr::mutate(datetime = z$datetime,
                oa = OmegaAragonite) %>% 
  dplyr::select(datetime, oa)
# select time series data and add Omega
d <- as_tibble(d) #%>%
  #dplyr::select(datetime, pressure, t_11m, s_mix, pco2, pH_sf, 
  #              ph_dur, pco2_calc, at_calc)
d <- left_join(d, oa_tbf, by = 'datetime')

# time series temp
ts_temp <- d %>%
  dplyr::filter(!is.na(t_11m)) %>% 
  ggplot() +
  geom_point(aes(x = datetime, y = t_11m), col = "blue", size = 0.05, na.rm = TRUE) +
  labs(title = "", x = "", y = "Temp. 11 m (°C)") +
  mytheme()
# time series sal
ts_sal <- d %>%
  dplyr::filter(!is.na(s_mix)) %>% 
  ggplot() +
  geom_point(aes(x = datetime, y = s_mix), col = "blue", size = 0.05, na.rm = TRUE) + 
  labs(title = "", x = "", y = "Salinity") +
  mytheme()
# time series pco2
ts_co2 <- d %>%
  #dplyr::filter(!is.na(pco2)) %>% 
  ggplot() +
  geom_point(aes(x = datetime, y = pco2), col = "blue", size = 0.05, na.rm = TRUE) +
  geom_point(aes(x = datetime, y = co2), col = "red", size = 0.05, na.rm = TRUE) +
  scale_x_datetime(limits = xlim) +
  labs(title = "", x = "", y = expression(paste("pC", O[2], " (", mu, "atm)"))) +
  mytheme()
# time series pH
ts_ph <- d %>%
  dplyr::select(datetime, pH_sf, ph_dur) %>%
  dplyr::filter(!is.na(pH_sf)) %>%
  pivot_longer(-datetime, names_to = "pH", values_to = "value") %>%
  ggplot() +
    geom_point(aes(x = datetime, y = value, colour = pH), size = 0.05, na.rm = TRUE) +
  scale_x_datetime(limits = xlim) +
  labs(title = "",
       x = "",
       y = expression(paste("p", H[T]))) +
  # annotate(
  #   geom = "text",
  #   x = as.POSIXct("2018-06-01"),
  #   y = 8.35,
  #   label = "seaFET",
  #   color = "blue", size = 2
  # ) +
  # annotate(
  #   geom = "text",
  #   x = as.POSIXct("2018-06-01"),
  #   y = 8.25,
  #   label = "Durafet",
  #   color = "red", size = 2
  # ) +
  mytheme() +
  theme(legend.position = "none")
# time series AT
ts_at <- d %>%
  dplyr::filter(!is.na(at_calc)) %>% 
  ggplot() +
  geom_point(aes(x = datetime, y = at_calc), col = "blue", size = 0.05, na.rm = TRUE) +
  scale_x_datetime(limits = xlim) +
  labs(title = "", x = "", y = expression(paste(italic(A)[T]," ","(",mu, mol, " ", kg^-1,")"))) +
  mytheme()
# time series Omega
ts_oa <- d %>%
  dplyr::filter(!is.na(oa)) %>% 
  ggplot() +
  geom_point(aes(x = datetime, y = oa), col = "blue", size = 0.05, na.rm = TRUE) +
  scale_x_datetime(limits = xlim) +
  labs(title = "", x = "Time", y = expression(paste(Omega[a]))) +
  mytheme()
  # theme(axis.title.x = element_text(face=face_font, size=size_labs, color="black"),
  #       axis.text.x = element_text(face=face_font, size=size_labs, color="black"))

#assemble time-series figure
cp <- cowplot::plot_grid(ts_sal, ts_temp, ts_co2, ts_ph, ts_at, ts_oa,
                         align = "v",
                         ncol = 1,
                         axis = "tblr"
                         #labels = "auto",
                         #label_size = 7,
                         #label_x = 0.9,
                         #label_y = 0.95
                         #width = 18,
                         #units = "cm"
                         )
cowplot::ggsave2(filename = "figures/essd/ts_gg.png", plot = cp, width = 3, height = 6, unit = "cm")
cowplot::save_plot(filename = "figures/essd/ts_gg2.png", plot = cp, base_width = 3, base_height = 6)

# Boxplots ####
mytheme_bp <- function(size_labs = 9, face_font="plain", ...) {
  mytheme() +
    theme(aspect.ratio = 1/2,
          plot.margin = margin(t = -2, r = 0.2, b = -2, l = 0, unit = "cm")
#          axis.text.x=element_blank()
)
}

d_long <- d %>% 
  dplyr::mutate(month = as.factor(lubridate::month(x = datetime, label = TRUE, abbr = TRUE))) %>% 
  tidyr::pivot_longer(-c(datetime, month), names_to = "variable", values_to = "value") %>%
  dplyr::filter(!is.na(value)) # remove NAs to avoid warnings

s_mix_bp <- ggplot(filter(d_long, variable=="s_mix"), aes(x = month, y = value, group=month)) +
  geom_violin(fill='lightblue', alpha=0.5, size=0.2, trim = TRUE) +
  geom_boxplot(notch=TRUE, size=0.2, fill="grey50", outlier.color = "blue", outlier.size = 0.3) +
  scale_x_discrete(breaks=c("Jan", "", "Mar", "", "May", "", "Jul", "", "Sep", "", "Nov", "")) +
  labs(title=NULL,x=NULL,y = "Salinity") +
  mytheme_bp() + 
  theme(axis.text.x=element_blank())
temp_11m_bp <- ggplot(filter(d_long, variable=="t_11m"), aes(x = month, y = value, group=month)) +
  geom_violin(fill='lightblue', alpha=0.5, size=0.2, trim = TRUE) +
  geom_boxplot(notch=TRUE, size=0.2, fill="grey50", outlier.color = "blue", outlier.size = 0.3) +
  scale_x_discrete(breaks=c("Jan", "", "Mar", "", "May", "", "Jul", "", "Sep", "", "Nov", "")) +
  labs(title=NULL,x=NULL,y = "Temp. 11 m (°C)") +
  mytheme_bp() + 
  theme(axis.text.x=element_blank())
pco2_bp <- ggplot(filter(d_long, variable=="pco2"), aes(x = month, y = value, group=month)) +
  geom_violin(fill='lightblue', alpha=0.5, size=0.2, trim = TRUE) +
  geom_boxplot(notch=TRUE, size=0.2, fill="grey50", outlier.color = "blue", outlier.size = 0.3) +
  scale_x_discrete(breaks=c("Jan", "", "Mar", "", "May", "", "Jul", "", "Sep", "", "Nov", "")) +
  labs(title=NULL,x=NULL,y = expression(paste("pC", O[2], " (", mu, "atm)"))) +
  mytheme_bp() + 
  theme(axis.text.x=element_blank())
ph_bp <- ggplot(filter(d_long, variable=="pH_sf"), aes(x = month, y = value, group=month)) +
  geom_violin(fill='lightblue', alpha=0.5, size=0.2, trim = TRUE) +
  geom_boxplot(notch=TRUE, size=0.2, fill="grey50", outlier.color = "blue", outlier.size = 0.3) +
  scale_x_discrete(breaks=c("Jan", "", "Mar", "", "May", "", "Jul", "", "Sep", "", "Nov", "")) +
  labs(title=NULL,x=NULL,y = expression(paste("p", H[T]))) +
  mytheme_bp() + 
  theme(axis.text.x=element_blank())
at_bp <- ggplot(filter(d_long, variable=="at_calc"), aes(x = month, y = value, group=month)) +
  geom_violin(fill='lightblue', alpha=0.5, size=0.2, trim = TRUE) +
  geom_boxplot(notch=TRUE, size=0.2, fill="grey50", outlier.color = "blue", outlier.size = 0.3) +
  scale_x_discrete(breaks=c("Jan", "", "Mar", "", "May", "", "Jul", "", "Sep", "", "Nov", "")) +
  labs(title=NULL,x=NULL,y = expression(paste(italic(A)[T]," ","(",mu, mol, " ", kg^-1,")"))) +
  mytheme_bp() +
  theme(axis.text.x = element_text(vjust = 1))
oa_bp <- ggplot(filter(d_long, variable=="oa"), aes(x = month, y = value, group=month)) +
  geom_violin(fill='lightblue', alpha=0.5, size=0.2, trim = TRUE) +
  geom_boxplot(notch=TRUE, size=0.2, fill="grey50", outlier.color = "blue", outlier.size = 0.3) +
  scale_x_discrete(breaks=c("Jan", "", "Mar", "", "May", "", "Jul", "", "Sep", "", "Nov", "")) +
  labs(title=NULL,x=NULL,y = expression(paste(Omega[a]))) +
  mytheme_bp() +
  theme(axis.text.x = element_text(vjust = 1))
g <- cowplot::plot_grid(s_mix_bp, temp_11m_bp, pco2_bp, ph_bp, at_bp, oa_bp, ncol=2,
                        align="v",
                        labels="AUTO")
ggsave(file="figures/essd/boxplots.png", g,  width = 18, height = 14, units = "cm")

# pCO2 ####
da <- dplyr::filter(d, !is.na(pco2_calc), !is.na(pco2))
lm2 <- lm2(x = da$pco2, y = da$pco2_calc)
#fit <- lmodel2(data = d,  pco2 ~ pco2_calc , nperm = 99)
p <-  ggplot(d, aes(x=pco2_calc, y= pco2), na.rm = TRUE) +
  geom_point(color = "blue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_abline(aes(intercept = lm2$fit$regression.results[2,2], slope = lm2$fit$regression.results[2,3]), colour = "blue") +
  labs(x = expression(paste("Calculated pC", O[2], " (", mu, "atm)")), 
       y = expression(paste("In situ pC", O[2], " (", mu, "atm)")),
       title = paste("Adj R2 = ", signif(lm2$fit$rsquare, 3),
                     "; Intercept =", signif(lm2$fit$regression.results[2,2], 3),
                     "; Slope =", signif(lm2$fit$regression.results[2,3], 3),
                     "; P =", signif(lm2$fit$P.param, 3),
                     "; RMSE = ", signif(lm2$rmse, 3))) +
  coord_fixed(ratio = 1 ,xlim=c(200, 400) , ylim=c(200, 400))+
  mytheme(size_labs = 10, plot.title = element_text(face=face_font, size=size_labs, color="black")) +
  theme(aspect.ratio = 1,
        plot.margin = margin(t = 0.5, unit = "cm"))
p
ggsave(file="figures/essd/pco2.png", p,  width = 14, height = 14, units = "cm")


# pH: seafet vs spectro ####
da <- dplyr::mutate(ph_calc = 
                      dplyr::filter(d, !is.na(ph_sf), !is.na(ph_s_sf_t_insi))
lm2 <- lm2(x = da$pco2_calc, y = da$pco2)
geom_point(color = "blue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_abline(aes(intercept = lm2$fit$regression.results[2,2], slope = lm2$fit$regression.results[2,3]), colour = "blue") +
  labs(x = expression(paste("Calculated pC", O[2], " (", mu, "atm)")), 
       y = expression(paste("In situ pC", O[2], " (", mu, "atm)")),
       title = paste("Adj R2 = ", signif(lm2$fit$rsquare, 3),
                     "; Intercept =", signif(lm2$fit$regression.results[2,2], 3),
                     "; Slope =", signif(lm2$fit$regression.results[2,3], 3),
                     "; P =", signif(lm2$fit$P.param, 3),
                     "; RMSE = ", signif(lm2$rmse, 3))) +
  coord_fixed(ratio = 1 ,xlim=c(200, 400) , ylim=c(200, 400))+
  mytheme(size_labs = 10, plot.title = element_text(face=face_font, size=size_labs, color="black")) +
  theme(aspect.ratio = 1,
        plot.margin = margin(t = 0.5, unit = "cm"))
p
ggsave(file="figures/essd/dur_spec.png", p,  width = 9, height = 9, units = "cm")


