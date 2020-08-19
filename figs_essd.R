## Figures AWIPEV-CO2 for ESSD paper

# set up ####
Sys.setlocale("LC_ALL", "en_US.UTF-8")
Sys.setenv(TZ='UTC') # on utilise UTC
rm(list = ls())
if (!require("tidyverse")) install.packages("tidyverse")
library(tidyverse)
# if (!require("robfilter")) install.packages("robfilter")
# library("robfilter")
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
# if (!require("dygraphs")) install.packages("dygraphs")
# library("dygraphs")
# if (!require("knitr")) install.packages("knitr")
# require("knitr")
if (!require("lmodel2")) install.packages("lmodel2")
library("lmodel2")
# if (!require("captioner")) install.packages("captioner")
# library("captioner")
if (!require("xts")) install.packages("xts")
library("xts")
# if (!require("seismicRoll")) install.packages("seismicRoll")
# library("seismicRoll")
if (!require("scales")) install.packages("scales")
library("scales")
if (!require("cowplot")) install.packages("cowplot")
library(cowplot)
if (!require("htmlwidgets")) install.packages("htmlwidgets")
library("htmlwidgets")
if (!require("hms")) install.packages("hms")
library("hms")

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
    theme(axis.text.x = element_text(face=face_font, size=size_labs, color="black"),
          axis.title.x = element_text(face=face_font, size=size_labs, margin = margin(0,0,0,0,"pt")),
          axis.text.y = element_text(face=face_font, color="black", size=size_labs),
          axis.title.y = element_text(face=face_font, size=size_labs),
          axis.ticks.x = element_line(size=0.1),
          axis.ticks.y = element_line(size=0.1),
          title = element_text(face=face_font, size=size_labs),
          axis.ticks.length = unit(1.1, "mm"),
          panel.grid.major = element_line(size = 0.25, color="grey50", linetype="dashed"),
          panel.grid.minor = element_line(size = rel(0.5), color = "grey50", linetype="dotted"),
          aspect.ratio = 1 / 3,
          plot.margin = margin(t = 0, r = 1, b = 0, l = 1, unit = "lines")
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
             skip = 181,
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
  geom_point(aes(x = datetime, y = t_11m), col = "blue", size = 0.15, na.rm = TRUE) +
  labs(title = "", x = "", y = "Temp. 11 m (°C)") +
  mytheme(size_labs = 8)
# time series sal
ts_sal <- d %>%
  dplyr::filter(!is.na(s_mix)) %>% 
  ggplot() +
  geom_point(aes(x = datetime, y = s_mix), col = "blue", size = 0.15, na.rm = TRUE) + 
  labs(title = "", x = "", y = "Salinity") +
  mytheme(size_labs = 8)
# time series pco2
ts_co2 <- d %>%
  #dplyr::filter(!is.na(pco2)) %>% 
  ggplot() +
  geom_point(aes(x = datetime, y = pco2), col = "blue", size = 0.15, na.rm = TRUE) +
  geom_point(aes(x = datetime, y = co2), col = "red", size = 0.15, na.rm = TRUE) +
  scale_x_datetime(limits = xlim) +
  labs(title = "", x = "", y = expression(paste("pC", O[2], " (", mu, "atm)"))) +
  mytheme(size_labs = 8)
# time series pH
ts_ph <- d %>%
  dplyr::select(datetime, pH_sf, ph_dur) %>%
  dplyr::filter(!is.na(pH_sf)) %>%
  pivot_longer(-datetime, names_to = "pH", values_to = "value") %>%
  ggplot() +
  geom_point(aes(x = datetime, y = value, colour = pH), size = 0.15, na.rm = TRUE) +
  scale_x_datetime(limits = xlim) +
  labs(title = "",
       x = "",
       y = expression(paste("p", H[T]))) +
  # scale_color_manual(
  #   values = c("blue", "red"),
  #   name = "",
  #   breaks = c("pHint_tot_sf", "ph_dur_corr"),
  #   labels = c("SeaFET", "Durafet")
  # ) +
  annotate(
    geom = "text",
    x = as.POSIXct("2018-06-01"),
    y = 8.35,
    label = "seaFET",
    color = "blue", size = 2
  ) +
  annotate(
    geom = "text",
    x = as.POSIXct("2018-06-01"),
    y = 8.25,
    label = "Durafet",
    color = "red", size = 2
  ) +
  mytheme(size_labs = 8) +
  theme(legend.position = "none")
# time series AT
ts_at <- d %>%
  dplyr::filter(!is.na(at_calc)) %>% 
  ggplot() +
  geom_point(aes(x = datetime, y = at_calc), col = "blue", size = 0.15, na.rm = TRUE) +
  scale_x_datetime(limits = xlim) +
  labs(title = "", x = "", y = expression(paste(italic(A)[T]," ","(",mu, mol.kg^-1,")"))) +
  mytheme(size_labs = 8) +
  theme(axis.title.x = element_text(face="plain", size=8))
# time series Omega
ts_oa <- d %>%
  dplyr::filter(!is.na(oa)) %>% 
  ggplot() +
  geom_point(aes(x = datetime, y = oa), col = "blue", size = 0.15, na.rm = TRUE) +
  scale_x_datetime(limits = xlim) +
  labs(title = "", x = "", y = expression(paste(Omega[a]))) +
  mytheme(size_labs = 8) +
  theme(axis.title.x = element_text(face="plain", size=8))
#assemble time-series figure
cp <- cowplot::plot_grid(ts_sal, ts_temp, ts_co2, ts_ph, ts_at, ts_oa,
                         align = "v",
                         ncol = 1
                         #width = 18,
                         #units = "cm"
                         )
cowplot::ggsave2(filename = "figures/essd/ts_gg.png", plot = cp, height = 20, units = "cm")
grid.arrange(ts_sal, ts_temp, ts_pco2, ts_ph, ts_at, ts_oa, ncol = 1)
# ts_sal / ts_temp / ts_at
# grid.arrange(ts_sal, ts_temp, ts_pco2, ts_ph, ts_at, ncol = 1)
# cowplot::save_plot(filename = "figures/ts_sp.pdf", plot = cp)
# multiplot(ts_sal, ts_temp, ts_pco2, ts_ph, ts_at, cols = 1)

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
lm2 <- lm2(x = da$pco2_calc, y = da$pco2)
#fit <- lmodel2(data = d,  pco2 ~ pco2_calc , nperm = 99)
p <-  ggplot(d, aes(x=pco2_calc, y= pco2), na.rm = TRUE) +
  geom_point(color = "blue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_abline(aes(intercept = fit$regression.results[2,2], slope = fit$regression.results[2,3]), colour = "blue") +
  labs(x = expression(paste("Calculated pC", O[2], " (", mu, "atm)")), 
       y = expression(paste("In situ pC", O[2], " (", mu, "atm)")),
       title = paste("Adj R2 = ", signif(fit$rsquare, 3),
                     "; Intercept =", signif(fit$regression.results[2,2], 3),
                     "; Slope =", signif( fit$regression.results[2,3], 3),
                     "; P =", signif(fit$P.param, 3),
                     "; RMSE = ", signif(rmse, 3))) +
  coord_fixed(ratio = 1 ,xlim=c(200, 400) , ylim=c(200, 400))+
#  labs(x = expression(paste("Calculated pC", O[2], " (", mu, "atm)")), 
#       y = expression(paste("In situ pC", O[2], " (", mu, "atm)"))) + 
  mytheme(size_labs = 10) +
  theme(aspect.ratio=1, plot.title = element_text(size=7))
p
ggsave(file="figures/essd/pco2.png", p,  width = 14, height = 14, units = "cm")


# pH: durafet vs spectro ####
da <- dplyr::filter(d, !is.na(ph_dur), !is.na(ph_s_dur_t_fb))
rmse <- rmse2(x = da$ph_s_dur_t_fb, y = da$ph_dur)
p <-  ggplot(da, aes(x = ph_s_dur_t_fb, y = ph_dur)) +
  geom_point(color = "blue", na.rm = TRUE) + 
  scale_color_discrete(guide = "none")+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_abline(aes(intercept = fit$regression.results[2,2], 
                  slope = fit$regression.results[2,3]), 
              colour = "blue") +
  labs(x="Spectrophotometric pH", y="Durafet pH",
       title = paste("Adj R2 = ", signif(fit$rsquare, 3),
                     "; Intercept =", signif(fit$regression.results[2,2], 3),
                     "; Slope =", signif( fit$regression.results[2,3], 3),
                     "; P =", signif(fit$P.param, 3),
                     "; RMSE = ", signif(rmse, 3))) +
  coord_fixed(ratio = 1) +
  mytheme(size_labs = 6) +
  theme(aspect.ratio = 1)
p
ggsave(file="figures/essd/dur_spec.png", p,  width = 9, height = 9, units = "cm")


# pH: SeaFET vs Durafet ####
d_no_na <- dplyr::filter(!is.na())
ph_sf_4 <- pHinsi(pH=d$ph_sf, ALK=d$at_calc, Tinsi=4, Tlab=d$t_11m, Pinsi=11/10, S=d$s_insi, Pt=0, Sit=0, 
                  k1k2 = "l", kf = "dg", ks = "d", pHscale = "T", b = "u74")
  d <- d %>%
  dplyr::mutate(ph_sf_4 = )
fit <- lmodel2(data = d,  ph_sf_4 ~ ph_dur_4, nperm = 99)
p <-  ggplot(d, aes(x = ph_dur, y = ph_s_dur_t_fb)) +
  geom_point(color = "blue", na.rm = TRUE) + 
  scale_color_discrete(guide = "none")+
  geom_abline(slope=1, intercept = 0, linetype = "dashed") +
  geom_abline(aes(intercept = fit$regression.results[2,2], 
                  slope = fit$regression.results[2,3]), 
              colour = "blue") +
  labs(x="Durafet pH", y="Spectrophotometric pH",
       title = paste("Spec pH vs Durafet pH\nAdj R2 = ", signif(fit$rsquare, 3),
                     "\nIntercept =", signif(fit$regression.results[2,2], 3),
                     "\nSlope =", signif( fit$regression.results[2,3], 3),
                     "\nP =", signif(fit$P.param, 3))) +
  coord_fixed(ratio = 1) +
  mytheme(size_labs = 6) +
  theme(aspect.ratio = 1)
ggsave(file="figures/essd/dur_spec.png", p,  width = 9, height = 9, units = "cm")
p



