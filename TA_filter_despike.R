# Script to experiment despike() on TA data because of spike for acid wash around midnight

if (Sys.getenv("LOGNAME") == "samir") path = "../../pCloud\ Sync/exp168_AWIPEV-CO2/fb_awipev-co2_server/"

load(file = paste0(path,"ny-alesund/data/NRT_data/all_nydata_hour.Rdata"))
d_hour <- d_hour %>%
  dplyr::select(datetime, AT_filtered)

# despike
d_hour <- d_hour %>%
  dplyr::mutate(AT_filtered2= despike(d_hour$AT_filtered, reference= "median", n=0.3, k=217, replace="NA"))%>%
#dplyr::mutate(AT_filtered2= despike(d_hour$AT_filtered, reference= "median", n=0.3, k=217, replace="NA"))
dplyr::mutate(AT_filtered3= despike(d_hour$AT_filtered2, reference= "median", n=0.3, k=217, replace="NA"))

# plot
at_fb_xts <- dplyr::select(d_hour, datetime,AT_filtered,AT_filtered2,AT_filtered3)
at_fb_xts <- as.xts(at_fb_xts, order.by = at_fb_xts$datetime)
dygraph(at_fb_xts, ylab="at") %>%
  dySeries("AT_filtered", color = "blue", strokeWidth = 0, label = "raw") %>%
  dySeries("AT_filtered2", color = "red", strokeWidth = 0, label = "filt") %>%
  dySeries("AT_filtered3", color = "green", strokeWidth = 0, label = "filt3") %>%
  dyHighlight(highlightCircleSize = 8,highlightSeriesBackgroundAlpha = 0.2,hideOnMouseOut = TRUE) %>%
  dyOptions(drawGrid = TRUE, drawPoints = TRUE, pointSize = 2,useDataTimezone = TRUE) %>%
  dyRangeSelector(height = 30)

tmp <- NULL
df <- expand.grid(n = seq(0.1, 1, by = 0.1), k =as.integer(seq(25, 221, by = 4)))
df$N <- NULL
for (i in 1:nrow(df)) {
  tmp <- despike(d_hour$AT_filtered, reference= "median", n=df$n[i], k=df$k[i], replace="NA")
  df$N[i] <- length(tmp[!is.na(tmp)])
  df$V[i] <- sd(tmp, na.rm = TRUE)
  }
summary(df)
#We want to keep data with highest N AND smallest Variance
df$percentN <- (df$N/max(df$N))*100
acceptable <- df%>%
  mutate(N = ifelse(percentN > 70 & V <80, N, NA),
         k = ifelse(percentN > 70 & V <80, k, NA))%>%
  filter(!is.na(k))

ggplot(data=df, aes(x = k, y = n)) +
  geom_point(aes(colour = N))
  
ggplot(data=df, aes(x = k, y = n)) +
  geom_point(aes(colour = V)) 

# plot
at_fb_xts <- dplyr::select(d_hour, datetime, AT_filtered,AT_filtered2,AT_filtered3)
at_fb_xts <- as.xts(at_fb_xts, order.by = at_fb_xts$datetime)
dygraph(at_fb_xts, ylab="at") %>%
  dySeries("AT_filtered", color = "blue", strokeWidth = 0, label = "raw") %>%
  dySeries("AT_filtered2", color = "red", strokeWidth = 0, label = "filt") %>%
  dySeries("AT_filtered3", color = "green", strokeWidth = 0, label = "filt3") %>%
  dyHighlight(highlightCircleSize = 8,highlightSeriesBackgroundAlpha = 0.2,hideOnMouseOut = TRUE) %>%
  dyOptions(drawGrid = TRUE, drawPoints = TRUE, pointSize = 2,useDataTimezone = TRUE) %>%
  dyRangeSelector(height = 30)  
    