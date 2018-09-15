
######################################################
#
# Included is data on daily DUIs as recorded by the
# Seattle and San Francisco city police departments
# from July 2010 to Jun 2014.
#
# Uber operated in SF over the entire period and
# launched in Seattle on eventdate 391 (08/11/2011) 
# in the data set. Data are ordered by day for each city.
#
# variables
# ---------
# eventdate : starts at 0; increases by 1 each day
# incidents : number of daily DUIs recorded
# post : dummy variable equal to 1 if date after Uber Seattle launch
# sea : dummy equal to 1 if seattle 
# sf : dummy equal to 1 if San Francisco
# dayofweek : Sun-Sat, factor
# marijuana : marijuana legal (eventdate 876 ())
#

require(tidyverse)
if (!require(lubridate)) {
  install.packages("lubridate"); library(lubridate)
}

df_raw <- read_csv("data/uber_dui.csv", na=c("#N/A"))
summary(df_raw) # 493 days had no incident data;
# For SFPD, only days w/ incidents reported;
# Not sure about Seattle PD

# To make things simpler, let's just assume all missing
# days had no DUI arrests in either city
df <- df_raw %>% 
  transmute(eventdate,
            date=as_date("2010-7-16") + eventdate,
            dayofweek,
            city=ifelse(sea==1, "Seattle", "SF"),
            incidents=ifelse(is.na(incidents), 0, incidents),
            post,
            marijuana
         )

summary(df)

## ------------------------------------------------------------------------
library(ggplot2)

ggplot(df, aes(x=date, y=incidents, color=city)) + 
  geom_line() +
  geom_vline(xintercept = as_date("2011-08-11")) +
  geom_vline(xintercept = as_date("2012-12-08"))

## ------------------------------------------------------------------------
df_sea <- df %>% filter(city=="Seattle")

if (!require(ggfortify)) {
  install.packages("ggfortify"); library(ggfortify)
}

# Visually inspect the weekly pattern
sea_ts <- ts(df_sea$incidents, start=2010, frequency = 7)
autoplot(stl(sea_ts, s.window = 'periodic'), ts.colour = 'blue')

# Visually inspect the monthly pattern
sea_ts <- ts(df_sea$incidents, start=2010, frequency = 30)
autoplot(stl(sea_ts, s.window = 'periodic'), ts.colour = 'blue')

if (!require(strucchange)) {
  install.packages("strucchange"); library(strucchange)
}

autoplot(breakpoints(incidents ~ eventdate, data=df_sea))

## ------------------------------------------------------------------------
# Base Seattle model
(mod1_sea <- lm(incidents ~ dayofweek + eventdate, data=df_sea)) %>% 
  summary()

# Dummy variable approach
(mod2_sea <- lm(incidents ~ dayofweek + eventdate + post, data=df_sea)) %>% 
  summary()

anova(mod1_sea, mod2_sea)

# Chow test on whether there is a structure change on day 392 (8/11/2011)
# when Uber entered Seattle

# NOTE: Data must be in time order!
df_sea <- df_sea %>% arrange(eventdate)
sctest(incidents ~ dayofweek + eventdate, type = "Chow", point = 392, data = df_sea)

## ------------------------------------------------------------------------
# First, let's establish the base model that
# includes both cities
(mod3 <- lm(incidents ~ post + dayofweek + eventdate * city, data=df)) %>% 
  summary()

# Specify the DID form of the base model


# Add a control for marijuana legalization



